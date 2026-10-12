// Imports
#include <iostream>
#include <algorithm>
#include <array>
#include <cstddef>
#include <json/value.h>
#include <stdexcept>
#include <string>

// Self-Imports
#include "o00_Globals.h"
#include "o02_Mesh.h"
#include "o02_MeshFVM.h"

///// Mesh Generation /////
template <size_t Dim> void MeshFVM<Dim>::calculateFaces(std::array<size_t, Dim> cVec, Json::Value refData, std::array<std::vector<double>, Dim>& nFaces){
    // Control
    std::string sAlgo = refData["algorithm"].asString();
    size_t iAxis = refData["axis"].asInt(), NSec = refData["N"].asInt(), cNode = cVec[iAxis];
    double x0 = refData["range"][0].asDouble(), x1 = refData["range"][1].asDouble(); double length = x1 - x0;

    // Face Positions
    if (sAlgo == "Bidirectional"){ // Bidirectional (A, xC)
        
        // Mesh Parameters
        double A = refData["strength"].asDouble(), xC = refData["centering"].asDouble() - x0; 
        if (A > 2.5 || A < -2.5){ throw std::invalid_argument("Check .json"); }
        if (xC <= 0 || xC >= length){ throw std::invalid_argument("Check .json"); }
        // A (+): Closer to center, A (-): Away from center
        
        // Generate Mesh
        for (int i = cNode; i < cNode+NSec+1; i++) {
            nFaces[iAxis][i] = x0 + (static_cast<double>(i)-cNode) * length / NSec + A * (xC - (static_cast<double>(i)-cNode) * length / NSec) * (1 - (static_cast<double>(i)-cNode)/NSec) * (static_cast<double>(i)-cNode) / NSec;
        }

    } else if (sAlgo == "PowerLaw"){ // Power-Law (Kappa)
        
        // Mesh Parameters
        double kappa = refData["kappa"].asDouble();
        if (kappa <= 0){ std::exit(0); }
        // kappa < 1: Shifts right, kappa > 1: Shifts left

        // Generate Mesh
        for (int i = cNode; i < cNode+NSec+1; i++) {
            nFaces[iAxis][i] = x0 + length * pow(((static_cast<double>(i)-cNode) / NSec), kappa);
        }

    } else if(sAlgo == "HyperSingle"){ // Hyperbolic Tangent (Single-Sided)
        
        // Mesh Parameters
        double delta = refData["delta"].asDouble(); double A{}, B{};
        if (delta == 0){ throw std::invalid_argument("Check .json"); }
        // delta: independent of sign (+/-), higher value pushes towards the left

        // Generate Mesh
        for (int i = cNode; i < cNode+NSec+1; i++){
            A = tanh(delta * ((static_cast<double>(i) - cNode) / NSec - 1)); B = tanh(delta);
            nFaces[iAxis][i] = x0 + length * (1 + A / B);
        }

    } else if(sAlgo == "HyperDouble"){ // Hyperbolic Tangent (Double-Sided)
        
        // Mesh Parameters
        double delta = refData["delta"].asDouble(); double A{}, B{};
        if (delta == 0){ throw std::invalid_argument("Check .json"); }
        // delta: independent of sign (+/-), higher value pushes towards the edges

        // Generate Mesh
        for (int i = cNode; i < cNode+NSec+1; i++){
            A = tanh(delta*((static_cast<double>(i) - cNode)/NSec - 0.5)); B = tanh(0.5 * delta);
            nFaces[iAxis][i] = x0 + 0.5 * length * (1 + A/B);
        }

    } else if (sAlgo == "Exponential") { // Exponential

        // Mesh Parameters
        double kappa = refData["kappa"].asDouble();
        if (kappa == 0){ throw  std::invalid_argument("Check .json"); }
        // kappa < 0: Shifts right, kappa > 0: Shifts left 
        
        // Generate Mesh
        for (int i = cNode; i < cNode+NSec+1; i++){
            nFaces[iAxis][i] = x0 + length * (exp(kappa * (static_cast<double>(i) - cNode) / NSec) - 1) / (exp(kappa) - 1);
        }

    } else { throw std::logic_error("Mesh algorithm not recognized.\n"); } 
}

template <size_t Dim> void MeshFVM<Dim>::generateMeshSolver(MeshSolver<Dim>& Msh, Json::Value qNode, Json::Value sections, Json::Value refinement, Json::Value obstacles) {
	// Control (nD)
	for(Json::Value::ArrayIndex i = 0; i < Msh.N.size(); i++) {Msh.N[i] = qNode[i].asInt(); Msh.totNodes *= Msh.N[i];}

	// Geometry Resize (nD)
	for (size_t i = 0; i < Msh.N.size(); i++) {Msh.Faces[i].resize(Msh.N[i]+1); Msh.Nodes[i].resize(Msh.N[i]); Msh.deltaX[i].resize(Msh.N[i]); Msh.dX[i].resize(Msh.N[i]+1);}

    // Faces (nD)
    std::array<size_t, Dim> cNode{}; 
    for (Json::Value::ArrayIndex i = 0; i < refinement.size(); i++){
        // Calculate
        calculateFaces(cNode, refinement[i], Msh.Faces);

        // Control
        cNode[refinement[i]["axis"].asInt()] += refinement[i]["N"].asInt();
    }

    // Nodes (nD)
    for (size_t i = 0; i < Msh.N.size(); i++){
        for (size_t j = 0; j < Msh.N[i]; j++){
            Msh.Nodes[i][j] = 0.5 * (Msh.Faces[i][j] + Msh.Faces[i][j+1]);
        }
    }

    // Deltas (nD)
    for (size_t i = 0; i < Msh.N.size(); i++){
        for (size_t j = 0; j < Msh.deltaX[i].size(); j++){
            // Delta X
            Msh.deltaX[i][j] = Msh.Faces[i][j+1] - Msh.Faces[i][j];

            // dX
            if (j == Msh.deltaX[i].size()-1){continue;}
            Msh.dX[i][j+1] = Msh.Nodes[i][j+1] - Msh.Nodes[i][j];
        }

        // dX
        Msh.dX[i].front() = Msh.Nodes[i].front() - Msh.Faces[i].front();
        Msh.dX[i].back() = Msh.Faces[i].back() - Msh.Nodes[i].back();
    }


    std::cout << "S Size: "; for (std::vector<double>& Sk : Msh.S) {std::cout << Sk.size() << " ";} std::cout << "\n";


    // Resize (nD)
    Msh.nMat.resize(Msh.totNodes); for (std::vector<double>& Sk : Msh.S) {Sk.resize(Msh.totNodes, 1);} Msh.Vp.resize(Msh.totNodes, 1); Msh.bObs.resize(Msh.totNodes, false);
    Msh.Phi.resize(Msh.totNodes); Msh.sPhi.resize(Msh.totNodes); Msh.oPhi.resize(Msh.totNodes);
    

    std::cout << "S Size: "; for (std::vector<double>& Sk : Msh.S) {std::cout << Sk.size() << " ";} std::cout << "\n";
    throw std::invalid_argument("Test exit");


    // Sections
    std::array<size_t, Dim> i0{}, i1{}; size_t valMat{}; double valSource{};
    for (Json::Value::ArrayIndex i = 0; i < sections.size(); i++) {
        // Positions
        for (Json::Value::ArrayIndex j = 0; j < Dim; j++) {
            i0[j] = std::find(Msh.Faces[j].begin(), Msh.Faces[j].end(), sections[i]["x0"][j].asDouble()) - Msh.Faces[j].begin();
            i1[j] = std::find(Msh.Faces[j].begin(), Msh.Faces[j].end(), sections[i]["x1"][j].asDouble()) - Msh.Faces[j].begin() - 1;
        }

        // Properties
        valMat = sections[i]["material"].asInt(); valSource = sections[i]["source"].asDouble();
        runLoopMesh<Dim>(Msh.N, [&](std::array<size_t, 3> iX, size_t Ny, size_t Nx) {
                size_t idx = calcIndex(iX[0], iX[1], Ny, iX[2], Nx); Msh.nMat[idx] = valMat; Msh.sPhi[idx] = valSource;
                for (size_t nD = 0; nD < Dim; nD++) {
                    Msh.Vp[idx] *= Msh.deltaX[nD][iX[nD]];
                    for (size_t nDD = 0; nDD < Dim; nDD++) { if (nD != nDD) { Msh.S[nD][idx] *= Msh.deltaX[nDD][iX[nDD]]; } }
                } }, i0, i1);
    }

    // Obstacles
    Obstacle<Dim> tempObs{};
    for (Json::Value::ArrayIndex i = 0; i < obstacles.size(); i++) {
        // Positions
        for (Json::Value::ArrayIndex j = 0; j < Dim; j++) {
            tempObs.i0[j] = std::find(Msh.Faces[j].begin(), Msh.Faces[j].end(), obstacles[i]["x0"][j].asDouble()) - Msh.Faces[j].begin();
            tempObs.i1[j] = std::find(Msh.Faces[j].begin(), Msh.Faces[j].end(), obstacles[i]["x1"][j].asDouble()) - Msh.Faces[j].end() - 1;
        }

        // bObs
        runLoopMesh<Dim>(Msh.N, [&](std::array<size_t, 3> iX, size_t Ny, size_t Nx) {
                    size_t idx = calcIndex(iX[0], iX[1], Ny, iX[2], Nx); Msh.bObs[idx] = true;
                }, tempObs.i0, tempObs.i1);

        // Control
        Msh.Obs.push_back(std::move(tempObs)); tempObs = {};
    }

    /* if constexpr (Dim == 1) {sections1D(Msh, sections, obstacles);} */ 
    /* else if constexpr (Dim == 2) {sections2D(Msh, sections, obstacles);} */
    /* else if constexpr (Dim == 3) {sections3D(Msh, sections, obstacles);} */

    // Coefficients (nD)
    Msh.matA.resize(Msh.totNodes); Msh.matB.resize(Msh.totNodes, 0); Msh.tempA.resize(Msh.totNodes); Msh.tempB.resize(Msh.totNodes, 0);
}

template <size_t Dim> void MeshFVM<Dim>::deriveMeshBase(MeshSolver<Dim> p, std::array<MeshBase<Dim>, Dim>& V) {
    // Dimension Loop
    for (size_t i = 0; i < Dim; i++) {
        // Nodes
        for (size_t j = 0; j < Dim; j++) { V[i].N[j] = p.N[j]; }
        V[i].N[i] += 1; for (size_t val : V[i].N) {V[i].totNodes *= val;}

        // Resize
        for (size_t j = 0; j < Dim; j++) {
            V[i].Faces[j].resize(V[i].N[j]+1); V[i].Nodes[j].resize(V[i].N[j]); V[i].deltaX[j].resize(V[i].N[j]); V[i].dX[j].resize(V[i].N[j]+1);
        }

        // Geometry
        for (size_t j = 0; j < Dim; j++) {
            // i = j
            if (i == j) { 
                V[i].Nodes[j] = p.Faces[j]; for (size_t k = 0; k < p.Nodes[j].size(); k++) { V[i].Faces[j][k+1] = p.Nodes[j][k]; }
                V[i].Faces[j].front() = V[i].Nodes[j].front() - p.Nodes[j].front(); V[i].Faces[j].back() = 2 * V[i].Nodes[j].back() - p.Nodes[j].back(); continue;
            }

            // i != j
            V[i].Nodes[j] = p.Nodes[j];
            for (size_t k = 0; k < V[i].Nodes[j].size(); k++) {V[i].Nodes[j][k] = p.Nodes[j][k]; V[i].Faces[j][k] = p.Faces[j][k];} V[i].Faces[j].back() = p.Faces[j].back();
        }

        // Deltas
        for (size_t j = 0; j < Dim; j++) {
            for (size_t k = 0; k < V[i].deltaX[j].size(); k++) {
                V[i].deltaX[j][k] = V[i].Faces[j][k+1] - V[i].Faces[j][k];
                if (k == V[i].deltaX[i].size()-1) { continue; } V[i].dX[j][k+1] = V[i].Nodes[j][k+1] - V[i].Nodes[j][k];
            }
            V[i].dX[j].front() = V[i].Nodes[j].front() - V[i].Faces[j].front();
            V[i].dX[j].back() = V[i].Faces[j].back() - V[i].Nodes[j].back();
        }

        // Value
        V[i].Phi.resize(V[i].totNodes); V[i].oPhi.resize(V[i].totNodes);

        // Geometry
        V[i].Vp.resize(V[i].totNodes, 1); for (size_t j = 0; j < Dim; j++) { V[i].S[j].resize(V[i].totNodes, 1); }
        runLoopMesh(V[i].N, [&](std::array<size_t, 3> iX, size_t Ny, size_t Nx) {
                size_t idx = calcIndex(iX[0], iX[1], Ny, iX[2], Nx);
                for (size_t nD = 0; nD < Dim; nD++) {
                    V[i].Vp[idx] *= V[i].deltaX[nD][iX[nD]];
                    for (size_t nDD = 0; nDD < Dim; nDD++) { if (nD != nDD) { V[i].S[nD][idx] *= V[i].deltaX[nDD][iX[nDD]];} }
                } });

        // Coefficients
        V[i].matA.resize(V[i].totNodes); V[i].matB.resize(V[i].totNodes, 0); V[i].oR.resize(V[i].totNodes, 0);
    }
}

///// Boundary Conditions /////
bool isFormula (const std::string value) {
    try { size_t i = 0; std::stod(value, &i); return i != value.length(); }
    catch (...) { return true; }
}

template <size_t Dim> void importInitialConditions(MeshSolver<Dim>& Msh, std::string sPath) {
    // PENDING -- MAKE THIS FUNCTION
    // Import .csv -- Check range and dimension (interpolate if needed) -- Store data
}

template <size_t Dim> void MeshFVM<Dim>::addBoundariesSolver(MeshSolver<Dim>& Msh, Material Mat, Parser& Prs, Json::Value boundaries, double dInit, std::string sInit) {
    // Initial Conditions
    if (!Mat.bPath) {std::fill(Msh.Phi.begin(), Msh.Phi.end(), dInit); std::fill(Msh.oPhi.begin(), Msh.oPhi.end(), dInit);} else {importInitialConditions(Msh, sInit);}

    // Control
    Msh.BC.resize(boundaries.size());

    // Boundaries
    for (Json::Value::ArrayIndex i = 0; i < boundaries.size(); i++) {

        // Position
        Msh.BC[i].side = boundaries[i]["side"].asInt(); if (Msh.BC[i].side != 0 && Msh.BC[i].side != 1) { throw std::invalid_argument("Boundary side not define correctly.\n");}
        for (Json::Value::ArrayIndex j = 0; j < Dim; j++) {
            Msh.BC[i].i0[j] = std::lower_bound(Msh.Nodes[i].begin(), Msh.Nodes[i].end(), boundaries[i]["x0"][j].asDouble() - epsFind) - Msh.Nodes[i].begin();
            Msh.BC[i].i1[j] = std::lower_bound(Msh.Nodes[i].begin(), Msh.Nodes[i].end(), boundaries[i]["x1"][j].asDouble() - epsFind) - Msh.Nodes[i].begin();
        }

        // Value
        if (isFormula(boundaries[i]["value"].asString())) {

            // Control
            Msh.BC[i].expression = boundaries[i]["value"].asString();
            Msh.BC[i].iExpr = Prs.registerExpression(Msh.BC[i].expression);

            // Parser
            if (Msh.BC[i].expression.find(" t ") != std::string::npos) { // Update time
                Msh.BC[i].bUpdate = true; Msh.BC[i].iEq = 0;
                Msh.BC[i].value = Prs.evaluateTime(Msh.BC[i].iEq, 0);
            } else if (Msh.BC[i].expression.find(" x ") != std::string::npos || Msh.BC[i].expression.find(" y ") != std::string::npos || Msh.BC[i].expression.find(" z ") != std::string::npos) {Msh.BC[i].iEq = 1;} // Update coordinates 
            else { throw std::invalid_argument("Equation not recognized: " + Msh.BC[i].expression);}

        } else {Msh.BC[i].value = boundaries[i]["value"].asDouble();}
        

        // Type
        if (boundaries[i]["type"] == "Dirichlet") {Msh.BC[i].type = 0;}
        else if (boundaries[i]["type"] == "Neumann") {Msh.BC[i].type = 1;}
        else if (boundaries[i]["type"] == "Robin") {Msh.BC[i].type = 2;

            // Value
            if (isFormula(boundaries[i]["alpha"].asString())) {
                // Control
                Msh.BC[i].bA = true; Msh.BC[i].expressionA = boundaries[i]["alpha"].asString();
                Msh.BC[i].iExprA = Prs.registerExpression(Msh.BC[i].expressionA);

                // Parser
                if (Msh.BC[i].expressionA.find(" t ") != std::string::npos) { // Update time
                    Msh.BC[i].iEqA = 0; Msh.BC[i].alpha = Prs.evaluateTime(Msh.BC[i].iEqA, 0);
                } else if (Msh.BC[i].expressionA.find(" x ") != std::string::npos || Msh.BC[i].expressionA.find(" y ") != std::string::npos || Msh.BC[i].expressionA.find(" z ") != std::string::npos) {Msh.BC[i].iEqA = 1; if constexpr (Dim == 1) { throw std::invalid_argument("Coordinate functions not supported for 1D\n"); }} // Update coordinates
                else { throw std::invalid_argument("Equation not recognized: " + Msh.BC[i].expressionA); } 
            }

        } else { throw std::invalid_argument("Boundary type not recognized: " + std::to_string(Msh.BC[i].type)); }


        // PENDING -- DEFINE AS LAMBDAS
        /* // Resize */
        /* if constexpr (Dim == 1) {sizeBoundary1D(Msh.BC[i], Msh, Prs);} */
        /* else if constexpr (Dim == 2) {sizeBoundary2D(Msh.BC[i], Msh, Prs);} */
        /* else if constexpr (Dim == 3) {sizeBoundary3D(Msh.BC[i], Msh, Prs);} */


    }
}

    // Was it worth it to pass the std::array<MeshBase<Dim>, Dim>& V ? Or should I just pass MeshBase<Dim>
    // Test with addBoundariesBase() -- The code doesn't know what index to access if I send them independently, would need to pass that as well -- in main() -> for (size_t i = 0; i < Dim; i++) { Msh.addBoundariesBase(i, V, Mat, Prs, boundaries, bInit, sInit) }

template <size_t Dim> void MeshFVM<Dim>::addBoundariesBase(size_t i, std::array<MeshBase<Dim>, Dim>& V, Material Mat, Parser& Prs, Json::Value boundaries, std::vector<double> dInit, std::vector<std::string> sInit) {

}

///// Compiler Instances /////
template class MeshFVM<1>;
template class MeshFVM<2>;
template class MeshFVM<3>;
