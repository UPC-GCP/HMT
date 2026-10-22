// Imports
#include <algorithm>
#include <array>
#include <cstddef>
#include <json/value.h>
#include <optional>
#include <stdexcept>
#include <string>

// Self-Imports
#include "o00_Globals.h"
#include "o02_Mesh.h"
#include "o02_MeshFVM.h"

///// Mesh Generation /////
template <size_t Dim> void MeshFVM<Dim>::calculateFaces(const std::array<size_t, Dim>& cVec, const Json::Value& refData, std::array<std::vector<double>, Dim>& nFaces){
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

template <size_t Dim> void MeshFVM<Dim>::generateMeshSolver(MeshSolver<Dim>& Msh, const Json::Value& qNode, const Json::Value& sections, const Json::Value& refinement, const Json::Value& obstacles) {
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

    // Resize (nD)
    Msh.nMat.resize(Msh.totNodes); for (std::vector<double>& Sk : Msh.S) {Sk.resize(Msh.totNodes, 1);} Msh.Vp.resize(Msh.totNodes, 1); Msh.bObs.resize(Msh.totNodes, false);
    Msh.Phi.resize(Msh.totNodes); Msh.sPhi.resize(Msh.totNodes); Msh.oPhi.resize(Msh.totNodes);
    
    // Sections
    std::array<size_t, Dim> i0{}, i1{}; size_t valMat{}; double valSource{};
    for (Json::Value::ArrayIndex i = 0; i < sections.size(); i++) {
        // Positions
        for (Json::Value::ArrayIndex j = 0; j < Dim; j++) {
            i0[j] = std::find(Msh.Faces[j].begin(), Msh.Faces[j].end(), sections[i]["x0"][j].asDouble()) - Msh.Faces[j].begin();
            i1[j] = std::find(Msh.Faces[j].begin(), Msh.Faces[j].end(), sections[i]["x1"][j].asDouble()) - Msh.Faces[j].begin();
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
            tempObs.i1[j] = std::find(Msh.Faces[j].begin(), Msh.Faces[j].end(), obstacles[i]["x1"][j].asDouble()) - Msh.Faces[j].begin();
        }

        // bObs
        runLoopMesh<Dim>(Msh.N, [&](std::array<size_t, 3> iX, size_t Ny, size_t Nx) {
                    size_t idx = calcIndex(iX[0], iX[1], Ny, iX[2], Nx); Msh.bObs[idx] = true;
                }, tempObs.i0, tempObs.i1);

        // Control
        Msh.Obs.push_back(std::move(tempObs)); tempObs = {};
    }

    // Coefficients (nD)
    Msh.matA.resize(Msh.totNodes); Msh.matB.resize(Msh.totNodes, 0); Msh.tempA.resize(Msh.totNodes); Msh.tempB.resize(Msh.totNodes, 0);
}

template <size_t Dim> void MeshFVM<Dim>::deriveMeshBase(const MeshSolver<Dim>& p, std::array<MeshBase<Dim>, Dim>& V) {
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

template <size_t Dim> void importInitialConditions(MeshBase<Dim>& Msh, std::string sPath) {
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
            if (Msh.BC[i].expression.find(" t ") != std::string::npos) { Msh.BC[i].bUpdate = true; Msh.BC[i].iEq = 0; Msh.BC[i].value = Prs.evaluateTime(Msh.BC[i].iEq, 0); } // Update time
            else if (Msh.BC[i].expression.find(" x ") != std::string::npos || Msh.BC[i].expression.find(" y ") != std::string::npos || Msh.BC[i].expression.find(" z ") != std::string::npos) { Msh.BC[i].iEq = 1; if constexpr (Dim == 1) { throw std::invalid_argument("Coordinate expressions not supported for 1D\n"); } } // Update coordinates 
            else { throw std::invalid_argument("Equation not recognized: " + Msh.BC[i].expression); }
        } else { Msh.BC[i].value = boundaries[i]["value"].asDouble(); }
        
        // Type
        if (boundaries[i]["type"] == "Dirichlet") { Msh.BC[i].type = 0; }
        else if (boundaries[i]["type"] == "Neumann") { Msh.BC[i].type = 1; }
        else if (boundaries[i]["type"] == "Robin") { Msh.BC[i].type = 2;
            // Value
            if (isFormula(boundaries[i]["alpha"].asString())) {
                // Control
                Msh.BC[i].bA = true; Msh.BC[i].expressionA = boundaries[i]["alpha"].asString();
                Msh.BC[i].iExprA = Prs.registerExpression(Msh.BC[i].expressionA);

                // Parser
                if (Msh.BC[i].expressionA.find(" t ") != std::string::npos) { Msh.BC[i].iEqA = 0; Msh.BC[i].alpha = Prs.evaluateTime(Msh.BC[i].iEqA, 0); } // Update time
                else if (Msh.BC[i].expressionA.find(" x ") != std::string::npos || Msh.BC[i].expressionA.find(" y ") != std::string::npos || Msh.BC[i].expressionA.find(" z ") != std::string::npos) { Msh.BC[i].iEqA = 1; if constexpr (Dim == 1) { throw std::invalid_argument("Coordinate functions not supported for 1D\n"); } } // Update coordinates
                else { throw std::invalid_argument("Equation not recognized: " + Msh.BC[i].expressionA); } 
            } else { Msh.BC[i].alpha = boundaries[i]["alpha"].asDouble(); }
        } else { throw std::invalid_argument("Boundary type not recognized: " + boundaries[i]["type"].asString()); }

        // Resize
        size_t tempSize{1}; 
        for (size_t nD = 0; nD < Dim; nD++) {
            //Control
            if (Msh.BC[i].i0[nD] != Msh.BC[i].i1[nD]) { continue; }
            for (size_t nDD = 0; nDD < Dim; nDD++) { if (nDD != nD) { tempSize *= Msh.N[nDD]; } }

            // Value
            if (Msh.BC[i].iEq == 0) { // Time
                // Resize
                Msh.BC[i].Phi.resize(tempSize, Msh.BC[i].value); Msh.BC[i].oPhi.resize(tempSize, Msh.BC[i].value);
                if (Msh.BC[i].type == 2) { Msh.BC[i].A.resize(tempSize, Msh.BC[i].alpha); Msh.BC[i].oA.resize(tempSize, Msh.BC[i].alpha); }
            } else if (Msh.BC[i].iEq == 1) { // Coordinates
                if constexpr (Dim > 1) {
                    // Reduce
                    std::array<size_t, Dim - 1> tempN{}, tempi0{}, tempi1{}, tempX{}; tempN = reduceDimension<Dim>(Msh.N, nD); tempi0 = reduceDimension<Dim>(Msh.BC[i].i0, nD); tempi1 = reduceDimension(Msh.BC[i].i1, nD);
                    // Resize & Parser
                    Msh.BC[i].Phi.resize(tempSize); Msh.BC[i].oPhi.resize(tempSize);
                    runLoopMesh<Dim - 1>(tempN, [&](std::array<size_t, 3> iX, size_t Ny, size_t Nx) {
                            size_t idx = calcIndex(iX[0], iX[1], Ny, iX[2], Nx);
                            std::array<double, 3> iXX{}; for (size_t k = 0, kk = 0; k < Dim; k++) { if (k == nD) { Msh.BC[i].side == 0 ? iXX[k] = Msh.Faces[nD].front() : iXX[k] = Msh.Faces[nD].back(); } iXX[k] = Msh.Nodes[k][iX[kk++]]; }
                            Msh.BC[i].Phi[idx] = Prs.evaluateCoordinates(Msh.BC[i].iExpr, iXX, nD);
                            }, tempi0, tempi1); Msh.BC[i].oPhi = Msh.BC[i].Phi;
                    
                    if (Msh.BC[i].type == 2) {
                        Msh.BC[i].A.resize(tempSize); Msh.BC[i].oA.resize(tempSize);
                        runLoopMesh<Dim - 1>(tempN, [&](std::array<size_t, 3> iX, size_t Ny, size_t Nx) {
                                size_t idx = calcIndex(iX[0], iX[1], Ny, iX[2], Nx);
                                std::array<double, 3> iXX{}; for (size_t k = 0, kk = 0; k < Dim; k++) { if (k == nD) { Msh.BC[i].side == 0 ? iXX[k] = Msh.Faces[nD].front() : iXX[k] = Msh.Faces[nD].back(); } iXX[k] = Msh.Nodes[k][iX[kk++]]; }
                                Msh.BC[i].A[idx] = Prs.evaluateCoordinates(Msh.BC[i].iExprA, iXX, nD);
                                }, tempi0, tempi1); Msh.BC[i].oA = Msh.BC[i].A;
                    }
                }
            } else { // Scalar
                // Resize
                Msh.BC[i].Phi.resize(tempSize, Msh.BC[i].value); Msh.BC[i].oPhi.resize(tempSize, Msh.BC[i].value);
                if (Msh.BC[i].type == 2) { Msh.BC[i].A.resize(tempSize, Msh.BC[i].alpha); Msh.BC[i].oA.resize(tempSize, Msh.BC[i].alpha); }
            }
        }
    }
}

    // Was it worth it to pass the std::array<MeshBase<Dim>, Dim>& V ? Or should I just pass MeshBase<Dim>
    // Test with addBoundariesBase() -- The code doesn't know what index to access if I send them independently, would need to pass that as well -- in main() -> for (size_t i = 0; i < Dim; i++) { Msh.addBoundariesBase(i, V, Mat, Prs, boundaries, bInit, sInit) }

template <size_t Dim> void MeshFVM<Dim>::addBoundariesBase(size_t iV, std::array<MeshBase<Dim>, Dim>& V, const Material& Mat, Parser& Prs, const Json::Value& boundaries, std::vector<double> dInit, std::vector<std::string> sInit) {
    // Initial Conditions
    if (!Mat.bPath) {std::fill(V[iV].Phi.begin(), V[iV].Phi.end(), dInit[iV]); std::fill(V[iV].oPhi.begin(), V[iV].oPhi.end(), dInit[iV]);} else {importInitialConditions(V[iV], sInit[iV]);}

    // Control
    V[iV].BC.resize(boundaries.size());

    // Boundaries
    for (Json::Value::ArrayIndex i = 0; i < boundaries.size(); i++) {
        // Position
        V[iV].BC[i].side = boundaries[i]["side"].asInt(); if (V[iV].BC[i].side != 0 && V[iV].BC[i].side != 1) { throw std::invalid_argument("Boundary side not define correctly.\n");}
        for (Json::Value::ArrayIndex j = 0; j < Dim; j++) {
            V[iV].BC[i].i0[j] = std::lower_bound(V[iV].Nodes[i].begin(), V[iV].Nodes[i].end(), boundaries[i]["x0"][j].asDouble() - epsFind) - V[iV].Nodes[i].begin();
            V[iV].BC[i].i1[j] = std::lower_bound(V[iV].Nodes[i].begin(), V[iV].Nodes[i].end(), boundaries[i]["x1"][j].asDouble() - epsFind) - V[iV].Nodes[i].begin();
        }

        // Value
        if (isFormula(boundaries[i]["value"].asString())) {
            // Control
            V[iV].BC[i].expression = boundaries[i]["value"].asString();
            V[iV].BC[i].iExpr = Prs.registerExpression(V[iV].BC[i].expression);

            // Parser
            if (V[iV].BC[i].expression.find(" t ") != std::string::npos) { V[iV].BC[i].bUpdate = true; V[iV].BC[i].iEq = 0; V[iV].BC[i].value = Prs.evaluateTime(V[iV].BC[i].iEq, 0); } // Update time
            else if (V[iV].BC[i].expression.find(" x ") != std::string::npos || V[iV].BC[i].expression.find(" y ") != std::string::npos || V[iV].BC[i].expression.find(" z ") != std::string::npos) { V[iV].BC[i].iEq = 1; if constexpr (Dim == 1) { throw std::invalid_argument("Coordinate expressions not supported for 1D\n"); } } // Update coordinates 
            else { throw std::invalid_argument("Equation not recognized: " + V[iV].BC[i].expression); }
        } else { V[iV].BC[i].value = boundaries[i]["value"].asDouble(); }
        
        // Type
        if (boundaries[i]["type"] == "Dirichlet") { V[iV].BC[i].type = 0; }
        else if (boundaries[i]["type"] == "Neumann") { V[iV].BC[i].type = 1; }
        else { throw std::invalid_argument("Boundary type not recognized: " + boundaries[i]["type"].asString()); }

        // Resize
        size_t tempSize{1}; 
        for (size_t nD = 0; nD < Dim; nD++) {
            //Control
            if (V[iV].BC[i].i0[nD] != V[iV].BC[i].i1[nD]) { continue; }
            for (size_t nDD = 0; nDD < Dim; nDD++) { if (nDD != nD) { tempSize *= V[iV].N[nDD]; } }

            // Value
            if (V[iV].BC[i].iEq == 0) { // Time
                // Resize
                V[iV].BC[i].Phi.resize(tempSize, V[iV].BC[i].value); V[iV].BC[i].oPhi.resize(tempSize, V[iV].BC[i].value);
                if (V[iV].BC[i].type == 2) { V[iV].BC[i].A.resize(tempSize, V[iV].BC[i].alpha); V[iV].BC[i].oA.resize(tempSize, V[iV].BC[i].alpha); }
            } else if (V[iV].BC[i].iEq == 1) { // Coordinates
                if constexpr (Dim > 1) {
                    // Reduce
                    std::array<size_t, Dim - 1> tempN{}, tempi0{}, tempi1{}, tempX{}; tempN = reduceDimension<Dim>(V[iV].N, nD); tempi0 = reduceDimension<Dim>(V[iV].BC[i].i0, nD); tempi1 = reduceDimension(V[iV].BC[i].i1, nD);
                    // Resize & Parser
                    V[iV].BC[i].Phi.resize(tempSize); V[iV].BC[i].oPhi.resize(tempSize);
                    runLoopMesh<Dim - 1>(tempN, [&](std::array<size_t, 3> iX, size_t Ny, size_t Nx) {
                            size_t idx = calcIndex(iX[0], iX[1], Ny, iX[2], Nx);
                            std::array<double, 3> iXX{}; for (size_t k = 0, kk = 0; k < Dim; k++) { if (k == nD) { V[iV].BC[i].side == 0 ? iXX[k] = V[iV].Faces[nD].front() : iXX[k] = V[iV].Faces[nD].back(); } iXX[k] = V[iV].Nodes[k][iX[kk++]]; }
                            V[iV].BC[i].Phi[idx] = Prs.evaluateCoordinates(V[iV].BC[i].iExpr, iXX, nD);
                            }, tempi0, tempi1); V[iV].BC[i].oPhi = V[iV].BC[i].Phi;
                }
            } else { // Scalar
                // Resize
                V[iV].BC[i].Phi.resize(tempSize, V[iV].BC[i].value); V[iV].BC[i].oPhi.resize(tempSize, V[iV].BC[i].value);
                if (V[iV].BC[i].type == 2) { V[iV].BC[i].A.resize(tempSize, V[iV].BC[i].alpha); V[iV].BC[i].oA.resize(tempSize, V[iV].BC[i].alpha); }
            }
        }
    }
}

///// Compiler Instances /////
template class MeshFVM<1>;
template class MeshFVM<2>;
template class MeshFVM<3>;
