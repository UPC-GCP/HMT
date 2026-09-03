/* #include <algorithm> */
#include <algorithm>
#include <cstddef>
/* #include <cstdint> */
#include <cstdio>
#include <cstdlib>
#include <iostream>
/* #include <utility> */
/* #include <iterator> */
#include <json/json.h>
#include <json/value.h>
#include <cmath>
#include <stdexcept>
/* #include <string> */
/* #include <strings.h> */
/* #include <vector> */

// Self-Imports
#include "o01_Material.h"
#include "o02_Mesh.h"


////////// Mesh Generation //////////
template <size_t Dim> void Mesh<Dim>::calculateFaces(std::array<size_t, Dim> cVec, Json::Value refData, std::array<std::vector<double>, Dim>& nFaces){
    // Control
    std::string sAlgo = refData["algorithm"].asString();
    size_t iAxis = refData["axis"].asInt(), NSec = refData["N"].asInt(), cNode = cVec[iAxis];
    double x0 = refData["range"][0].asDouble(), x1 = refData["range"][1].asDouble(); double length = x1 - x0;

    // Face Positions
    if (sAlgo == "Bidirectional"){ // Bidirectional (A, xC)
        
        // Mesh Parameters
        double A = refData["strength"].asDouble(), xC = refData["centering"].asDouble() - x0; 
        if (A > 2.5 || A < -2.5){std::cerr << "Bidirectional mesh parameter out of range. (strength)\n"; throw std::invalid_argument("Check .json");}
        if (xC <= 0 || xC >= length){std::cerr << "Bidirectional mesh parameter out of range. (centering)\n"; throw std::invalid_argument("Check .json");}
        // A (+): Closer to center, A (-): Away from center
        
        // Generate Mesh
        for (int i = cNode; i < cNode+NSec+1; i++) {
            nFaces[iAxis][i] = x0 + (static_cast<double>(i)-cNode) * length / NSec + A * (xC - (static_cast<double>(i)-cNode) * length / NSec) * (1 - (static_cast<double>(i)-cNode)/NSec) * (static_cast<double>(i)-cNode) / NSec;
        }

    } else if (sAlgo == "PowerLaw"){ // Power-Law (Kappa)
        
        // Mesh Parameters
        double kappa = refData["kappa"].asDouble();
        if (kappa <= 0){std::cerr << "Power-Law mesh parameter out of range.\n"; std::exit(0);}
        // kappa < 1: Shifts right, kappa > 1: Shifts left

        // Generate Mesh
        for (int i = cNode; i < cNode+NSec+1; i++) {
            nFaces[iAxis][i] = x0 + length * pow(((static_cast<double>(i)-cNode) / NSec), kappa);
        }

    } else if(sAlgo == "HyperSingle"){ // Hyperbolic Tangent (Single-Sided)
        
        // Mesh Parameters
        double delta = refData["delta"].asDouble(); double A{}, B{};
        if (delta == 0){std::cerr << "Hyperbolic tangent mesh parameter out of range.\n"; throw std::invalid_argument("Check .json");}
        // delta: independent of sign (+/-), higher value pushes towards the left

        // Generate Mesh
        for (int i = cNode; i < cNode+NSec+1; i++){
            A = tanh(delta * ((static_cast<double>(i) - cNode) / NSec - 1)); B = tanh(delta);
            nFaces[iAxis][i] = x0 + length * (1 + A / B);
        }

    } else if(sAlgo == "HyperDouble"){ // Hyperbolic Tangent (Double-Sided)
        
        // Mesh Parameters
        double delta = refData["delta"].asDouble(); double A{}, B{};
        if (delta == 0){std::cerr << "Hyperbolic tangent mesh parameter out of range.\n"; throw std::invalid_argument("Check .json");}
        // delta: independent of sign (+/-), higher value pushes towards the edges

        // Generate Mesh
        for (int i = cNode; i < cNode+NSec+1; i++){
            A = tanh(delta*((static_cast<double>(i) - cNode)/NSec - 0.5)); B = tanh(0.5 * delta);
            nFaces[iAxis][i] = x0 + 0.5 * length * (1 + A/B);
        }

    } else if (sAlgo == "Exponential") { // Exponential

        // Mesh Parameters
        double kappa = refData["kappa"].asDouble();
        if (kappa == 0){std::cerr << "Exponential mesh parameter out of range.\n"; throw  std::invalid_argument("Check .json");}
        // kappa < 0: Shifts right, kappa > 0: Shifts left 
        
        // Generate Mesh
        for (int i = cNode; i < cNode+NSec+1; i++){
            nFaces[iAxis][i] = x0 + length * (exp(kappa * (static_cast<double>(i) - cNode) / NSec) - 1) / (exp(kappa) - 1);
        }

    } else {std::cerr << "Algorithm not recognized.\n";}
}

void sections1D(MeshSolver<1>& Msh, Json::Value sections, Json::Value obstacles){
    // Sections
    std::vector<size_t> Pos0(Msh.N.size()), Pos1(Msh.N.size()); size_t idx{};
    for (Json::Value::ArrayIndex i = 0; i < sections.size(); i++){
        // Find Positions (nD)
        for (int j = 0; j < Msh.N.size(); j++){
            Pos0[j] = std::lower_bound(Msh.Faces[j].begin(), Msh.Faces[j].end(), sections[i]["x0"][j].asDouble()) - Msh.Faces[j].begin();
            Pos1[j] = std::lower_bound(Msh.Faces[j].begin(), Msh.Faces[j].end(), sections[i]["x1"][j].asDouble()) - Msh.Faces[j].begin();
        }

        // Internal Nodes (Non-nD)
        for (size_t j = Pos0[0]; j < Pos1[0]; j++){
            // Control
            idx = calcIndex(j);

            // Material
            Msh.nMat[idx] = sections[i]["material"].asInt(); Msh.sPhi[idx] = sections[i]["source"].asDouble();

            // Geometry
            Msh.S[0][idx] = Msh.deltaX[0][j]; Msh.Vp[idx] = Msh.deltaX[0][j]; 
        }
    }
    
    // Obstacles
    if (obstacles.size() != 0) {std::cerr << "Obstacles not admitted for 1D.\n"; throw std::invalid_argument("Check .json\n");}
}

void sections2D(MeshSolver<2>& Msh, Json::Value sections, Json::Value obstacles){
    // Sections
    std::vector<size_t> Pos0(Msh.N.size()), Pos1(Msh.N.size()); size_t idx{};
    for (Json::Value::ArrayIndex i = 0; i < sections.size(); i++){
        // Find Positions (nD)
        for (int j = 0; j < Msh.N.size(); j++){
            Pos0[j] = std::lower_bound(Msh.Faces[j].begin(), Msh.Faces[j].end(), sections[i]["x0"][j].asDouble()) - Msh.Faces[j].begin();
            Pos1[j] = std::lower_bound(Msh.Faces[j].begin(), Msh.Faces[j].end(), sections[i]["x1"][j].asDouble()) - Msh.Faces[j].begin();
        }

        // Internal Nodes (Non-nD)
        for (size_t j = Pos0[0]; j < Pos1[0]; j++){
            for (size_t k = Pos0[1]; k < Pos1[1]; k++){
                // Control
                idx = calcIndex(j, k, Msh.N[1]);

                // Material
                Msh.nMat[idx] = sections[i]["material"].asInt(); Msh.sPhi[idx] = sections[i]["source"].asDouble();

                // Geometry
                Msh.S[0][idx] = Msh.deltaX[1][k]; Msh.S[1][idx] = Msh.deltaX[0][j]; Msh.Vp[idx] = Msh.deltaX[0][j] * Msh.deltaX[1][k];
            }
        }
    }

    // Obstacles
    Obstacle<2> tempObs{}; 
    for (Json::Value::ArrayIndex i = 0; i < obstacles.size(); i++){
        // Find Positions (nD)
        for (int j = 0; j < Msh.N.size(); j++){
            tempObs.i0[j] = std::lower_bound(Msh.Nodes[j].begin(), Msh.Nodes[j].end(), obstacles[i]["x0"][j].asDouble()) - Msh.Nodes[j].begin();
            tempObs.i1[j] = std::lower_bound(Msh.Nodes[j].begin(), Msh.Nodes[j].end(), obstacles[i]["x1"][j].asDouble()) - Msh.Nodes[j].begin();
        }

        // bObs
        for (size_t j = tempObs.i0[0]; j < tempObs.i1[0]; j++){
            for (size_t k = tempObs.i0[1]; k < tempObs.i1[1]; k++){
                    idx = calcIndex(j, k, Msh.N[1]); Msh.bObs[idx] = true;
            }
        }

        // Control
        Msh.Obs.push_back(std::move(tempObs)); tempObs = {}; 
    }
}

void sections3D(MeshSolver<3>& Msh, Json::Value sections, Json::Value obstacles){
    // Sections
    std::vector<size_t> Pos0(Msh.N.size()), Pos1(Msh.N.size()); size_t idx{}, valMat{}; double valSource{};
    for (Json::Value::ArrayIndex i = 0; i < sections.size(); i++){
        // Control
        valMat = sections[i]["material"].asInt(); valSource = sections[i]["source"].asDouble();

        // Find Positions (nD)
        for (int j = 0; j < Msh.N.size(); j++){
            Pos0[j] = std::lower_bound(Msh.Faces[j].begin(), Msh.Faces[j].end(), sections[i]["x0"][j].asDouble()) - Msh.Faces[j].begin();
            Pos1[j] = std::lower_bound(Msh.Faces[j].begin(), Msh.Faces[j].end(), sections[i]["x1"][j].asDouble()) - Msh.Faces[j].begin();
        }

        // Internal Nodes (Non-nD)
        for (size_t l = Pos0[2]; l < Pos1[2]; l++){
            for (size_t j = Pos0[0]; j < Pos1[0]; j++){
                for (size_t k = Pos0[1]; k < Pos1[1]; k++){
                    // Control
                    idx = calcIndex(j, k, Msh.N[1], l, Msh.N[0]);

                    // Material
                    Msh.nMat[idx] = valMat; Msh.sPhi[idx] = valSource;

                    // Geometry
                    Msh.S[0][idx] = Msh.deltaX[1][k] * Msh.deltaX[2][l]; Msh.S[1][idx] = Msh.deltaX[0][j] * Msh.deltaX[2][l]; Msh.S[2][idx] = Msh.deltaX[0][j] * Msh.deltaX[1][k]; Msh.Vp[idx] = Msh.deltaX[0][j] * Msh.deltaX[1][k] * Msh.deltaX[2][l]; 
                }
            }
        }

    }

    Obstacle<3> tempObs{}; 
    for (Json::Value::ArrayIndex i = 0; i < obstacles.size(); i++){
        // Find Positions (nD)
        for (int j = 0; j < Msh.N.size(); j++){
            tempObs.i0[j] = std::lower_bound(Msh.Nodes[j].begin(), Msh.Nodes[j].end(), obstacles[i]["x0"][j].asDouble()) - Msh.Nodes[j].begin();
            tempObs.i1[j] = std::lower_bound(Msh.Nodes[j].begin(), Msh.Nodes[j].end(), obstacles[i]["x1"][j].asDouble()) - Msh.Nodes[j].begin();
        } 

        // bObs
        for (size_t l = tempObs.i0[2]; l < tempObs.i1[2]; l++){
            for (size_t j = tempObs.i0[0]; j < tempObs.i1[0]; j++){
                for (size_t k = tempObs.i0[1]; k < tempObs.i1[1]; k++){
                    idx = calcIndex(j, k, Msh.N[1], l, Msh.N[0]); Msh.bObs[idx] = true;
                }
            }
        }

        // Control
        Msh.Obs.push_back(std::move(tempObs)); tempObs = {}; 
    }
}

template <size_t Dim> void Mesh<Dim>::generateMeshSolver(MeshSolver<Dim>& Msh, Material Mat, Json::Value qNode, Json::Value sections, Json::Value refinement, Json::Value obstacles){
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
    Msh.nMat.resize(Msh.totNodes); for (std::vector<double>& Sk : Msh.S) {Sk.resize(Msh.totNodes);} Msh.Vp.resize(Msh.totNodes); Msh.bObs.resize(Msh.totNodes);
    Msh.Phi.resize(Msh.totNodes); Msh.sPhi.resize(Msh.totNodes); Msh.oPhi.resize(Msh.totNodes);

    // Sections
    if constexpr (Dim == 1) {sections1D(Msh, sections, obstacles);} 
    else if constexpr (Dim == 2) {sections2D(Msh, sections, obstacles);}
    else if constexpr (Dim == 3) {sections3D(Msh, sections, obstacles);}

    // Coefficients (nD)
    Msh.matA.resize(Msh.totNodes); Msh.matB.resize(Msh.totNodes, 0); Msh.tempA.resize(Msh.totNodes); Msh.tempB.resize(Msh.totNodes, 0);
}

///// Boundary Conditions /////

/* bool isFormula(std::string value){ */
/*     // Stringstream */
/*     std::stringstream ss; ss << value; */

/*     // Check */
/*     float num = 0; ss >> num; */

/*     // Return */
/*     if (ss.good()) {return true;} */
/*     else if (num == 0 && value[0] != 0) {return true;} */
/*     else {return false;} */
/* } */

bool isFormula (const std::string value) {
    try { size_t i = 0; std::stod(value, &i); return i != value.length(); }
    catch (...) { return true; }
}

template <size_t Dim> void importInitialConditions(MeshSolver<Dim>& Msh, std::string sPath) {
    // Import .csv
    // Control (Check dimensions -- Just compare totNodes)
    // Store Data
}

void sizeBoundary1D(Boundary<1>& BC, MeshBase<1> Msh, Parser& Prs) {
    if (BC.i0[0] == BC.i1[0]) {

        // X Boundary
        if (BC.bUpdate) {

            // Time
            if (BC.iEq == 0) {
                BC.Phi.resize(1, BC.value); BC.oPhi.resize(1, BC.value);
                if (BC.type == 2) {BC.A.resize(1, BC.alpha); BC.oA.resize(1, BC.alpha);}
            }

        } else {
            // Scalar
            BC.Phi.resize(1, BC.value); BC.oPhi.resize(1, BC.value);
            if (BC.type == 2) {BC.A.resize(1, BC.alpha); BC.oA.resize(1, BC.alpha);}
        }

    }
}

void sizeBoundary2D(Boundary<2>& BC, MeshBase<2> Msh, Parser& Prs) {
    if (BC.i0[0] == BC.i1[0]) {
        
        // X Boundary 
        if (BC.bUpdate) {

            // Time
            if (BC.iEq == 0) {
                BC.Phi.resize(Msh.N[1], BC.value); BC.oPhi.resize(Msh.N[1], BC.value);
                if (BC.type == 2) {BC.A.resize(Msh.N[1], BC.alpha); BC.oA.resize(Msh.N[1], BC.alpha);}
            }

            // Coordinates
            else if (BC.iEq == 1) {
                BC.Phi.resize(Msh.N[1]); BC.Phi.resize(Msh.N[1]);
                for (size_t j = BC.i0[1]; j < BC.i1[1]; j++) {BC.Phi[j] = Prs.evaluateCoordinates(BC.iExpr, BC.side == 0 ? Msh.Faces[0].front() : Msh.Faces[0].back(), Msh.Nodes[1][j]);} BC.oPhi = BC.Phi;

                if (BC.type == 2) {
                    BC.A.resize(Msh.N[1]); BC.oA.resize(Msh.N[1]);
                    for (size_t j = BC.i0[1];  j < BC.i1[1]; j++) {BC.A[j] = Prs.evaluateCoordinates(BC.iExprA, BC.side == 0 ? Msh.Faces[0].front() : Msh.Faces[0].back(), Msh.Nodes[1][j]); BC.oA = BC.A;
                    }
                }
            }

        } else {
            // Scalar 
            BC.Phi.resize(Msh.N[1], BC.value); BC.oPhi.resize(Msh.N[1], BC.value);
            if (BC.type == 2) {BC.A.resize(Msh.N[1], BC.alpha); BC.oA.resize(Msh.N[1], BC.alpha);}
        }

    } else if (BC.i0[1] == BC.i1[1]) {
        
        // Y Boundary
        if (BC.bUpdate) {

            // Time
            if (BC.iEq == 0) {
                BC.Phi.resize(Msh.N[0], BC.value); BC.oPhi.resize(Msh.N[0], BC.value);
                if (BC.type == 2) {BC.A.resize(Msh.N[0], BC.alpha); BC.oA.resize(Msh.N[0], BC.alpha);}
            }

            // Coordinates
            else if (BC.iEq == 1) {
                BC.Phi.resize(Msh.N[0]); BC.oPhi.resize(Msh.N[0]);
                for (size_t i = BC.i0[0]; i < BC.i1[0]; i++) {BC.Phi[i] = Prs.evaluateCoordinates(BC.iExpr, Msh.Nodes[0][i], BC.side == 0 ? Msh.Faces[1].front() : Msh.Faces[1].back());} BC.oPhi = BC.Phi;

                if (BC.type == 2) {
                    BC.A.resize(Msh.N[0]); BC.oA.resize(Msh.N[0]);
                    for (size_t i = BC.i0[0]; i < BC.i1[0]; i++) {BC.A[i] = Prs.evaluateCoordinates(BC.iExprA, Msh.Nodes[0][i], BC.side == 0 ? Msh.Faces[1].front() : Msh.Faces[1].back());} BC.oA = BC.A;
                }
            }

        } else {
            // Scalar 
            BC.Phi.resize(Msh.N[0], BC.value); BC.oPhi.resize(Msh.N[0], BC.value);
            if (BC.type == 2) {BC.A.resize(Msh.N[0], BC.alpha); BC.oA.resize(Msh.N[0], BC.alpha);}
        }

    }
}

void sizeBoundary3D(Boundary<3>& BC, MeshBase<3> Msh, Parser& Prs) {
    if (BC.i0[0] == BC.i1[0]) {
        
        // X Boundary
        if (BC.bUpdate) {

            // Time
            if (BC.iEq == 0) {
                BC.Phi.resize(Msh.N[1] * Msh.N[2], BC.value); BC.oPhi.resize(Msh.N[1] * Msh.N[2], BC.value);
                if (BC.type == 2) {BC.A.resize(Msh.N[1] * Msh.N[2], BC.alpha); BC.oA.resize(Msh.N[1] * Msh.N[2], BC.alpha);}
            }

            // Coordinates
            else if (BC.iEq == 1) {
                BC.Phi.resize(Msh.N[1] * Msh.N[2]); BC.oPhi.resize(Msh.N[1] * Msh.N[2]);
                for (size_t j = BC.i0[1]; j < BC.i1[1]; j++) {
                    for (size_t k = BC.i0[2]; k < BC.i1[2]; k++) {
                        BC.Phi[calcIndex(j, k, Msh.N[2])] = Prs.evaluateCoordinates(BC.iExpr, BC.side == 0 ? Msh.Faces[0].front() : Msh.Faces[0].back(), Msh.Nodes[1][j], Msh.Nodes[2][k]);
                    } BC.oPhi = BC.Phi;
                }

                if (BC.type == 2) {
                    BC.A.resize(Msh.N[1] * Msh.N[2]); BC.oA.resize(Msh.N[1] * Msh.N[2]);
                    for (size_t j = BC.i0[1]; j < BC.i1[1]; j++) {
                        for (size_t k = BC.i0[2]; k < BC.i1[2]; k++) {
                            BC.Phi[calcIndex(j, k, Msh.N[2])] = Prs.evaluateCoordinates(BC.iExpr, BC.side = 0 ? Msh.Faces[0].front() : Msh.Faces[0].back(), Msh.Nodes[1][j], Msh.Nodes[2][k]);
                        } BC.oPhi = BC.Phi;
                    }
                }
            }

        } else {
            // Scalar 
            BC.Phi.resize(Msh.N[1] * Msh.N[2], BC.value); BC.oPhi.resize(Msh.N[1] * Msh.N[2], BC.value);
            if (BC.type == 2) {BC.A.resize(Msh.N[1] * Msh.N[2], BC.alpha); BC.oA.resize(Msh.N[1] * Msh.N[2], BC.alpha);}
        }

    } else if (BC.i0[1] == BC.i1[1]) {
        
        // Y Boundary
        if (BC.bUpdate) {

            // Time
            if (BC.iEq == 0) {
                BC.Phi.resize(Msh.N[0] * Msh.N[2], BC.value); BC.oPhi.resize(Msh.N[0] * Msh.N[2], BC.value);
                if (BC.type == 2) {BC.A.resize(Msh.N[0] * Msh.N[2], BC.alpha); BC.oA.resize(Msh.N[0] * Msh.N[2], BC.alpha);}
            }

            // Coordinates
            else if (BC.iEq == 1) {
                BC.Phi.resize(Msh.N[0] * Msh.N[2]); BC.oPhi.resize(Msh.N[0] * Msh.N[2]);
                for (size_t i = BC.i0[0]; i < BC.i1[0]; i++) {
                    for (size_t k = BC.i0[2]; k < BC.i1[2]; k++) {
                        BC.Phi[calcIndex(i, k, Msh.N[2])] = Prs.evaluateCoordinates(BC.iExpr, Msh.Nodes[0][i], BC.side == 0 ? Msh.Faces[1].front() : Msh.Faces[1].back(), Msh.Nodes[2][k]);
                    } BC.oPhi = BC.Phi;
                }

                if (BC.type == 2) {
                    BC.A.resize(Msh.N[0] * Msh.N[2]); BC.oA.resize(Msh.N[0] * Msh.N[2]);
                    for (size_t i = BC.i0[0]; i < BC.i1[0]; i++) {
                        for (size_t k = BC.i0[2]; k < BC.i1[2]; k++) {
                            BC.A[calcIndex(i, k, Msh.N[2])] = Prs.evaluateCoordinates(BC.iExprA, Msh.Nodes[0][i], BC.side == 0 ? Msh.Faces[1].front() : Msh.Faces[2].back(), Msh.Nodes[2][k]);
                        } BC.oA = BC.A;
                    }
                }
            }

        } else {
            // Scalar
            BC.Phi.resize(Msh.N[0] * Msh.N[2], BC.value); BC.oPhi.resize(Msh.N[0] * Msh.N[2], BC.value);
            if (BC.type == 2) {BC.A.resize(Msh.N[0] * Msh.N[2], BC.alpha); BC.oA.resize(Msh.N[0] * Msh.N[2], BC.alpha);}
        }

    } else if (BC.i0[2] == BC.i1[2]) {
        
        // Z Boundary
        if (BC.bUpdate) {

            // Time
            if (BC.iEq == 0) {
                BC.Phi.resize(Msh.N[0] * Msh.N[1], BC.value); BC.oPhi.resize(Msh.N[0] * Msh.N[1], BC.value);
                if (BC.type == 2) {BC.A.resize(Msh.N[0] * Msh.N[1], BC.alpha); BC.oA.resize(Msh.N[0] * Msh.N[1], BC.alpha);}
            }

            // Coordinates
            else if (BC.iEq == 1) {
                BC.Phi.resize(Msh.N[0] * Msh.N[1]); BC.oPhi.resize(Msh.N[0] * Msh.N[1]);
                for (size_t i = BC.i0[0]; i < BC.i1[0]; i++) {
                    for (size_t j = BC.i0[1]; j < BC.i1[1]; j++) {
                        BC.Phi[calcIndex(i, j, Msh.N[1])] = Prs.evaluateCoordinates(BC.iExpr, Msh.Nodes[0][i], Msh.Nodes[1][j], BC.side == 0 ? Msh.Faces[2].front() : Msh.Faces[2].back());
                    }
                }

                if (BC.type == 2) {
                    BC.A.resize(Msh.N[0] * Msh.N[1]); BC.oA.resize(Msh.N[0] * Msh.N[1]);
                    for (size_t i = BC.i0[0]; i < BC.i1[0]; i++) {
                        for (size_t j = BC.i0[1]; j < BC.i1[1]; j++) {
                            BC.A[calcIndex(i, j, Msh.N[1])] = Prs.evaluateCoordinates(BC.iExpr, Msh.Nodes[0][i], Msh.Nodes[1][j], BC.side == 0 ? Msh.Faces[2].front() : Msh.Faces[2].back());
                        }
                    }
                }
            }

        } else {
            // Scalar
            BC.Phi.resize(Msh.N[0] * Msh.N[1], BC.value); BC.oPhi.resize(Msh.N[0] * Msh.N[1], BC.value);
            if (BC.type == 2) {BC.A.resize(Msh.N[0] * Msh.N[1], BC.alpha); BC.oA.resize(Msh.N[0] * Msh.N[1], BC.alpha);}
        }

    }
}

template <size_t Dim> void Mesh<Dim>::addBoundariesSolver(MeshSolver<Dim>& Msh, Material Mat, Parser& Prs, Json::Value boundaries, double dInit, std::string sInit) {
    // Initial Conditions
    if (Mat.bPath) {std::fill(Msh.Phi.begin(), Msh.Phi.end(), dInit); std::fill(Msh.oPhi.begin(), Msh.oPhi.end(), dInit);} else {importInitialConditions(Msh, sInit);}

    // Control
    Msh.BC.resize(boundaries.size());

    // Boundaries
    for (Json::Value::ArrayIndex i = 0; i < boundaries.size(); i++) {

        // Position
        Msh.BC[i].side = boundaries[i]["side"].asInt(); if (Msh.BC[i].side != 0 && Msh.BC[i].side != 1) {std::cerr << "Side not defined properly\n"; throw std::invalid_argument("Check .json\n");}
        for (Json::Value::ArrayIndex j = 0; j < Dim; j++) {
            Msh.BC[i].i0[j] = std::lower_bound(Msh.Nodes[i].begin(), Msh.Nodes[i].end(), boundaries[i]["x0"][j].asDouble() - epsFind) - Msh.Nodes[i].begin();
            Msh.BC[i].i1[j] = std::lower_bound(Msh.Nodes[i].begin(), Msh.Nodes[i].end(), boundaries[i]["x1"][j].asDouble() - epsFind) - Msh.Nodes[i].begin();
        }

        std::cout << "Boundary: " << i << " - "; for (size_t val : Msh.BC[i].i0) {std::cout << val << " ";} std::cout << " - "; for (size_t val : Msh.BC[i].i1) {std::cout << val << " ";} std::cout << "\n";

        // Value
        if (isFormula(boundaries[i]["value"].asString())) {

            std::cout << "Formula: " << boundaries[i]["value"].asString().size() << " " << boundaries[i]["value"] << "\n";

            // Control
            Msh.BC[i].bUpdate = true; Msh.BC[i].expression = boundaries[i]["value"].asString();
            Msh.BC[i].iExpr = Prs.registerExpression(Msh.BC[i].expression);

            // Parser
            if (Msh.BC[i].expression.find(" t ") != std::string::npos) { // Update time
                Msh.BC[i].iEq = 0; Msh.BC[i].value = Prs.evaluateTime(Msh.BC[i].iEq, 0);
            } else if (Msh.BC[i].expression.find(" x ") != std::string::npos || Msh.BC[i].expression.find(" y ") != std::string::npos || Msh.BC[i].expression.find(" z ") != std::string::npos) {Msh.BC[i].iEq = 1;} // Update coordinates 
            else {std::cerr << "Equation not recognized: " << Msh.BC[i].expression << "\n"; throw std::invalid_argument("Check .json");}
        } else {Msh.BC[i].value = boundaries[i]["value"].asDouble();}

        // Type
        if (boundaries[i]["type"] == "Dirichlet") {Msh.BC[i].type = 0;}
        else if (boundaries[i]["type"] == "Neumann") {Msh.BC[i].type = 1;}
        else if (boundaries[i]["type"] == "Robin") {Msh.BC[i].type = 2;

            std::cout << "Alpha\n";

            // Value
            if (isFormula(boundaries[i]["alpha"].asString())) {
                // Control
                Msh.BC[i].bA = true; Msh.BC[i].expressionA = boundaries[i]["alpha"].asString();
                Msh.BC[i].iExprA = Prs.registerExpression(Msh.BC[i].expressionA);

                // Parser
                if (Msh.BC[i].expressionA.find(" t ") != std::string::npos) { // Update time
                    Msh.BC[i].iEqA = 0; Msh.BC[i].alpha = Prs.evaluateTime(Msh.BC[i].iEqA, 0);
                } else if (Msh.BC[i].expressionA.find(" x ") != std::string::npos || Msh.BC[i].expressionA.find(" y ") != std::string::npos || Msh.BC[i].expressionA.find(" z ") != std::string::npos) {Msh.BC[i].iEqA = 1; if constexpr (Dim == 1) {std::cerr << "Coordinate functions not supported for 1D\n"; throw std::invalid_argument("Check .json");}} // Update coordinates
                else {std::cerr << "Equation not recognized: " << Msh.BC[i].expressionA << "\n"; throw std::invalid_argument("Check .json");} 
            }

        } else {std::cerr << "Boundary type not recognized.\n"; throw std::invalid_argument("Check .json\n");}

        // Resize
        if constexpr (Dim == 1) {sizeBoundary1D(Msh.BC[i], Msh, Prs);}
        if constexpr (Dim == 2) {sizeBoundary2D(Msh.BC[i], Msh, Prs);}
        else if constexpr (Dim == 3) {sizeBoundary3D(Msh.BC[i], Msh, Prs);}

    }

}

// Compiler Instances
template class Mesh<1>;
template class Mesh<2>;
template class Mesh<3>;
