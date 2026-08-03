/* #include <algorithm> */
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
/* #include <utility> */
#include <vector>
#include <json/json.h>
#include <json/value.h>
#include <cmath>

// Self-Imports
#include "o02_Mesh3D.h"


bool isFormula(std::string value){

    // Stringstream
    std::stringstream ss; ss << value;

    // Check
    float num = 0; ss >> num;

    // Return
    if (ss.good()){
        return true;
    } else if (num == 0 && value[0] != 0){
        return true;
    } else {
        return false;
    }

}


Mesh::Mesh() {

    // UPDATE
    // Will no longer update algorithm and parameters for the entire domain
    // Have it defined for each region and only read during evaluation

    /* // Algorithm Selection */
    /* this->algorithm = algorithm; */

    /* // Bidirectional Non-uniform */
    /* if (true) {this->A = A; this->xC = xC;} else {std::cerr << "Bidirectional Non-uniform";} */

    /* // Unidirectional Non-uniform */
    /* if (true) {this->kappa = kappa;} else {std::cerr << "Unidirectional Non-uniform";} */

    /* // Hyperbolic Tangent (Single / Double) */
    /* if (true) {this->delta = delta;} else {std::cerr << "Hyperbolic Tangent";} */
    
}


void Mesh::calculateFaces(std::vector<size_t> cVec, Json::Value refData, std::vector<std::vector<double>>& nFaces){

    // General
    std::string sAlgo = refData["algorithm"].asString();
    size_t iAxis = refData["axis"].asInt(), NSec = refData["N"].asInt(); size_t cNode = cVec[iAxis];
    double x0 = refData["range"][0].asDouble(), x1 = refData["range"][1].asDouble(); double length = x1 - x0;

    // Face Positions
    if (sAlgo == "Bidirectional"){ // Face Positions 0: Bidirectional (A, xC)
        
        // Mesh Parameters
        double A = refData["strength"].asDouble(), xC = refData["centering"].asDouble() - x0; 
        if (A > 2.5 || A < -2.5){std::cerr << "Bidirectional mesh parameter out of range. (strength)\n"; std::exit(0);}
        if (xC <= 0 || xC >= length){std::cerr << "Bidirectional mesh parameter out of range. (centering)\n"; std::exit(0);}
        // A (+): Closer to center, A (-): Away from center
        
        // Generate Mesh
        for (int i = cNode; i < cNode+NSec+1; i++) {
            nFaces[iAxis][i] = x0 + (static_cast<double>(i)-cNode) * length / NSec + A * (xC - (static_cast<double>(i)-cNode) * length / NSec) * (1 - (static_cast<double>(i)-cNode)/NSec) * (static_cast<double>(i)-cNode) / NSec;
        }

    } else if (sAlgo == "PowerLaw"){ // Face Positions 1: Power-Law (Kappa)
        
        // Mesh Parameters
        double kappa = refData["kappa"].asDouble();
        if (kappa <= 0){std::cerr << "Power-Law mesh parameter out of range.\n"; std::exit(0);}
        // kappa < 1: Shifts right, kappa > 1: Shifts left

        // Generate Mesh
        for (int i = cNode; i < cNode+NSec+1; i++) {
            nFaces[iAxis][i] = x0 + length * pow(((static_cast<double>(i)-cNode) / NSec), kappa);
        }

    } else if(sAlgo == "HyperSingle"){ // Face Positions 2: Hyperbolic Tangent (Single-Sided)
        
        // Mesh Parameters
        double delta = refData["delta"].asDouble(); double A{}, B{};
        if (delta == 0){std::cerr << "Hyperbolic tangent mesh parameter out of range.\n"; std::exit(0);}
        // delta: independent of sign (+/-), higher value pushes towards the left

        // Generate Mesh
        for (int i = cNode; i < cNode+NSec+1; i++){
            A = tanh(delta * ((static_cast<double>(i) - cNode) / NSec - 1)); B = tanh(delta);
            nFaces[iAxis][i] = x0 + length * (1 + A / B);
        }

    } else if(sAlgo == "HyperDouble"){ // Face Positions 3: Hyperbolic Tangent (Double-Sided)
        
        // Mesh Parameters
        double delta = refData["delta"].asDouble(); double A{}, B{};
        if (delta == 0){std::cerr << "Hyperbolic tangent mesh parameter out of range.\n"; std::exit(0);}
        // delta: independent of sign (+/-), higher value pushes towards the edges

        // Generate Mesh
        for (int i = cNode; i < cNode+NSec+1; i++){
            A = tanh(delta*((static_cast<double>(i) - cNode)/NSec - 0.5)); B = tanh(0.5 * delta);
            nFaces[iAxis][i] = x0 + 0.5 * length * (1 + A/B);
        }

    } else if (sAlgo == "Exponential") { // Face Positions 4: Exponential

        // Mesh Parameters
        double kappa = refData["kappa"].asDouble();
        if (kappa == 0){std::cerr << "Exponential mesh parameter out of range.\n"; std::exit(0);}
        // kappa < 0: Shifts right, kappa > 0: Shifts left 
        
        // Generate Mesh
        for (int i = cNode; i < cNode+NSec+1; i++){
            nFaces[iAxis][i] = x0 + length * (exp(kappa * (static_cast<double>(i) - cNode) / NSec) - 1) / (exp(kappa) - 1);
        }

    } 

}


void Mesh::generateMesh(MeshSolver& Msh, double Phi0, Json::Value qNode, Json::Value sections, Json::Value refinement, Json::Value obs){

	// Control (nD)
	Msh.N.resize(qNode.size());
	for(Json::Value::ArrayIndex i = 0; i < Msh.N.size(); i++){
		Msh.N[i] = qNode[i].asInt();
	}
    for (size_t val : Msh.N){Msh.totNodes *= val;}

	// Geometry Resize (nD)
	Msh.Faces.resize(Msh.N.size()); Msh.Nodes.resize(Msh.N.size()); Msh.deltaX.resize(Msh.N.size()); Msh.dX.resize(Msh.N.size()); 
	for (size_t i = 0; i < Msh.N.size(); i++){
		Msh.Faces[i].resize(Msh.N[i]+1); Msh.Nodes[i].resize(Msh.N[i]); Msh.deltaX[i].resize(Msh.N[i]); Msh.dX[i].resize(Msh.N[i]+1); 
	}

    // Faces Loop (nD)
    std::vector<size_t> cNode; cNode.resize(Msh.N.size(), 0);
    for (int i = 0; i < refinement.size(); i++){
        // Calculate
        calculateFaces(cNode, refinement[i], Msh.Faces);

        // Control
        cNode[refinement[i]["axis"].asInt()] += refinement[i]["N"].asInt();
    }

    // CV Positions (nD)
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



    /* // Resize (Non-nD) */
    /* Msh.nMat.resize(Msh.N[0]); Msh.Phi.resize(Msh.N[0]); Msh.sPhi.resize(Msh.N[0]); Msh.Sw.resize(Msh.N[0]); Msh.Se.resize(Msh.N[0]); Msh.Ss.resize(Msh.N[0]); Msh.Sn.resize(Msh.N[0]); Msh.Vp.resize(Msh.N[0]); Msh.oPhi.resize(Msh.N[0]); Msh.bObs.resize(Msh.N[0]); */
    /* for (size_t i = 0; i < Msh.N[0]; i++){ */
        /* Msh.nMat[i].resize(Msh.N[1], 0); Msh.Phi[i].resize(Msh.N[1], Phi0); Msh.sPhi[i].resize(Msh.N[1], 0); Msh.Sw[i].resize(Msh.N[1], 0); Msh.Se[i].resize(Msh.N[1], 0); Msh.Ss[i].resize(Msh.N[1], 0); Msh.Sn[i].resize(Msh.N[1], 0); Msh.Vp[i].resize(Msh.N[1], 0); Msh.oPhi[i].resize(Msh.N[1], Phi0); Msh.bObs[i].resize(Msh.N[1]); */
    /* } */

    /* // Sections */
    /* std::vector<int> Pos0(Msh.N.size()), Pos1(Msh.N.size()); */
    /* for (Json::Value::ArrayIndex i = 0; i < sections.size(); i++){ */
        /* // Find Positions (nD) */
        /* for (int j = 0; j < Msh.N.size(); j++){ */
            /* Pos0[j] = std::lower_bound(Msh.Faces[j].begin(), Msh.Faces[j].end(), sections[i]["x0"][j].asDouble()) - Msh.Faces[j].begin(); */
            /* Pos1[j] = std::lower_bound(Msh.Faces[j].begin(), Msh.Faces[j].end(), sections[i]["x1"][j].asDouble()) - Msh.Faces[j].begin(); */
        /* } */

        /* // Internal Nodes (Non-nD) */
        /* for (int j = Pos0[0]; j < Pos1[0]; j++){ */
            /* for (int k = Pos0[1]; k < Pos1[1]; k++){ */
                /* // Material */
                /* Msh.nMat[j][k] = sections[i]["material"].asInt(); Msh.sPhi[j][k] = sections[i]["source"].asDouble(); */

                /* /1* // Geometry *1/ */
                /* Msh.Sw[j][k] = Msh.deltaX[1][k] * W; Msh.Se[j][k] = Msh.deltaX[1][k] * W; Msh.Ss[j][k] = Msh.deltaX[0][j] * W; Msh.Sn[j][k] = Msh.deltaX[0][j] * W; */
                /* Msh.Vp[j][k] = Msh.deltaX[0][j] * Msh.deltaX[1][k] * W; */
            /* } */
        /* } */
    /* } */

    /* // Obstacles */
    /* Obstacle tempObs{}; tempObs.i0.resize(Msh.N.size()); tempObs.i1.resize(Msh.N.size()); */
    /* for (Json::Value::ArrayIndex i = 0; i < obs.size(); i++){ */
        /* // Find Positions (nD) */
        /* for (int j = 0; j < Msh.N.size(); j++){ */
            /* tempObs.i0[j] = std::lower_bound(Msh.Nodes[j].begin(), Msh.Nodes[j].end(), obs[i]["x0"][j].asDouble()) - Msh.Nodes[j].begin(); */
            /* tempObs.i1[j] = std::lower_bound(Msh.Nodes[j].begin(), Msh.Nodes[j].end(), obs[i]["x1"][j].asDouble()) - Msh.Nodes[j].begin(); */
        /* } */ 

        /* // bObs */
        /* for (size_t j = tempObs.i0[0]; j < tempObs.i1[0]; j++){ */
            /* for (size_t k = tempObs.i0[1]; k < tempObs.i1[1]; k++){ */
                /* Msh.bObs[j][k] = true; */
            /* } */
        /* } */

        /* // Control */
        /* obstacles.push_back(std::move(tempObs)); tempObs = {}; */ 
        /* tempObs.i0.resize(Msh.N.size()); tempObs.i1.resize(Msh.N.size()); */
    /* } */

    /* // Coefficients (nD) */
    /* Msh.matA.resize(Msh.totNodes); Msh.matB.resize(Msh.totNodes, 0); Msh.tempA.resize(Msh.totNodes); Msh.tempB.resize(Msh.totNodes, 0); */

}
