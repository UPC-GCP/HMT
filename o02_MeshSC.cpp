// Imports
/* #include <algorithm> */
/* #include <cstddef> */
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
/* #include <iterator> */
#include <string>
#include <utility>
#include <vector>
#include <json/json.h>
#include <cmath>

// Self-Imports
#include "o02_MeshSC.h"


Mesh::Mesh(int algo, double W, double A, double xC, double kStr, double delta) {

    // Geometry
    this->W = W;

    // Mesh Parameters
    algorithm = algo;
    strength = A; centering = xC; kStrength = kStr; this->delta = delta;
    
}


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

void Mesh::calculateFaces(int cNode, int NSec, double x0, double x1, std::vector<double>& fVec) {

    // General
    double length = x1 - x0;

    // Face Positions
    if (algorithm == 0){
        // Face Positions 0: Bidirectional Non-uniform (A, xC)
        for (int i = cNode; i < cNode+NSec+1; i++) {
            // fVec[i] = x0 + (i-cNode) * length / NSec + strength * (centering - (i-cNode) * length / NSec) * (1 - (i-cNode)/NSec) * (i-cNode) / NSec;
		std::cerr << "Bidirectional Non-Uniform not currently working.\n";
        }
    } else if (algorithm == 1){
        // Face Positions 1: Unidirectional Non-uniform (Kappa)
        for (int i = cNode; i < cNode+NSec+1; i++) {
            fVec[i] = x0 + pow(((i-cNode) * length / NSec), kStrength);
        }
    } else if(algorithm == 2){
        // Face Positions 2: Hyperbolic Tangent (Single Side)
        double A, B;
        for (int i = cNode; i < cNode+NSec+1; i++){
            A = tanh(delta * ((static_cast<double>(i) - cNode) / NSec - 1)); B = tanh(delta);
            fVec[i] = x0 + length * (1 + A / B);
        }
    } else if(algorithm == 3){
        // Face Positions 3: Hyperbolic Tangent (Double-Sided)
        double A, B;
        for (int i = cNode; i < cNode+NSec+1; i++){
            A = tanh(delta*((static_cast<double>(i) - cNode)/NSec - 0.5));
            B = tanh(0.5 * delta);
            fVec[i] = x0 + 0.5 * length * (1 + A/B);
        }
    }

}


void Mesh::generateMesh(MeshSolver& Msh, double Phi0, Json::Value qNode, Json::Value sections, Json::Value refinement, Json::Value obs){

	// Control (nD)
	Msh.N.resize(qNode.size());
	for(Json::Value::ArrayIndex i = 0; i < Msh.N.size(); i++){
		Msh.N[i] = qNode[i].asInt();
	}
    for (int val : Msh.N){Msh.totNodes *= val;}

    std::cout << "Test1\n";

	// Geometry Resize (nD)
	Msh.Faces.resize(Msh.N.size()); Msh.Nodes.resize(Msh.N.size()); Msh.deltaX.resize(Msh.N.size()); Msh.dX.resize(Msh.N.size()); Msh.bObs.resize(Msh.N.size());

    std::cout << "Test2\n";

	for (size_t i = 0; i < Msh.N.size(); i++){
		Msh.Faces[i].resize(Msh.N[i]+1); Msh.Nodes[i].resize(Msh.N[i]); Msh.deltaX[i].resize(Msh.N[i]); Msh.dX[i].resize(Msh.N[i]+1); Msh.bObs[i].resize(Msh.N[i], false);
	}
    
    std::cout << "Test3\n";
    // Faces Loop (nD)
    std::vector<size_t> cNode; cNode.resize(Msh.N.size(), 0);
    for (int i = 0; i < refinement.size(); i++){
        // Calculate
        calculateFaces(cNode[refinement[i]["axis"].asInt()], refinement[i]["N"].asInt(), refinement[i]["range"][0].asDouble(), refinement[i]["range"][1].asDouble(), Msh.Faces[refinement[i]["axis"].asInt()]);

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

    // Resize (Non-nD)
    Msh.nMat.resize(Msh.N[0]); Msh.Phi.resize(Msh.N[0]); Msh.sPhi.resize(Msh.N[0]); Msh.Sw.resize(Msh.N[0]); Msh.Se.resize(Msh.N[0]); Msh.Ss.resize(Msh.N[0]); Msh.Sn.resize(Msh.N[0]); Msh.Vp.resize(Msh.N[0]); Msh.oPhi.resize(Msh.N[0]);
    for (size_t i = 0; i < Msh.N[0]; i++){
        Msh.nMat[i].resize(Msh.N[1], 0); Msh.Phi[i].resize(Msh.N[1], Phi0); Msh.sPhi[i].resize(Msh.N[1], 0); Msh.Sw[i].resize(Msh.N[1], 0); Msh.Se[i].resize(Msh.N[1], 0); Msh.Ss[i].resize(Msh.N[1], 0); Msh.Sn[i].resize(Msh.N[1], 0); Msh.Vp[i].resize(Msh.N[1], 0); Msh.oPhi[i].resize(Msh.N[1], Phi0);
    }

    std::cout << "Test4\n";

    // Sections
    std::vector<int> Pos0(Msh.N.size()), Pos1(Msh.N.size());

    std::cout << "Test5\n";

    for (Json::Value::ArrayIndex i = 0; i < sections.size(); i++){
        
        std::cout << i << ": " << Msh.nMat.size() << " " << Msh.nMat[0].size() << "\n";

        // Find Positions (nD)
        for (int j = 0; j < Msh.N.size(); j++){
            Pos0[j] = std::find(Msh.Faces[j].begin(), Msh.Faces[j].end(), sections[i]["x0"][j].asDouble()) - Msh.Faces[j].begin();
            Pos1[j] = std::find(Msh.Faces[j].begin(), Msh.Faces[j].end(), sections[i]["x1"][j].asDouble()) - Msh.Faces[j].begin();
        }

        std::cout << i << ": [" << Pos0[0] << "," << Pos0[1] << "] [" << Pos1[0] << "," << Pos1[1] << "]\n";

        // Internal Nodes (Non-nD)
        for (int j = Pos0[0]; j < Pos1[0]; j++){
            for (int k = Pos0[1]; k < Pos1[1]; k++){
                // Material
                Msh.nMat[j][k] = sections[i]["material"].asInt(); Msh.sPhi[j][k] = sections[i]["source"].asDouble();

                /* // Geometry */
                Msh.Sw[j][k] = Msh.deltaX[1][k] * W; Msh.Se[j][k] = Msh.deltaX[1][k] * W; Msh.Ss[j][k] = Msh.deltaX[0][j] * W; Msh.Sn[j][k] = Msh.deltaX[0][j] * W;
                Msh.Vp[j][k] = Msh.deltaX[0][j] * Msh.deltaX[1][k] * W;
            }
        }

        std::cout << "Values initialized with no issues\n";

    }

    std::cout << "Test6\n";



    // Obstacles
    Obstacle tempObs{}; 
    for (Json::Value::ArrayIndex i = 0; i < obs.size(); i++){
        // This needs to find positions and store them in obstacles
        // bObs = true for all nodes in range
        
        // Find Positions (nD)
        for (int j = 0; j < Msh.N.size(); j++){
            tempObs.i0[j] = std::lower_bound(Msh.Nodes[j].begin(), Msh.Nodes[j].end(), obs[i]["x0"][j].asDouble()) - Msh.Nodes[j].begin();
            tempObs.i1[j] = std::lower_bound(Msh.Nodes[j].begin(), Msh.Nodes[j].end(), obs[i]["x1"][j].asDouble()) - Msh.Nodes[j].begin();
        } obstacles.push_back(std::move(tempObs)); tempObs = {};

        // bObs



    }





    // Coefficients (nD)
    Msh.matA.resize(Msh.totNodes); Msh.matB.resize(Msh.totNodes, 0); Msh.tempA.resize(Msh.totNodes); Msh.tempB.resize(Msh.totNodes, 0);

}
