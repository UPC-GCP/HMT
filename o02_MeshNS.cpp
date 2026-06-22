// Imports
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
/* #include <memory> */
/* #include <pthread.h> */
#include <iterator>
#include <string>
#include <vector>
#include <json/json.h>
#include <cmath>

// Self-Imports
/* #include "o01_Material.h" */
#include "o02_MeshNS.h"
/* #include "o09_ExpressionParser.h" */


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

void Mesh::generateMesh(MeshSolver& Msh, Material Mat, Json::Value qNode, Json::Value sections, Json::Value refinement){

	// Control (nD)
	Msh.N.resize(qNode.size());
	for(Json::Value::ArrayIndex i = 0; i < Msh.N.size(); i++){
		Msh.N[i] = qNode[i].asInt();
	}
    for (int val : Msh.N){Msh.totNodes *= val;}

	// Geometry Resize (nD)
	Msh.Faces.resize(Msh.N.size()); Msh.Nodes.resize(Msh.N.size()); Msh.deltaX.resize(Msh.N.size()); Msh.dX.resize(Msh.N.size());
	for (size_t i = 0; i < Msh.N.size(); i++){
		Msh.Faces[i].resize(Msh.N[i]+1); Msh.Nodes[i].resize(Msh.N[i]); Msh.deltaX[i].resize(Msh.N[i]); Msh.dX[i].resize(Msh.N[i]+1);
	}
    
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
        Msh.nMat[i].resize(Msh.N[1], 0); Msh.Phi[i].resize(Msh.N[1], Mat.Phi0); Msh.sPhi[i].resize(Msh.N[1], 0); Msh.Sw[i].resize(Msh.N[1], 0); Msh.Se[i].resize(Msh.N[1], 0); Msh.Ss[i].resize(Msh.N[1], 0); Msh.Sn[i].resize(Msh.N[1], 0); Msh.Vp[i].resize(Msh.N[1], 0); Msh.oPhi[i].resize(Msh.N[1], Mat.Phi0);
    }

    // Sections
    std::vector<int> Pos0(Msh.N.size()), Pos1(Msh.N.size());
    for (Json::Value::ArrayIndex i = 0; i < sections.size(); i++){
        // Find Positions (nD)
        for (int j = 0; j < Msh.N.size(); j++){
            Pos0[j] = std::find(Msh.Faces[j].begin(), Msh.Faces[j].end(), sections[i]["x0"][j].asDouble()) - Msh.Faces[j].begin();
            Pos1[j] = std::find(Msh.Faces[j].begin(), Msh.Faces[j].end(), sections[i]["x1"][j].asDouble()) - Msh.Faces[j].begin();
        }

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
    }

    // Coefficients (nD)
    Msh.matA.resize(Msh.totNodes); Msh.matB.resize(Msh.totNodes, 0); Msh.tempA.resize(Msh.totNodes); Msh.tempB.resize(Msh.totNodes, 0);

}



void Mesh::calculateMeshGeometry(MeshBase& Msh, double valInit){

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

    /* Msh.sPhi.resize(Msh.N[0]); */  // missing for MeshSolver
    /* Msh.sPhi[i].resize(Msh.N[1], 0); */ // inside the loop
    // Resize (Non-nD)
    Msh.Phi.resize(Msh.N[0]); Msh.Sw.resize(Msh.N[0]); Msh.Se.resize(Msh.N[0]); Msh.Ss.resize(Msh.N[0]); Msh.Sn.resize(Msh.N[0]); Msh.Vp.resize(Msh.N[0]); Msh.oPhi.resize(Msh.N[0]);
    for (size_t i = 0; i < Msh.N[0]; i++){
        Msh.Phi[i].resize(Msh.N[1], valInit); Msh.Sw[i].resize(Msh.N[1], 0); Msh.Se[i].resize(Msh.N[1], 0); Msh.Ss[i].resize(Msh.N[1], 0); Msh.Sn[i].resize(Msh.N[1], 0); Msh.Vp[i].resize(Msh.N[1], 0); Msh.oPhi[i].resize(Msh.N[1], valInit);
    }

    // Geometry
    for (size_t j = 0; j < Msh.N[0]; j++){
        for (size_t k = 0; k < Msh.N[1]; k++){
                Msh.Sw[j][k] = Msh.deltaX[1][k] * W; Msh.Se[j][k] = Msh.deltaX[1][k] * W; Msh.Ss[j][k] = Msh.deltaX[0][j] * W; Msh.Sn[j][k] = Msh.deltaX[0][j] * W;
                Msh.Vp[j][k] = Msh.deltaX[0][j] * Msh.deltaX[1][k] * W;
        }
    }

    // Coefficients (nD)
    Msh.matA.resize(Msh.totNodes); Msh.matB.resize(Msh.totNodes, 0); 
    /* Msh.tempA.resize(Msh.totNodes); Msh.tempB.resize(Msh.totNodes, 0); */ 

}


void Mesh::generateMeshVelocity(Material Mat, MeshSolver p, MeshBase& u, MeshBase& v){
    
	// Control (nD)
	u.N.resize(p.N.size()); v.N.resize(p.N.size());
    u.N[0] = p.N[0] + 1; u.N[1] = p.N[1]; for (int val : u.N){u.totNodes *= val;}
    v.N[0] = p.N[0]; v.N[1] = p.N[1] + 1; for (int val : v.N){v.totNodes *= val;}
    
    // Geometry Resize (nD)
    u.Faces.resize(p.N.size()); u.Nodes.resize(p.N.size()); u.deltaX.resize(p.N.size()); u.dX.resize(p.N.size());
    v.Faces.resize(p.N.size()); v.Nodes.resize(p.N.size()); v.deltaX.resize(p.N.size()); v.dX.resize(p.N.size());
	for (size_t i = 0; i < p.N.size(); i++){
		u.Faces[i].resize(u.N[i]+1); u.Nodes[i].resize(u.N[i]); u.deltaX[i].resize(u.N[i]); u.dX[i].resize(u.N[i]);
        v.Faces[i].resize(v.N[i]+1); v.Nodes[i].resize(v.N[i]); v.deltaX[i].resize(v.N[i]); v.dX[i].resize(v.N[i]);
	}

    // Geometry uNodes (non-nD)
    for (size_t i = 0; i < u.Nodes[0].size(); i++){u.Nodes[0][i] = p.Faces[0][i]; u.Faces[0][i+1] = p.Nodes[0][i];} 
    u.Faces[0][0] = u.Nodes[0].front() - u.Faces[0][1]; u.Faces[0].back() = u.Nodes[0].back() + u.Faces[0][1];
    for (size_t i = 0; i < u.Nodes[1].size(); i++){u.Nodes[1][i] = p.Nodes[1][i]; u.Faces[1][i] = p.Faces[1][i];}
    u.Faces[1].back() = p.Faces[1].back();

    // Geometry vNodes (non-nD)
    for (size_t i = 0; i < v.Nodes[0].size(); i++){v.Nodes[0][i] = p.Nodes[0][i]; v.Faces[0][i] = p.Faces[0][i];}
    v.Faces[0].back() = p.Faces[0].back();
    for (size_t i = 0; i < v.Nodes[1].size(); i++){v.Nodes[1][i] = p.Faces[1][i]; v.Faces[1][i+1] = p.Nodes[1][i];}
    v.Faces[1][0] = v.Nodes[1].front() - v.Faces[1][1]; v.Faces[1].back() = v.Nodes[1].back() + v.Faces[1][1];

    // Calculate Geometry
    calculateMeshGeometry(u, Mat.VF0[0]);
    calculateMeshGeometry(v, Mat.VF0[1]);

}


void Mesh::addBoundaryConditions(Json::Value boundaries, Material Mat, ExpressionParser& Prs){

    // Resize
    std::string sType{};
    boundaryConditions.resize(boundaries.size());

    for (Json::Value::ArrayIndex i = 0; i < boundaries.size(); i++){

        // Control
        boundaryConditions[i].x0.resize(p.N.size()); boundaryConditions[i].x1.resize(p.N.size());

        if (boundaries[i]["type"] == "Dirichlet") {

            // Control
            boundaryConditions[i].type = 0;
            boundaryConditions[i].side = boundaries[i]["side"].asInt();

            // Position
            for (int j = 0; j < p.N.size(); j++){
                boundaryConditions[i].x0[j] = boundaries[i]["x0"][j].asDouble();
                boundaryConditions[i].x1[j] = boundaries[i]["x1"][j].asDouble();
            }

            // Value
            if (isFormula(boundaries[i]["value"].asString())){
                
                boundaryConditions[i].value = 0;
                boundaryConditions[i].bUpdate = true;
                boundaryConditions[i].iExpr = Prs.registerExpression(boundaries[i]["value"].asString()); // FutureWork: Store expressions as strings and check if same expression is already stored to reuse the index instead of saving it once more
                
                sType = boundaries[i]["value"].asString();
                if (sType.find(" t ") != std::string::npos){
                    boundaryConditions[i].iEq = 0; boundaryConditions[i].value = Prs.evaluateTime(boundaryConditions[i].iEq, 0); // updateTime
                } else if (sType.find(" x ") != std::string::npos || sType.find(" y ") != std::string::npos){
                    boundaryConditions[i].iEq = 1; // updateCoords
                } else {boundaryConditions[i].iEq = 0;} // Generic uses updateTime

            } else {
                boundaryConditions[i].value = boundaries[i]["value"].asDouble();
            }

            // Store Values
            if (boundaries[i]["x0"][0].asDouble() == boundaries[i]["x1"][0].asDouble()){
                boundaryConditions[i].Phi.resize(p.N[1], boundaryConditions[i].value);
                boundaryConditions[i].oPhi.resize(p.N[1], boundaryConditions[i].value); 
            } else if (boundaries[i]["x0"][1].asDouble() == boundaries[i]["x1"][1].asDouble()){
                boundaryConditions[i].Phi.resize(p.N[0], boundaryConditions[i].value);
                boundaryConditions[i].oPhi.resize(p.N[0], boundaryConditions[i].value);
            }

        } else if (boundaries[i]["type"] == "Neumann") {
            
            // Control
            boundaryConditions[i].type = 1;

            // Position
            for (int j = 0; j < p.N.size(); j++){
                boundaryConditions[i].x0[j] = boundaries[i]["x0"][j].asDouble();
                boundaryConditions[i].x1[j] = boundaries[i]["x1"][j].asDouble();
            }

            // Value
            boundaryConditions[i].value = boundaries[i]["value"].asDouble();
            boundaryConditions[i].side = boundaries[i]["side"].asInt();

            // Store Values
            if (boundaries[i]["x0"][0].asDouble() == boundaries[i]["x1"][0].asDouble()){
                boundaryConditions[i].Phi.resize(p.N[1], boundaryConditions[i].value);
                boundaryConditions[i].oPhi.resize(p.N[1], boundaryConditions[i].value);
            } else if (boundaries[i]["x0"][1].asDouble() == boundaries[i]["x1"][1].asDouble()){
                boundaryConditions[i].Phi.resize(p.N[0], boundaryConditions[i].value);
                boundaryConditions[i].oPhi.resize(p.N[0], boundaryConditions[i].value);
            }

        } else if (boundaries[i]["type"] == "Robin") {
            
            // Control
            boundaryConditions[i].type = 2;
            
            // Position
            for (int j = 0; j < p.N.size(); j++){
                boundaryConditions[i].x0[j] = boundaries[i]["x0"][j].asDouble();
                boundaryConditions[i].x1[j] = boundaries[i]["x1"][j].asDouble();
            }

            // Value
            boundaryConditions[i].value = boundaries[i]["value"].asDouble();
            boundaryConditions[i].side = boundaries[i]["side"].asInt();
            boundaryConditions[i].alpha = boundaries[i]["alpha"].asDouble();

        } else {
            std::cerr << "Error: Invalid boundary condition type " << boundaries[i]["type"].asString() << std::endl;
        }

    }

}


void Mesh::addBoundariesVelocity(Json::Value boundaries, Material Mat){

    // Resize
    boundaryVelocity.resize(boundaries.size());
    std::vector<int> Pos0{}, Pos1{}; Pos0.resize(p.N.size()); Pos1.resize(p.N.size());

    for (Json::Value::ArrayIndex i = 0; i < boundaries.size(); i++){

        // Control
        boundaryVelocity[i].i0.resize(p.N.size()); boundaryVelocity[i].i1.resize(p.N.size());

        if (boundaries[i]["type"] == "Dirichlet") {

            // Control
            boundaryVelocity[i].type = 0;
            boundaryVelocity[i].side = boundaries[i]["side"].asInt();
            boundaryVelocity[i].uVal = boundaries[i]["uValue"].asDouble();
            boundaryVelocity[i].vVal = boundaries[i]["vValue"].asDouble();
            
            // Positions
            Pos0[0] = boundaries[i]["x0"][0].asInt(); Pos0[1] = boundaries[i]["x0"][1].asInt();
            Pos1[0] = boundaries[i]["x1"][0].asInt(); Pos1[1] = boundaries[i]["x1"][1].asInt();

            // xBoundary: xIndex = u.Nodes[0], yCoord = u.Faces[1]
            // yBoundary: yIndex = v.Nodes[1], xCoord = v.Faces[0]

            // Indexes    
            if (Pos0[0] == Pos1[0]){

                // First 
                boundaryVelocity[i].i0[0] = std::distance(u.Nodes[0].begin(), std::find(u.Nodes[0].begin(), u.Nodes[0].end(), boundaries[i]["x0"][0].asDouble())); 
                boundaryVelocity[i].i0[1] = std::distance(u.Faces[1].begin(), std::find(u.Faces[1].begin(), u.Faces[1].end(), boundaries[i]["x0"][1].asDouble()));

                // Last
                boundaryVelocity[i].i1[0] = std::distance(u.Nodes[0].begin(), std::find(u.Nodes[0].begin(), u.Nodes[0].end(), boundaries[i]["x1"][0].asDouble()));
                boundaryVelocity[i].i1[1] = std::distance(u.Faces[1].begin(), std::find(u.Faces[1].begin(), u.Faces[1].end(), boundaries[i]["x1"][1].asDouble())) - 1;

            } else if (Pos0[1] == Pos1[1]){

                // First 
                boundaryVelocity[i].i0[0] = std::distance(v.Faces[0].begin(), std::find(v.Faces[0].begin(), v.Faces[0].end(), boundaries[i]["x0"][0].asDouble()));
                boundaryVelocity[i].i0[1] = std::distance(v.Nodes[1].begin(), std::find(v.Nodes[1].begin(), v.Nodes[1].end(), boundaries[i]["x1"][1].asDouble()));

                // Last
                boundaryVelocity[i].i1[0] = std::distance(v.Faces[0].begin(), std::find(v.Faces[0].begin(), v.Faces[0].end(), boundaries[i]["x1"][0].asDouble())) - 1;
                boundaryVelocity[i].i1[1] = std::distance(v.Nodes[1].begin(), std::find(v.Nodes[1].begin(), v.Nodes[1].end(), boundaries[i]["x1"][1].asDouble()));

            } else {std::cerr << "Boundary range not specified correctly.\n";}

        } else if (boundaries[i]["type"] == "Neumann"){

            // Control
            boundaryVelocity[i].type = 1;
            boundaryVelocity[i].side = boundaries[i]["side"].asInt();
            boundaryVelocity[i].uVal = boundaries[i]["uValue"].asDouble();
            boundaryVelocity[i].vVal = boundaries[i]["vValue"].asDouble();

            // Positions
            // PENDING IMPLEMENT LATER

        } else {std::cerr << "Velocity only accepts Dirichlet at current build.\n";}

    }

}

