// Imports
#include "json/value.h"
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
/* #include <memory> */
/* #include <pthread.h> */
#include <string>
#include <vector>
#include <json/json.h>
#include <cmath>
#include <algorithm>

// Self-Imports
#include "o01_Material.h"
#include "o02_Mesh.h"
#include "o09_ExpressionParser.h"

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


void Mesh::generateMesh(Material& Mat, Json::Value qNode, Json::Value sections, Json::Value refinement){

	// Control (nD)
	N.resize(qNode.size());
	for(Json::Value::ArrayIndex i = 0; i < N.size(); i++){
		N[i] = qNode[i].asInt();
	}
    for (int val : N){totNodes *= val;}

	// Geometry (nD)
	Faces.resize(N.size()); Nodes.resize(N.size()); DeltaX.resize(N.size()); dX.resize(N.size());
	for (size_t i = 0; i < N.size(); i++){
		Faces[i].resize(N[i]+1); Nodes[i].resize(N[i]); DeltaX[i].resize(N[i]); dX[i].resize(N[i]+1);
	}
    
    // Faces Loop (nD)
    std::vector<size_t> cNode; cNode.resize(N.size(), 0);
    for (int i = 0; i < refinement.size(); i++){
        // Calculate
        calculateFaces(cNode[refinement[i]["axis"].asInt()], refinement[i]["N"].asInt(), refinement[i]["range"][0].asDouble(), refinement[i]["range"][1].asDouble(), Faces[refinement[i]["axis"].asInt()]);

        // Control
        cNode[refinement[i]["axis"].asInt()] += refinement[i]["N"].asInt();
    }

    // CV Positions (nD)
    for (size_t i = 0; i < N.size(); i++){
        for (size_t j = 0; j < N[i]; j++){
            Nodes[i][j] = 0.5 * (Faces[i][j] + Faces[i][j+1]);
        }
    }
    
    // Deltas (nD)
    for (size_t i = 0; i < N.size(); i++){
        for (size_t j = 0; j < DeltaX[i].size(); j++){
            // Delta X
            DeltaX[i][j] = Faces[i][j+1] - Faces[i][j];

            // dX
            if (j == DeltaX[i].size()-1){continue;}
            dX[i][j+1] = Nodes[i][j+1] - Nodes[i][j];
        }

        // dX
        dX[i].front() = Nodes[i].front() - Faces[i].front();
        dX[i].back() = Faces[i].back() - Nodes[i].back();
    }

    // Resize (Non-nD)
    nMat.resize(N[0]); vPhi.resize(N[0]); sPhi.resize(N[0]); Sw.resize(N[0]); Se.resize(N[0]); Ss.resize(N[0]); Sn.resize(N[0]); Vp.resize(N[0]);
    for (size_t i = 0; i < N[0]; i++){
        nMat[i].resize(N[1], 0); vPhi[i].resize(N[1], Mat.Phi0); sPhi[i].resize(N[1], 0); Sw[i].resize(N[1], 0); Se[i].resize(N[1], 0); Ss[i].resize(N[1], 0); Sn[i].resize(N[1], 0); Vp[i].resize(N[1], 0);
    }

    // Sections
    std::vector<int> Pos0(N.size()), Pos1(N.size());
    for (Json::Value::ArrayIndex i = 0; i < sections.size(); i++){
        // Find Positions (nD)
        for (int j = 0; j < N.size(); j++){
            Pos0[j] = std::find(Faces[j].begin(), Faces[j].end(), sections[i]["x0"][j].asDouble()) - Faces[j].begin();
            Pos1[j] = std::find(Faces[j].begin(), Faces[j].end(), sections[i]["x1"][j].asDouble()) - Faces[j].begin();
        }

        // Internal Nodes (Non-nD)
        for (int j = Pos0[0]; j < Pos1[0]; j++){
            for (int k = Pos0[1]; k < Pos1[1]; k++){
                // Material
                nMat[j][k] = sections[i]["material"].asInt(); sPhi[j][k] = sections[i]["source"].asDouble();

                /* // Geometry */
                Sw[j][k] = DeltaX[1][k] * W; Se[j][k] = DeltaX[1][k] * W; Ss[j][k] = DeltaX[0][j] * W; Sn[j][k] = DeltaX[0][j] * W;
                Vp[j][k] = DeltaX[0][j] * DeltaX[1][k] * W;
            }
        }
    }

    // Coefficients (nD)
    matA.resize(totNodes); bp.resize(totNodes, 0); tempA.resize(totNodes); tempB.resize(totNodes, 0);

}

void Mesh::addBoundaryConditions(Json::Value boundaries, ExpressionParser& Prs){

    // Resize
    boundaryConditions.resize(boundaries.size());
    std::vector<double> Pos0, Pos1; std::string sType{};

    for (Json::Value::ArrayIndex i = 0; i < boundaries.size(); i++){

        boundaryConditions[i].x0.resize(N.size()); boundaryConditions[i].x1.resize(N.size());

        if (boundaries[i]["type"] == "Dirichlet") {

            // Control
            boundaryConditions[i].type = 0;
            boundaryConditions[i].side = boundaries[i]["side"].asInt();

            // Position
            for (int j = 0; j < N.size(); j++){
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
                    boundaryConditions[i].iEq = 0; // updateTime
                } else if (sType.find(" x ") != std::string::npos || sType.find(" y ") != std::string::npos){
                    boundaryConditions[i].iEq = 1; // updateCoords
                } else {boundaryConditions[i].iEq = 0;} // Generic uses updateTime

            } else {
                boundaryConditions[i].value = boundaries[i]["value"].asDouble();
            }

        } else if (boundaries[i]["type"] == "Neumann") {
            
            // Control
            boundaryConditions[i].type = 1;

            // Position
            for (int j = 0; j < N.size(); j++){
                boundaryConditions[i].x0[j] = boundaries[i]["x0"][j].asDouble();
                boundaryConditions[i].x1[j] = boundaries[i]["x1"][j].asDouble();
            }

            // Value
            boundaryConditions[i].value = boundaries[i]["value"].asDouble();
            boundaryConditions[i].side = boundaries[i]["side"].asInt();

        } else if (boundaries[i]["type"] == "Robin") {
            
            // Control
            boundaryConditions[i].type = 2;
            
            // Position
            for (int j = 0; j < N.size(); j++){
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


void Mesh::addVelocityField(Json::Value nField, ExpressionParser& Prs){

    // Control
    if (N.size() != nField.size()){std::cerr << "Velocity field does not match dimensions.\n";}

    // Resize
    vConv.resize(N.size()); vExpr.resize(N.size());
    for (size_t k = 0; k < N.size(); k++){
        vConv[k].resize(N[0]+1);
        for (size_t i = 0; i < N[0]+1; i++){
            vConv[k][i].resize(N[1]+1);
        }
    }

    // Parser
    for (Json::Value::ArrayIndex i = 0; i < vConv.size(); i++){
        if (isFormula(nField[i]["value"].asString())){
            vExpr[i] = Prs.registerExpression(nField[i]["value"].asString());
        }
    }
    
    // Calculate Field
    for (size_t k = 0; k < N.size(); k++){
        for (size_t i = 0; i < N[0]+1; i++){
            for (size_t j = 0; j < N[1]+1; j++){
                vConv[k][i][j].Vw = Prs.evaluateCoordinates(vExpr[k], Faces[0][i], Nodes[1][j]);
                vConv[k][i][j].Ve = Prs.evaluateCoordinates(vExpr[k], Faces[0][i+1], Nodes[1][j]);
                vConv[k][i][j].Vs = Prs.evaluateCoordinates(vExpr[k], Nodes[0][i], Faces[1][j]);
                vConv[k][i][j].Vn = Prs.evaluateCoordinates(vExpr[k], Nodes[0][i], Faces[1][j+1]);
            }
        }
    }

}
