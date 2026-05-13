// Imports
#include "json/value.h"
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <pthread.h>
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


bool Mesh::isFormula(std::string value){

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
	N.resize(qNode.size(), 2);
	for(Json::Value::ArrayIndex i = 0; i < N.size(); i++){
		N[i] = qNode[i].asInt() + 2;
	}
    for (int val : N){totNodes *= val;}

	// Geometry (nD)
	Faces.resize(N.size()); Nodes.resize(N.size()); ndelta.resize(N.size()); nd.resize(N.size());
	for (size_t i = 0; i < N.size(); i++){
		Faces[i].resize(N[i]-1); Nodes[i].resize(N[i]); ndelta[i].resize(N[i]); nd[i].resize(N[i]);
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
        for (size_t j = 1; j < N[i]-1; j++){
            Nodes[i][j] = 0.5 * (Faces[i][j] + Faces[i][j-1]);
        }
        Nodes[i].front() = Faces[i].front(); Nodes[i].back() = Faces[i].back();
    }

    // Deltas (nD)
    for (size_t i = 0; i < N.size(); i++){
        for (size_t j = 0; j < nd[i].size()-1; j++){
            nd[i][j] = Nodes[i][j+1] - Nodes[i][j];
            if (j > 0){ndelta[i][j] = Faces[i][j] - Faces[i][j-1];}
        }
    }


    /* std::cout << "deltaX:\n"; */
    /* for (std::vector<double> vec : ndelta){ */
    /*     for (double val : vec){std::cout << val << " ";} std::cout << "\n"; */
    /* } */

    /* std::cout << "dX:\n"; */
    /* for (std::vector<double> vec : nd){ */
    /*     for (double val : vec){std::cout << val << " ";} std::cout << "\n"; */
    /* } */

    // Resize (Non-nD)
    nMat.resize(N[0]); sPhi.resize(N[0]); vPhi.resize(N[0]); Sw.resize(N[0]); Se.resize(N[0]); Ss.resize(N[0]); Sn.resize(N[0]); Vp.resize(N[0]);
    for (size_t i = 0; i < N[0]; i++){
        nMat[i].resize(N[1], 0); sPhi[i].resize(N[1], 0); vPhi[i].resize(N[1], 0); Sw[i].resize(N[1], 0); Se[i].resize(N[1], 0); Ss[i].resize(N[1], 0); Sn[i].resize(N[1], 0); Vp[i].resize(N[1], 0);
    }

    // Sections
    std::vector<int> xPos(N.size()), yPos(N.size());
    for (Json::Value::ArrayIndex i = 0; i < sections.size(); i++){

        // Find Positions (nD)
        for (int j = 0; j < N.size(); j++){
            xPos[j] = std::find(Faces[j].begin(), Faces[j].end(), sections[i]["x0"][j].asDouble()) - Faces[j].begin();
            yPos[j] = std::find(Faces[j].begin(), Faces[j].end(), sections[i]["x1"][j].asDouble()) - Faces[j].begin();
        }

        // Internal Nodes (Non-nD)
        for (int j = xPos[0]+1; j < yPos[0]+1; j++){
            for (int k = xPos[1]+1; k < yPos[1]+1; k++){
                // Material
                nMat[j][k] = sections[i]["material"].asInt(); sPhi[j][k] = sections[i]["source"].asDouble();

                // Geometry
                Sw[j][k] = ndelta[1][k] * W; Se[j][k] = ndelta[1][k] * W; Ss[j][k] = ndelta[0][j] * W; Sn[j][k] = ndelta[0][j] * W;
                Vp[j][k] = ndelta[0][j] * ndelta[1][k] * W;
            }
        }

    }

    // Boundary Nodes (Non-nD)
    for (int i = 0; i < nMat.size(); i++){
        nMat[i].front() = nMat[i][1]; nMat[i].back() = nMat[i][nMat[i].size()-2];
        Se[i].front() = Se[i][1]; Sw[i].back() = Sw[i][Sw[i].size()-2];
    }
    nMat.front() = nMat[1]; nMat.back() = nMat[nMat.size()-2];
    Sn.front() = Sn[1]; Ss.back() = Ss[Ss.size()-2];

    // Coefficients (nD)
    matA.resize(totNodes); bp.resize(totNodes, 0);


/*     std::cout << "Sw:\n"; */
/*     for (std::vector<double> vec : Sw){ */
/*         for (double val : vec){std::cout << val << " ";} std::cout << "\n"; */
/*     } */

/*     std::cout << "Se:\n"; */
/*     for (std::vector<double> vec : Se){ */
/*         for (double val : vec){std::cout << val << " ";} std::cout << "\n"; */
/*     } */

/*     std::cout << "Ss:\n"; */
/*     for (std::vector<double> vec : Ss){ */
/*         for (double val : vec){std::cout << val << " ";} std::cout << "\n"; */
/*     } */

/*     std::cout << "Sn:\n"; */
/*     for (std::vector<double> vec : Sn){ */
/*         for (double val : vec){std::cout << val << " ";} std::cout << "\n"; */
/*     } */

}

void Mesh::addBoundaryConditions(Json::Value boundaries, ExpressionParser& Prs){

    // Resize
    boundaryConditions.resize(boundaries.size());
    std::vector<double> Pos0, Pos1;

    for (Json::Value::ArrayIndex i = 0; i < boundaries.size(); i++){

        boundaryConditions[i].x0.resize(N.size()); boundaryConditions[i].x1.resize(N.size());

        if (boundaries[i]["type"] == "Dirichlet") {

            // Control
            boundaryConditions[i].type = 0;

            // Position
            for (int j = 0; j < N.size(); j++){
                boundaryConditions[i].x0[j] = boundaries[i]["x0"][j].asDouble();
                boundaryConditions[i].x1[j] = boundaries[i]["x1"][j].asDouble();
            }

            // Value
            if (isFormula(boundaries[i]["value"].asString())){
                boundaryConditions[i].value = 0;
                boundaryConditions[i].bUpdate = true;
                boundaryConditions[i].iExpr = Prs.registerExpression(boundaries[i]["value"].asString());
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

        } else if (boundaries[i]["type"] == "Convection") {
            
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
    vField.vConv.resize(N.size()); vField.iExpr.resize(N.size());
    for (size_t i = 0; i < N.size(); i++){
        vField.vConv[i].resize(N[0]);
        for (size_t j = 0; j < N[0]; j++){
            vField.vConv[i][j].resize(N[1], 0);
        }
    }
    
    // Parser
    for (Json::Value::ArrayIndex i = 0; i < vField.vConv.size(); i++){
        if (isFormula(nField[i]["value"].asString())){
            vField.iExpr[i] = Prs.registerExpression(nField[i]["value"].asString());
        }
    }

    // Calculate Field
    for (size_t k = 0; k < N.size(); k++){
        for (size_t i = 0; i < N[0]; i++){
            for (size_t j = 0; j < N[1]; j++){
                vField.vConv[k][i][j] = Prs.evaluateCoordinates(vField.iExpr[k], Nodes[0][i] , Nodes[1][j]);
            }
        }
    }

}
