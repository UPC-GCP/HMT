// Imports
#include <iostream>
#include <vector>
#include <json/json.h>
#include <cmath>
#include <algorithm>

// Self-Imports
#include "Material.h"
#include "Mesh.h"
#include "ExpressionParser.h"

Mesh::Mesh(int algo, double W, double H, double A, double xC, double kStr, double delta) {

    // Geometry
    this->W = W; this->H = H;

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

std::vector<double> Mesh::splitString(std::string strSplit, char delimiter){
    
    // Control
    std::vector<double> vRes;
    std::stringstream ss(strSplit);
    std::string tempVal;

    while (std::getline(ss, tempVal, delimiter)){
        vRes.push_back(std::stod(tempVal, 0));
    }

    return vRes;

}

void Mesh::newAddBoundaryConditions(Json::Value boundaries, ExpressionParser& Prs){

    // Resize
    newBoundaryConditions.resize(boundaries.size());
    std::vector<double> Pos0, Pos1;

    for (Json::Value::ArrayIndex i = 0; i < boundaries.size(); i++){

        newBoundaryConditions[i].x0.resize(N.size()); newBoundaryConditions[i].x1.resize(N.size());

        if (boundaries[i]["type"] == "Dirichlet") {

            // Control
            newBoundaryConditions[i].type = 0;

            // Position
            for (int j = 0; j < N.size(); j++){
                newBoundaryConditions[i].x0[j] = boundaries[i]["x0"][j].asDouble();
                newBoundaryConditions[i].x1[j] = boundaries[i]["x1"][j].asDouble();
            }

            // Value
            if (isFormula(boundaries[i]["value"].asString())){
                newBoundaryConditions[i].value = 0;
                newBoundaryConditions[i].bUpdate = true;
                newBoundaryConditions[i].iExpr = Prs.registerExpression(boundaries[i]["value"].asString());
            } else {
                newBoundaryConditions[i].value = boundaries[i]["value"].asDouble();
            }

        } else if (boundaries[i]["type"] == "Neumann") {
            
            // Control
            newBoundaryConditions[i].type = 1;

            // Position
            for (int j = 0; j < N.size(); j++){
                newBoundaryConditions[i].x0[j] = boundaries[i]["x0"][j].asDouble();
                newBoundaryConditions[i].x1[j] = boundaries[i]["x1"][j].asDouble();
            }

            // Value
            newBoundaryConditions[i].value = boundaries[i]["value"].asDouble();
            newBoundaryConditions[i].side = boundaries[i]["side"].asInt();

        } else if (boundaries[i]["type"] == "Convection") {
            
            // Control
            newBoundaryConditions[i].type = 2;
            
            // Position
            for (int j = 0; j < N.size(); j++){
                newBoundaryConditions[i].x0[j] = boundaries[i]["x0"][j].asDouble();
                newBoundaryConditions[i].x1[j] = boundaries[i]["x1"][j].asDouble();
            }

            // Value
            newBoundaryConditions[i].value = boundaries[i]["value"].asDouble();
            newBoundaryConditions[i].side = boundaries[i]["side"].asInt();
            newBoundaryConditions[i].alpha = boundaries[i]["alpha"].asDouble();

        } else {
            std::cerr << "Error: Invalid boundary condition type " << boundaries[i]["type"].asString() << std::endl;
        }

    }

}

void Mesh::newCalculateFaces(int cNode, int NSec, double x0, double x1, std::vector<double>& fVec) {

    // General
    double length = x1 - x0;

    // Face Positions
    if (algorithm == 0){
        // Face Positions 0: Bidirectional Non-uniform (A, xC)
        for (int i = cNode; i < cNode+NSec+1; i++) {
            fVec[i] = x0 + (i-cNode) * length / NSec + strength * (centering - (i-cNode) * length / NSec) * (1 - (i-cNode)/NSec) * (i-cNode) / NSec;
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

void Mesh::newGenerateMesh(Material& Mat, Json::Value qNode, Json::Value sections, Json::Value refinement){

    // Control (nD)
    N.resize(qNode.size(), 2);
    for (Json::Value::ArrayIndex i = 0; i < N.size(); i++){
        N[i] = qNode[i].asInt() + 2;
    }
    for (int val : N){totNodes *= val;}

    // Geometry (nD)
    Faces.resize(N.size()); Nodes.resize(N.size()); ndelta.resize(N.size()); nd.resize(N.size());
    for (int i = 0; i < N.size(); i++){
        Faces[i].resize(N[i] - 1); Nodes[i].resize(N[i]); ndelta[i].resize(N[i]); nd[i].resize(N[i]);
    }
    
    // Mesh Loop (nD)
    std::vector<int> cNode; cNode.resize(N.size(), 0);
    for (int i = 0; i < refinement.size(); i++){

        // Faces
        newCalculateFaces(cNode[refinement[i]["axis"].asInt()], refinement[i]["N"].asInt(), refinement[i]["range"][0].asDouble(), refinement[i]["range"][1].asDouble(), Faces[refinement[i]["axis"].asInt()]);

        // Control
        cNode[refinement[i]["axis"].asInt()] += refinement[i]["N"].asInt();

    }
    
    // CV Position (nD)
    for (int i = 0; i < N.size(); i++){
        for (int j = 1; j < N[i]-1; j++){
            Nodes[i][j] = 0.5 * (Faces[i][j] + Faces[i][j-1]);
        }
        Nodes[i].front() = Faces[i].front(); Nodes[i].back() = Faces[i].back();
    }

    // Deltas (nD)
    for (size_t i = 0; i < nd.size(); i++){
        for (size_t j = 0; j < nd[i].size()-1; j++){
            nd[i][j] = Nodes[i][j+1] - Nodes[i][j];
            if (j > 0){
                ndelta[i][j] = Faces[i][j] - Faces[i][j-1];
            }
        }
    }

    // for (size_t i = 0; i < nd.size(); i++){
    //     std::cout << "nDelta " << i << ": " << ndelta[i].size() << " - ";
    //     for (double val : ndelta[i]){
    //         std::cout << val << " ";
    //     } std::cout << "\n";
    // }

    // std::exit(0);

    // Resize
    nMat.resize(N[0]); nQv.resize(N[0]); nT.resize(N[0]); nSw.resize(N[0]); nSe.resize(N[0]); nSs.resize(N[0]); nSn.resize(N[0]); nVp.resize(N[0]);

    // Resize (Non-nD)
    for (size_t i = 0; i < nMat.size(); i++){
        nMat[i].resize(N[1], 0); nQv[i].resize(N[1], 0); nT[i].resize(N[1], Mat.T0); nSw[i].resize(N[1], 0); nSe[i].resize(N[1], 0); nSs[i].resize(N[1], 0); nSn[i].resize(N[1], 0); nVp[i].resize(N[1], 0);
    }

    // Sections Loop
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
                nMat[j][k] = sections[i]["material"].asInt(); nQv[j][k] = sections[i]["qV"].asDouble();

                // Geometry
                nSw[j][k] = ndelta[1][k] * W; nSe[j][k] = ndelta[1][k] * W; nSs[j][k] = ndelta[0][j] * W; nSn[j][k] = ndelta[0][j] * W;
                nVp[j][k] = ndelta[0][j] * ndelta[1][k] * W;
            }
        }

    }

    // Boundary Nodes (Non-nD)
    for (int i = 0; i < nMat.size(); i++){
        nMat[i].front() = nMat[i][1]; nMat[i].back() = nMat[i][nMat[i].size()-2];
        nSe[i].front() = nSe[i][1]; nSw[i].back() = nSw[i][nSw[i].size()-2];
    }
    nMat.front() = nMat[1]; nMat.back() = nMat[nMat.size()-2];
    nSn.front() = nSn[1]; nSs.back() = nSs[nSs.size()-2];

    // Coefficients (nD)
    matA.resize(totNodes); bp.resize(totNodes, 0);

}
