// Imports
#include <cerrno>
#include <iostream>
#include <vector>
#include <fstream>
#include <filesystem>

/* #include "exprtk.hpp" */
#include "o01_Material.h"
#include "o02_Mesh.h"
#include "o03_Discretizer.h"
#include "o05_Probe.h"
#include "o09_Medic.h"

std::ofstream createDiagnostic(std::filesystem::path fName){

    // Open File 
    std::ofstream file(fName);
    if (!file.is_open()){
        std::cerr << "Failed to open file. \n";
    }

    // File Header
    file << "Time";

    return file;

}

Medic::Medic(Mesh Msh, Probe& Prb, bool bExit){

    // Control
    if (!bExit){return;}

    // Path
    pathBase = Prb.pathBase; int k{};
    std::filesystem::path newPath(pathBase);

    // Energy Balance
    file = createDiagnostic(newPath / "Medic_Balance.csv");
    for (size_t i = 0; i < Msh.N[0]; i++){
        for (size_t j = 0; j < Msh.N[1]; j++){
            file << "," << Msh.Nodes[0][i] << " " << Msh.Nodes[1][j];
        }
    } file << "\n";

    // Residue
    fileR = createDiagnostic(newPath / "Medic_Residue.csv");
    for (size_t i = 1; i < Msh.N[0]-1; i++){
        for (size_t j = 1; j < Msh.N[1]-1; j++){
            fileR << "," << Msh.Nodes[0][i] << " " << Msh.Nodes[1][j];
        }
    } fileR << "\n";

}


void checkBoundaries(Material Mat, Mesh Msh, Discretizer Dsc, double t, int xPos, int yPos, int side){

    // Control
    double lamb{};

    // Loop
    for (Boundary bC : Msh.boundaryConditions){

        // NEED TO FIND THE RIGHT BOUNDARY AND IMPOSE THE VALUE OF THE BOUNDARY
        // HOW TO IDENTIFY WHICH SIDE BOUNDARY I AM DOING WITHOUT HAVING TO PUT ALL 12 CASES FOR EACH TYPE OF BOUNDARY CONDITION
        if (xPos >= bC.x0[0] && xPos <=  bC.x1[0] && yPos >= bC.x0[1] && yPos <= bC.x1[1]){

        }

    }

}


void Medic::getDiagnostic(Material Mat, Mesh Msh, Discretizer Dsc, std::vector<std::vector<double>> oldPhi, double t){

    // Control
    double tempErr, impErr, expErr; int k{};
    double lamb{}, lambw{}, lambe{}, lambs{}, lambn{};
    file << t;

    // Nodes Loop
    for (int i = 0; i < Msh.N[0]; i++){
        for (int j = 0; j < Msh.N[1]; j++){
            
            // Control
            k = i * Msh.N[1] + j; impErr = 0; expErr = 0;

            // West Boundary
            if (i != 0){
                impErr += Msh.tempA[k].aw * (Msh.vPhi[i-1][j] - Msh.vPhi[i][j]); 
                expErr += Msh.tempA[k].aw * (oldPhi[i-1][j] - oldPhi[i][j]);
            } else {}

            // East Boundary
            if (i != Msh.N[0]-1){
                impErr += Msh.tempA[k].ae * (Msh.vPhi[i+1][j] - Msh.vPhi[i][j]);
                expErr += Msh.tempA[k].ae * (oldPhi[i+1][j] - oldPhi[i][j]);
            }

            // South Boundary
            if (j != 0){
                impErr += Msh.tempA[k].as * (Msh.vPhi[i][j-1] - Msh.vPhi[i][j]);
                expErr += Msh.tempA[k].as * (oldPhi[i][j-1] - oldPhi[i][j]);
            }

            // North Boundary
            if (j != Msh.N[1]-1){
                impErr += Msh.tempA[k].an * (Msh.vPhi[i][j+1] - Msh.vPhi[i][j]);
                expErr += Msh.tempA[k].an * (oldPhi[i][j+1] - oldPhi[i][j]);
            }

            // Calculate // PENDING - VERIFY EXPRESSION, THIS ONE IGNORES TEMPB
            tempErr = Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp * Msh.Vp[i][j] * (Msh.vPhi[i][j] - oldPhi[i][j]) / Dsc.dt - Dsc.beta * impErr - (1 - Dsc.beta) * expErr + Msh.sPhi[i][j] * Msh.Vp[i][j];

            // Print to File
            file << "," << tempErr;

        }
    } file << "\n";

}


void Medic::getSystemResidual(Material Mat, Mesh Msh, Discretizer Dsc, double t){
    
    // Control
    double tempRes; int k{};
    fileR << t;

    // ALL INTERIOR NODES NOW
    // NEED TO LOOP OVER EVERYTHING AND DO ENERGY BALANCE USING TEMPA AND TEMPB

    // Interior Nodes
    for (int i = 0; i < Msh.N[0]; i++){
        for (int j = 0; j < Msh.N[1]; j++){
            
            // Control
            k = i * Msh.N[1] + j;

            /* std::cout << i << " " << j << " " << k << ": "; */

            // Calculate
            tempRes = 0;
            /* tempRes = Msh.matA[k].ap * Msh.vPhi[i][j] + Msh.matA[k].aw * Msh.vPhi[i-1][j] + Msh.matA[k].ae * Msh.vPhi[i+1][j] + Msh.matA[k].as * Msh.vPhi[i][j-1] + Msh.matA[k].an * Msh.vPhi[i][j+1] - Msh.bp[k]; */

            /* std::cout << tempRes << "\n"; */
            
            /* fileR << "," << tempRes; */

        }
    } // fileR << "\n";
}


void Medic::getGlobalBalance(Material Mat, Mesh Msh, Discretizer Dsc){

    // Internal Heat Generation
    double sumQ = 0;
    for (size_t i = 1; i < Msh.N[0]-1; i++){
        for (size_t j = 1; j < Msh.N[1]-1; j++){
            sumQ += Msh.sPhi[i][j] * Msh.Vp[i][j];
        }
    }

    /* // THIS ONLY TAKES INTO ACCOUNT THE IMPLICIT PART OF THE EXPRESSION, NEED TO INCLUDE THE EXPLICIT PART AS WELL */ 

    /* // Outward flux */
    /* double sumBC = 0; */
    
    /* // xBoundaries (west, east) */
    /* for (size_t i = 1; i < Msh.N[1]-1; i++){ */
    /*     sumBC += Mat.vMat[Msh.nMat[1][i]].gamma * Msh.Sw[1][i] * (Msh.vPhi[1][i] - Msh.vPhi[0][i]) / (Msh.ndelta[0][1] * 0.5); */
    /*     sumBC += Mat.vMat[Msh.nMat[Msh.N[0]-2][i]].gamma * Msh.Se[Msh.N[0]-2][i] * (Msh.vPhi[Msh.N[0]-2][i] - Msh.vPhi[Msh.N[0]-1][i]) / (Msh.ndelta[0][Msh.N[0]-2] * 0.5); */
    /* } */

    /* // yBoundaries (south, north) */
    /* for (size_t i = 1; i < Msh.N[0]-1; i++){ */
    /*     sumBC += Mat.vMat[Msh.nMat[i][1]].gamma * Msh.Ss[i][1] * (Msh.vPhi[i][1] - Msh.vPhi[i][0]) / (Msh.ndelta[1][1] * 0.5); */
    /*     sumBC += Mat.vMat[Msh.nMat[i][Msh.N[1]-2]].gamma * Msh.Sn[i][Msh.N[1]-2] * (Msh.vPhi[i][Msh.N[1]-2] - Msh.vPhi[i][Msh.N[1]-2]) / (Msh.ndelta[1][Msh.N[1]-2] * 0.5); */
    /* } */

    /* std::cout << "Global Energy Balance: " << sumQ << " " << sumBC << " " << sumQ - sumBC << "\n"; */

}


Medic::~Medic(){

	// Close Files
	file.close();
	fileR.close();

}
