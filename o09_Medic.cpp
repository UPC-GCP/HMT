// Imports
#include <cerrno>
#include <iostream>
#include <vector>
#include <fstream>
#include <filesystem>

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

Medic::Medic(Mesh Msh, Probe& Prb){

    // Path
    pathBase = Prb.pathBase; int k{};
    std::filesystem::path newPath(pathBase);

    // Energy Balance
    file = createDiagnostic(newPath / "Medic_Balance.csv");
    for (size_t i = 1; i < Msh.N[0]-1; i++){
        for (size_t j = 1; j < Msh.N[1]-1; j++){
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

void Medic::getDiagnostic(Material Mat, Mesh Msh, Discretizer Dsc, std::vector<std::vector<double>> oldTemp, double t){

    // Control
    double tempErr, impErr, expErr; int k{};
    double lamb{}, lambw{}, lambe{}, lambs{}, lambn{};
    file << t;

    // Nodes Loop
    for (int i = 1; i < Msh.N[0]-1; i++){
        for (int j = 1; j < Msh.N[1]-1; j++){
            
            // Control
            k = i * Msh.N[1] + j;

            // Harmonic Mean
            lamb = Mat.vMat[Msh.nMat[i][j]].gamma;
            lambw = Dsc.calcHarmonicMean(Msh.nd[0][i-1], {lamb, Mat.vMat[Msh.nMat[i-1][j]].gamma}, {Msh.ndelta[0][i], Msh.ndelta[0][i-1]});
            lambe = Dsc.calcHarmonicMean(Msh.nd[0][i], {lamb, Mat.vMat[Msh.nMat[i+1][j]].gamma}, {Msh.ndelta[0][i], Msh.ndelta[0][i+1]});
            lambs = Dsc.calcHarmonicMean(Msh.nd[1][j-1], {lamb, Mat.vMat[Msh.nMat[i][j-1]].gamma}, {Msh.ndelta[1][j], Msh.ndelta[1][j-1]});
            lambn = Dsc.calcHarmonicMean(Msh.nd[1][j], {lamb, Mat.vMat[Msh.nMat[i][j+1]].gamma}, {Msh.ndelta[1][j], Msh.ndelta[1][j+1]});

            // Calculate Error
	    impErr = - lambw * Msh.Sw[i][j] * (Msh.vPhi[i][j] - Msh.vPhi[i-1][j]) / Msh.nd[0][i-1] - lambe * Msh.Se[i][j] * (Msh.vPhi[i][j] - Msh.vPhi[i+1][j]) / Msh.nd[0][i] - lambs * Msh.Ss[i][j] * (Msh.vPhi[i][j] - Msh.vPhi[i][j-1]) / Msh.nd[1][j-1] - lambn * Msh.Sn[i][j] * (Msh.vPhi[i][j] - Msh.vPhi[i][j+1]) / Msh.nd[1][j];
	    expErr = - lambw * Msh.Sw[i][j] * (oldTemp[i][j] - oldTemp[i-1][j]) / Msh.nd[0][i-1] - lambe * Msh.Se[i][j] * (oldTemp[i][j] - oldTemp[i+1][j]) / Msh.nd[0][i] - lambs * Msh.Ss[i][j] * (oldTemp[i][j] - oldTemp[i][j-1]) / Msh.nd[1][j-1] - lambn * Msh.Sn[i][j] * (oldTemp[i][j] - oldTemp[i][j+1]) / Msh.nd[1][j];
            tempErr = Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp * Msh.Vp[i][j] * (Msh.vPhi[i][j] - oldTemp[i][j]) / Dsc.dt - Dsc.beta * impErr - (1 - Dsc.beta) * expErr - Msh.sPhi[i][j] * Msh.Vp[i][j]; 

            // Print to File
            file << "," << tempErr;

        }
    } file << "\n";

}


void Medic::getGlobalBalance(Material Mat, Mesh Msh, Discretizer Dsc){

    // Internal Heat Generation
    double sumQ = 0;
    for (size_t i = 1; i < Msh.N[0]-1; i++){
        for (size_t j = 1; j < Msh.N[1]-1; j++){
            sumQ += Msh.sPhi[i][j] * Msh.Vp[i][j];
        }
    }

    // THIS ONLY TAKES INTO ACCOUNT THE IMPLICIT PART OF THE EXPRESSION, NEED TO INCLUDE THE EXPLICIT PART AS WELL 

    // Outward flux
    double sumBC = 0;
    
    // xBoundaries (west, east)
    for (size_t i = 1; i < Msh.N[1]-1; i++){
        sumBC += Mat.vMat[Msh.nMat[1][i]].gamma * Msh.Sw[1][i] * (Msh.vPhi[1][i] - Msh.vPhi[0][i]) / (Msh.ndelta[0][1] * 0.5);
        sumBC += Mat.vMat[Msh.nMat[Msh.N[0]-2][i]].gamma * Msh.Se[Msh.N[0]-2][i] * (Msh.vPhi[Msh.N[0]-2][i] - Msh.vPhi[Msh.N[0]-1][i]) / (Msh.ndelta[0][Msh.N[0]-2] * 0.5);
    }

    // yBoundaries (south, north)
    for (size_t i = 1; i < Msh.N[0]-1; i++){
        sumBC += Mat.vMat[Msh.nMat[i][1]].gamma * Msh.Ss[i][1] * (Msh.vPhi[i][1] - Msh.vPhi[i][0]) / (Msh.ndelta[1][1] * 0.5);
        sumBC += Mat.vMat[Msh.nMat[i][Msh.N[1]-2]].gamma * Msh.Sn[i][Msh.N[1]-2] * (Msh.vPhi[i][Msh.N[1]-2] - Msh.vPhi[i][Msh.N[1]-2]) / (Msh.ndelta[1][Msh.N[1]-2] * 0.5);
    }

    std::cout << "Global Energy Balance: " << sumQ << " " << sumBC << " " << sumQ - sumBC << "\n";

}


void Medic::getSystemResidual(Material Mat, Mesh Msh, Discretizer Dsc, double t){
    
    // Control
    double tempRes; int k{};
    fileR << t;

    // Interior Nodes
    for (int i = 1; i < Msh.N[0]-1; i++){
        for (int j = 1; j < Msh.N[1]-1; j++){
            
            // Control
            k = i * Msh.N[1] + j;
            
            // Calculate
            tempRes = Msh.matA[k].ap * Msh.vPhi[i][j] + Msh.matA[k].aw * Msh.vPhi[i-1][j] + Msh.matA[k].ae * Msh.vPhi[i+1][j] + Msh.matA[k].as * Msh.vPhi[i][j-1] + Msh.matA[k].an * Msh.vPhi[i][j+1] - Msh.bp[k];
            
            fileR << "," << tempRes;

        }
    } fileR << "\n";
}

Medic::~Medic(){

	// Close Files
	file.close();
	fileR.close();

}
