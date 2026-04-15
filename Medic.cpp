// Imports
#include <iostream>
#include <vector>
#include <fstream>
#include <filesystem>

#include "Material.h"
#include "Mesh.h"
#include "Discretizer.h"
#include "Probe.h"
#include "Medic.h"

std::ofstream createDiagnostic(std::string fName){

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

    // Energy Balance
    file = createDiagnostic(Prb.pathBase + "Medic_Balance.csv");
    for (size_t i = 1; i < Msh.N[0]-1; i++){
        for (size_t j = 1; j < Msh.N[1]-1; j++){
            file << "," << Msh.Nodes[0][i] << " " << Msh.Nodes[1][j];
        }
    } file << "\n";

    // Residue
    fileR = createDiagnostic(Prb.pathBase + "Medic_Residue.csv");
    for (size_t i = 1; i < Msh.N[0]-1; i++){
        for (size_t j = 1; j < Msh.N[1]-1; j++){
            fileR << "," << Msh.Nodes[0][i] << " " << Msh.Nodes[1][j];
        }
    } fileR << "\n";

}

void Medic::getDiagnostic(Material Mat, Mesh Msh, Discretizer Dsc, std::vector<std::vector<double>> oldTemp, double t){

    // Control
    double tempErr; int k{};
    double lamb{}, lambw{}, lambe{}, lambs{}, lambn{};
    file << t;

    // Nodes Loop
    for (int i = 1; i < Msh.N[0]-1; i++){
        for (int j = 1; j < Msh.N[1]-1; j++){
            
            // Control
            k = i * Msh.N[1] + j;

            // Harmonic Mean
            lamb = Mat.vMat[Msh.nMat[i][j]].lambda;
            lambw = Dsc.calcHarmonicMean(Msh.nd[0][i-1], {lamb, Mat.vMat[Msh.nMat[i-1][j]].lambda}, {Msh.ndelta[0][i], Msh.ndelta[0][i-1]});
            lambe = Dsc.calcHarmonicMean(Msh.nd[0][i], {lamb, Mat.vMat[Msh.nMat[i+1][j]].lambda}, {Msh.ndelta[0][i], Msh.ndelta[0][i+1]});
            lambs = Dsc.calcHarmonicMean(Msh.nd[1][j-1], {lamb, Mat.vMat[Msh.nMat[i][j-1]].lambda}, {Msh.ndelta[1][j], Msh.ndelta[1][j-1]});
            lambn = Dsc.calcHarmonicMean(Msh.nd[1][j], {lamb, Mat.vMat[Msh.nMat[i][j+1]].lambda}, {Msh.ndelta[1][j], Msh.ndelta[1][j+1]});

            // Calculate Error
            tempErr = Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp * Msh.nVp[i][j] * (Msh.nT[i][j] - oldTemp[i][j]) / Dsc.dt - lambw * Msh.nSw[i][j] * (Msh.nT[i][j] - Msh.nT[i-1][j]) / Msh.nd[0][i-1] - lambe * Msh.nSe[i][j] * (Msh.nT[i][j] - Msh.nT[i+1][j]) / Msh.nd[0][i] - lambs * Msh.nSs[i][j] * (Msh.nT[i][j] - Msh.nT[i][j-1]) / Msh.nd[1][j-1] - lambn * Msh.nSn[i][j] * (Msh.nT[i][j] - Msh.nT[i][j+1]) / Msh.nd[1][j] + Msh.nQv[i][j] * Msh.nVp[i][j];
            
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
            sumQ += Msh.nQv[i][j] * Msh.nVp[i][j];
        }
    }

    // Outward flux
    double sumBC = 0;
    
    // xBoundaries (west, east)
    for (size_t i = 1; i < Msh.N[1]-1; i++){
        sumBC +=  Mat.vMat[Msh.nMat[1][i]].lambda * Msh.nSw[1][i] * (Msh.nT[1][i] - Msh.nT[0][i]) / (Msh.ndelta[0][1] * 0.5);
        sumBC += Mat.vMat[Msh.nMat[Msh.N[0]-2][i]].lambda * Msh.nSe[Msh.N[0]-2][i] * (Msh.nT[Msh.N[0]-2][i] - Msh.nT[Msh.N[0]-1][i]) / (Msh.ndelta[0][Msh.N[0]-2] * 0.5);
    }

    // yBoundaries (south, north)
    for (size_t i = 1; i < Msh.N[0]-1; i++){
        sumBC += Mat.vMat[Msh.nMat[i][1]].lambda * Msh.nSs[i][1] * (Msh.nT[i][1] - Msh.nT[i][0]) / (Msh.ndelta[1][1] * 0.5);
        sumBC += Mat.vMat[Msh.nMat[i][Msh.N[1]-2]].lambda * Msh.nSn[i][Msh.N[1]-2] * (Msh.nT[i][Msh.N[1]-2] - Msh.nT[i][Msh.N[1]-2]) / (Msh.ndelta[1][Msh.N[1]-2] * 0.5);
    }

    std::cout << "Global Energy Balance: " << sumQ << " " << sumBC << " " << sumQ - sumBC << "\n";

}


void Medic::getSystemResidual(Material Mat, Mesh Msh, Discretizer Dsc){
    
    // Control
    double tempRes; int k{};

    // Interior Nodes
    for (int i = 1; i < Msh.N[0]-1; i++){
        for (int j = 1; j < Msh.N[1]-1; j++){
            
            // Control
            k = i * Msh.N[1] + j;
            
            // Calculate
            tempRes = Msh.matA[k].ap * Msh.nT[i][j] + Msh.matA[k].aw * Msh.nT[i-1][j] + Msh.matA[k].ae * Msh.nT[i+1][j] + Msh.matA[k].as * Msh.nT[i][j-1] + Msh.matA[k].an * Msh.nT[i][j+1] - Msh.bp[k];
            
            fileR << "," << tempRes;

        }
    } fileR << "\n";
}