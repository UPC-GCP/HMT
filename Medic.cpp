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

    // Create File
    pathBase = Prb.pathBase;
    file = createDiagnostic(Prb.pathBase + "Medic_Diagnostics.csv");
    
    // Headers
    int k{};
    for (size_t i = 1; i < Msh.N[0]-1; i++){
        for (size_t j = 1; j < Msh.N[1]-1; j++){
            file << "," << Msh.Nodes[0][i] << " " << Msh.Nodes[1][j];
        }
    } file << "\n";

}

void Medic::getDiagnostic(Material Mat, Mesh Msh, Discretizer Dsc, std::vector<std::vector<double>> oldTemp, double t){

    // Control
    double tempErr; int k{};
    file << t;

    // Nodes Loop
    for (int i = 1; i < Msh.N[0]-1; i++){
        for (int j = 1; j < Msh.N[1]-1; j++){
            
            // Control
            k = i * Msh.N[1] + j;

            // Calculate Error
            tempErr = Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp * Msh.nVp[i][j] * (Msh.nT[i][j] - oldTemp[i][j]) / Dsc.dt + Msh.matA[k].aw * (Msh.nT[i][j] - Msh.nT[i-1][j]) + Msh.matA[k].ae * (Msh.nT[i][j] - Msh.nT[i+1][j]) + Msh.matA[k].as * (Msh.nT[i][j] - Msh.nT[i][j-1]) + Msh.matA[k].an * (Msh.nT[i][j] - Msh.nT[i][j+1]) + Msh.bp[k];
            
            // Print to File
            file << "," << tempErr;

        }
    } file << "\n";

}