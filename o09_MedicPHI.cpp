// Imports
#include <iostream>
#include <vector>
#include <fstream>
#include <filesystem>
#include <cmath>

#include "o01_Material.h"
#include "o02_Mesh.h"
#include "o03_Discretizer.h"
#include "o05_Probe.h"
#include "o09_Medic.h"
#include "o09_ExpressionParser.h"

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
    pathBase = Prb.pathBase; 
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
    fileR << ",Residue\n";

}


void Medic::getDiagnostic(Material Mat, Mesh Msh, Discretizer Dsc, ExpressionParser& Prs, double t){

    // Control
    int k{}; double tempVal{};
    file << t;

    // Nodes Loop
    for (int i = 0; i < Msh.N[0]; i++){
        for (int j = 0; j < Msh.N[1]; j++){
            
            // Control
            k = i * Msh.N[1] + j; tempVal = Msh.matA[k].ap * Msh.vPhi[i][j];

            // West Boundary
            if (i != 0){tempVal += Msh.matA[k].aw * Msh.vPhi[i-1][j];}

            // East Boundary
            if (i != Msh.N[0]-1){tempVal += Msh.matA[k].ae * Msh.vPhi[i+1][j];}

            // South Boundary
            if (j != 0){tempVal += Msh.matA[k].as * Msh.vPhi[i][j-1];}

            // North Boundary
            if (j != Msh.N[1]-1){tempVal += Msh.matA[k].an * Msh.vPhi[i][j+1];}
            
            // Print to File
            file << "," << (tempVal - Msh.bp[k]);

        }
    } file << "\n";

}


void Medic::getSystemResidual(Material Mat, Mesh Msh, Discretizer Dsc, double t){
    
    // Control
    double tempRes = 0, tempVal{}; int k{};
    fileR << t;

    // Nodes Loop
    for (int i = 0; i < Msh.N[0]; i++){
        for (int j = 0; j < Msh.N[1]; j++){
            
            // Control
            k = i * Msh.N[1] + j; tempVal = 0;

            // Calculate
            tempVal = Msh.matA[k].ap * Msh.vPhi[i][j] + (i > 0 ? Msh.matA[k].aw * Msh.vPhi[i-1][j] : 0) + (i < Msh.N[0]-1 ? Msh.matA[k].ae * Msh.vPhi[i+1][j] : 0) + (j > 0 ? Msh.matA[k].as * Msh.vPhi[i][j-1] : 0) + (j < Msh.N[1]-1 ? Msh.matA[k].an * Msh.vPhi[i][j+1] : 0);
            tempRes += std::pow(tempVal - Msh.bp[k], 2);

        }
    }

    // Write
    fileR << "," << std::sqrt(tempRes) << "\n";

}


void Medic::getGlobalBalance(Material Mat, Mesh Msh, Discretizer Dsc){

    // Internal Heat Generation
    double sumQ = 0;
    for (size_t i = 1; i < Msh.N[0]-1; i++){
        for (size_t j = 1; j < Msh.N[1]-1; j++){
            sumQ += Msh.sPhi[i][j] * Msh.Vp[i][j];
        }
    }

}


Medic::~Medic(){

	// Close Files
	file.close();
	fileR.close();

}
