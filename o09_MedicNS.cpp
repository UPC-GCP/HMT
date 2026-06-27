// Imports
#include <cerrno>
#include <cstddef>
#include <iostream>
#include <vector>
#include <fstream>
#include <filesystem>

#include "o01_Material.h"
#include "o02_MeshNS.h"
#include "o03_DiscretizerNS.h"
#include "o05_ProbeNS.h"
#include "o09_MedicNS.h"
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
    pathBase = Prb.pathBase; int k{};
    std::filesystem::path newPath(pathBase);

    // Energy Balance
    file = createDiagnostic(newPath / "Medic_Divergence.csv");
    for (size_t i = 0; i < Msh.p.N[0]; i++){
        for (size_t j = 0; j < Msh.p.N[1]; j++){
            file << "," << Msh.p.Nodes[0][i] << " " << Msh.p.Nodes[1][j];
        }
    } file << "\n";

    // Residue
    fileR = createDiagnostic(newPath / "Medic_Residue.csv"); 
    for (size_t i = 1; i < Msh.p.N[0]-1; i++){
        for (size_t j = 1; j < Msh.p.N[1]-1; j++){
            fileR << "," << Msh.p.Nodes[0][i] << " " << Msh.p.Nodes[1][j];
        }
    } fileR << "\n";

}


Medic::~Medic(){

	// Close Files
	file.close();
	fileR.close();

}


void Medic::getDiagnostic(Material Mat, Mesh Msh, Discretizer Dsc, ExpressionParser& Prs, double t){

    // Control
    double rho = Mat.vMat[Msh.p.nMat[0][0]].rho;
    file << t;

    // Write Values
    for (size_t i = 1; i < Msh.p.N[0]-1; i++){
        for (size_t j = 1; j < Msh.p.N[1]-1; j++){
            file << "," << (rho / Dsc.dt) * (Msh.u.Phi[i+1][j] * Msh.p.Se[i][j] - Msh.u.Phi[i][j] * Msh.p.Sw[i][j] + Msh.v.Phi[i][j+1] * Msh.p.Sn[i][j] - Msh.v.Phi[i][j] * Msh.p.Ss[i][j]);
        }
    } file << "\n";

}


void Medic::getSystemResidual(Material Mat, Mesh Msh, Discretizer Dsc, double t){

    // Control
    double tempRes; int k{};
    fileR << t;

    // Write Values
    for (size_t i = 1; i < Msh.p.N[0]-1; i++){
        for (size_t j = 1; j < Msh.p.N[1]-1; j++){
            k = i * Msh.p.N[1] + j;
            fileR << "," << Msh.p.matA[k].ap * Msh.p.Phi[i][j] + Msh.p.matA[k].aw * Msh.p.Phi[i-1][j] + Msh.p.matA[k].ae * Msh.p.Phi[i+1][j] + Msh.p.matA[k].as * Msh.p.Phi[i][j-1] + Msh.p.matA[k].an * Msh.p.Phi[i][j+1] - Msh.p.matB[k];
        }
    } fileR << "\n";
    

}


void Medic::getGlobalBalance(Material Mat, Mesh Msh, Discretizer Dsc){

    // Control
    double totDiv = 0, rho = Mat.vMat[Msh.p.nMat[0][0]].rho;

    // Calculate
    for (size_t i = 1; i < Msh.p.N[0]-1; i++){
        for (size_t j = 1; j < Msh.p.N[1]-1; j++){
            totDiv += Msh.u.Phi[i+1][j] * Msh.p.Se[i][j] - Msh.u.Phi[i][j] * Msh.p.Sw[i][j] + Msh.v.Phi[i][j+1] * Msh.p.Sn[i][j] - Msh.v.Phi[i][j] * Msh.p.Ss[i][j];
        }
    }

    // Control
    std::cout << "Global imbalance: " << (rho / Dsc.dt) * totDiv << "\n";

}

