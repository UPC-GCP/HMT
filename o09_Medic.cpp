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


double addBoundaries(Material Mat, Mesh Msh, Discretizer Dsc, ExpressionParser& Prs, int i, int j, double t, std::vector<std::vector<double>> oPhi, std::vector<std::vector<double>> obPhi){


    // Change of plans
    // This now uses the coefficients in tempA. Just need to save and use bcPhi
    // Can only use tempA for Dirichlet. Neumann needs to calculate it but both store bcPhi.

    // Control
    double gamma{}, tempVal{}, a0{}, tempStore{}, Fp{}, Dp{}, Pp{}; int k{};
    double Fw{}, Fe{}, Fs{}, Fn{}, Dw{}, De{}, Ds{}, Dn{}, Pw{}, Pe{}, Ps{}, Pn{};
    std::vector<int> Pos0{}, Pos1{}; Pos0.resize(Msh.N.size()); Pos1.resize(Msh.N.size());

    // Loop
    for (Boundary bC : Msh.boundaryConditions){

        // Positions
        for (size_t i = 0; i < Msh.N.size(); i++){
            Pos0[i] = std::lower_bound(Msh.Nodes[i].begin(), Msh.Nodes[i].end(), bC.x0[i] - Dsc.epsFind) - Msh.Nodes[i].begin();
            Pos1[i] = std::lower_bound(Msh.Nodes[i].begin(), Msh.Nodes[i].end(), bC.x1[i] - Dsc.epsFind) - Msh.Nodes[i].begin();
        }

        // Control 
        if (!(i >= Pos0[0] && i <=  Pos1[0] && j >= Pos0[1] && j <= Pos1[1])){continue;}
        k = i * Msh.N[1] + j;

        // Boundaries
        if (bC.type == 0){

            // Update Value
            if (bC.bUpdate && bC.iEq == 0){tempStore = Prs.evaluateTime(bC.iExpr, t - Dsc.dt);} else {tempStore = bC.value;}

            // Dirichlet
            if (Pos0[0] == Pos1[0]){

                // xBoundary
                if (bC.side == 0){
                    // West Boundary
                    tempVal = Dsc.beta * Msh.tempA[k].aw * (bC.value - Msh.vPhi[i][j]) + (1 - Dsc.beta) * Msh.tempA[k].aw * (obPhi[i][j] - oPhi[i][j]);
                } else if (bC.side == 1){
                    // East Boundary
                    tempVal = Dsc.beta * Msh.tempA[k].ae * (bC.value - Msh.vPhi[i][j]) + (1 - Dsc.beta) * Msh.tempA[k].ae * (obPhi[i][j] - oPhi[i][j]);
                }

            } else if (Pos0[1] == Pos1[1]){

                // yBoundary
                if (bC.side == 0){
                    // South Boundary
                    tempVal = Dsc.beta * Msh.tempA[k].as * (bC.value - Msh.vPhi[i][j]) + (1 - Dsc.beta) * Msh.tempA[k].as * (obPhi[i][j] - oPhi[i][j]);
                } else if (bC.side == 1){
                    // North Boundary
                    tempVal = Dsc.beta * Msh.tempA[k].an * (bC.value - Msh.vPhi[i][j]) + (1 - Dsc.beta) * Msh.tempA[k].an * (obPhi[i][j] - oPhi[i][j]);
                }

            }

        } else if (bC.type == 1){

            // Neumann
            if (Pos0[0] == Pos1[0]){

                // xBoundary
                if (bC.side == 0){
                    // West Boundary
                    tempVal = bC.value * Msh.Sw[i][j];
                } else if (bC.side == 1){
                    // East Boundary
                    tempVal = bC.value * Msh.Se[i][j];
                }

            } else if (Pos0[1] == Pos1[1]){

                // yBoundary
                if (bC.side == 0){
                    // South Boundary
                    tempVal = bC.value * Msh.Ss[i][j];
                } else if (bC.side == 1){
                    // North Boundary
                    tempVal = bC.value * Msh.Sn[i][j];
                }

            }

        } else if (bC.type == 2){

            // Robin

        }

    }

    return tempVal;

}


void Medic::getDiagnostic(Material Mat, Mesh Msh, Discretizer Dsc, ExpressionParser& Prs, std::vector<std::vector<double>> oldPhi, std::vector<std::vector<double>> bcPhi, double t){

    // Control
    double tempErr{}, impErr{}, expErr{}, extraErr{}; int k{};
    double lamb{}, lambw{}, lambe{}, lambs{}, lambn{};
    file << t;

    // Nodes Loop
    for (int i = 0; i < Msh.N[0]; i++){
        for (int j = 0; j < Msh.N[1]; j++){
            
            // Control
            k = i * Msh.N[1] + j; impErr = 0; expErr = 0; extraErr = 0;

            // West Boundary
            if (i != 0){
                impErr += Msh.tempA[k].aw * (Msh.vPhi[i-1][j] - Msh.vPhi[i][j]); 
                expErr += Msh.tempA[k].aw * (oldPhi[i-1][j] - oldPhi[i][j]);
            } else {extraErr += addBoundaries(Mat, Msh, Dsc, Prs, i, j, t, oldPhi, bcPhi);}

            // East Boundary
            if (i != Msh.N[0]-1){
                impErr += Msh.tempA[k].ae * (Msh.vPhi[i+1][j] - Msh.vPhi[i][j]);
                expErr += Msh.tempA[k].ae * (oldPhi[i+1][j] - oldPhi[i][j]);
            } else {extraErr += addBoundaries(Mat, Msh, Dsc, Prs, i, j, t, oldPhi, bcPhi);}

            // South Boundary
            if (j != 0){
                impErr += Msh.tempA[k].as * (Msh.vPhi[i][j-1] - Msh.vPhi[i][j]);
                expErr += Msh.tempA[k].as * (oldPhi[i][j-1] - oldPhi[i][j]);
            } else {extraErr += addBoundaries(Mat, Msh, Dsc, Prs, i, j, t, oldPhi, bcPhi);}

            // North Boundary
            if (j != Msh.N[1]-1){
                impErr += Msh.tempA[k].an * (Msh.vPhi[i][j+1] - Msh.vPhi[i][j]);
                expErr += Msh.tempA[k].an * (oldPhi[i][j+1] - oldPhi[i][j]);
            } else {extraErr += addBoundaries(Mat, Msh, Dsc, Prs, i, j, t, oldPhi, bcPhi);}

            // Calculate
            tempErr = Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp * Msh.Vp[i][j] * (Msh.vPhi[i][j] - oldPhi[i][j]) / Dsc.dt - Dsc.beta * impErr - (1 - Dsc.beta) * expErr + Msh.sPhi[i][j] * Msh.Vp[i][j] + extraErr;

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
