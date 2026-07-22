// Imports
#include <iostream>
#include <vector>
#include <string>
#include <json/json.h>
#include <cmath>
#include <math.h>
#include <ctime>

// Self-Imports
#include "o04_SolverSC.h"
#include "o09_libArithmeticSC.h"

Solver::Solver(std::string scheme, double maxIterations, double tolNum, double tolTime, std::string fName, std::string fSol){
    
    // Data
    maxIter = maxIterations; this->tolNum = tolNum; this->tolTemp = tolTime;

}

double Solver::calcErr(std::vector<std::vector<double>> matA, std::vector<std::vector<double>> matB, std::vector<std::vector<bool>> bObs){

    // Control
    size_t Nx = matA.size(), Ny = matA[0].size(), k{};
    std::vector<double> errVec; errVec.resize(Nx*Ny, 0);
    
    // Error
    for (size_t i = 0; i < Nx; i++){
        for (size_t j = 0; j < Ny; j++){
            if (bObs[i][j]){continue;} k = i * Ny + j; 
            errVec[k] = abs(matA[i][j] - matB[i][j]);
        }
    }

    // Norm
    double rsNew = operDotProd(errVec, errVec);
    
    return rsNew;

}

bool Solver::newSolve(std::vector<Matrix> matA, std::vector<std::vector<double>>& x, std::vector<double> matB, std::vector<std::vector<bool>> bObs){
    std::cout << "Virtual function \n";
    return false;
}
