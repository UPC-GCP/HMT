// Imports
#include <iostream>
#include <vector>
#include <string>
#include <json/json.h>
#include <cmath>
#include <math.h>
#include <numeric>
#include <ctime>
#include <algorithm>

// Self-Imports
#include "Material.h"
#include "Mesh.h"
#include "Discretizer.h"
#include "Solver.h"
#include "GS.h"
#include "libArithmetic.h"

void GS::newSolve(std::vector<Matrix> matA, std::vector<std::vector<double>>& x, std::vector<double> matB, std::vector<std::vector<int>> ignoreBC){

    // PENDING UPDATE TO CHANGES FROM OPERPRODMATVEC
    // QUE LATA TENER QUE HACER ESTO

    // // Control
    // int n = matB.size(), m = x.size(), l = x[0].size(), iPos; double alpha, rsNew, beta, tempErr = 1;
    // std::vector<std::vector<double>> xOld;
    // std::vector<double> vRes(n);

    // // Gauss-Seidel Loop
    // for (int k = 0; k < maxIter; k++){

    //     // std::cout << "GS Loop: " << k << "\n";

    //     // Nodes
    //     for (int i = 0; i < m; i++){
    //         for (int j = 0; j < l; j++){

    //             // Control
    //             if (std::count(ignoreBC.begin(), ignoreBC.end(), std::vector<int>{i, j})){continue;}
    //             xOld = x;
    //             iPos = i * l + j;

    //             // std::cout << "Nodes Loop: " << i << " " << j << " " << iPos << "\n";

    //             // Update Values
    //             if (i > 0 && i < m-1 && j > 0 && j < l-1){
    //                 // Interior Nodes
    //                 x[i][j] = (- matA[iPos].aw * x[i-1][j] - matA[iPos].ae * x[i+1][j] - matA[iPos].as * x[i][j-1] - matA[iPos].an * x[i][j+1] + matB[iPos]) / matA[iPos].ap;
    //             } else if (i == 0){
    //                 // West Boundary - ae, ap
    //                 x[i][j] = (- matA[iPos].ae * x[i+1][j] + matB[iPos]) / matA[iPos].ap;
    //             } else if (i == n-1){
    //                 // East Boundary - aw, ap
    //                 x[i][j] = (- matA[iPos].aw * x[i-1][j] + matB[iPos]) / matA[iPos].ap;
    //             } else if (j == 0){
    //                 // South Boundary - an, ap
    //                 x[i][j] = (- matA[iPos].an * x[i][j+1] + matB[iPos]) / matA[iPos].ap;
    //             } else if (j == m-1){
    //                 // North Boundary - as, ap
    //                 x[i][j] = (- matA[iPos].as * x[i][j-1] + matB[iPos]) / matA[iPos].ap;
    //             }

    //         }
    //     }

    //     // Error
    //     vRes = operCombLinVec(newProdMatVec(matA, x), matB, 1, -1);
    //     tempErr = std::sqrt(operDotProd(vRes, vRes));
    //     if (tempErr < tolNum){break;}

    // }

}