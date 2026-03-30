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
#include "CG.h"
#include "libArithmetic.h"

void CG::newSolve(std::vector<Matrix> matA, std::vector<std::vector<double>>& x, std::vector<double> matB, std::vector<std::vector<int>> ignoreBC){

    // Control
    int n = matB.size(), m = x.size(), l = x[0].size(), iPos; double alpha, rsNew, beta;
    std::vector<double> r = matB, Ap(n);

    // Residual
    std::vector<double> Ax = newProdMatVec(matA, x);
    for (int i = 0; i < n; i++){r[i] -= Ax[i];}

    // Direction
    std::vector<double> p = r;
    double rsOld = operDotProd(r, r);

    // Conjugate Gradient Loop
    for (int k = 0; k < maxIter; k++){

        // Step Size
        Ap = operProdMatVec(matA, p);
        alpha = rsOld / operDotProd(p, Ap);

        // Update Values
        for (int i = 0; i < m; i++){
            for (int j = 0; j < l; j++){
                iPos = i * l + j;
                x[i][j] += alpha * p[iPos]; r[iPos] -= alpha * Ap[iPos];
            }
        }
        
        // Control
        rsNew = operDotProd(r, r);

        // Error
        if (std::sqrt(rsNew) < tolNum){lastIter = k; lastRes = rsNew; break;}

        // Direction
        beta = rsNew / rsOld;
        for (int i = 0; i < n; i++){p[i] = r[i] + beta * p[i];}

        // Control
        rsOld = rsNew;

    }

}