// Imports
#include <iostream>
#include <vector>
#include <cmath>

// Self-Imports
#include "o02_Mesh.h"
#include "o04_GS.h"
#include "o09_libArithmetic.h"

bool GS::newSolve(std::vector<Matrix> matA, std::vector<std::vector<double>>& x, std::vector<double> matB, std::vector<std::vector<int>> ignoreBC){

    // Dimensions
    int m = x.size(), l = x[0].size(), n = matB.size();

    // Gauss-Seidel Iterations
    for (int k = 0; k < maxIter; k++){

        // Forward sweep (row-major)
        for (int i = 0; i < m; i++){
            for (int j = 0; j < l; j++){
                int idx = i * l + j;
                double rhs = matB[idx];
                if (i > 0)   rhs -= matA[idx].aw * x[i-1][j];  // updated west
                if (i < m-1) rhs -= matA[idx].ae * x[i+1][j];  // old east
                if (j > 0)   rhs -= matA[idx].as * x[i][j-1];  // updated south
                if (j < l-1) rhs -= matA[idx].an * x[i][j+1];  // old north
                if (std::abs(matA[idx].ap) > 0.0)
                    x[i][j] = rhs / matA[idx].ap;
            }
        }

        // Compute residual  r = b - A*x
        std::vector<double> xFlat(n);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < l; j++)
                xFlat[i*l + j] = x[i][j];

        std::vector<double> Ax = operProdMatVec(matA, xFlat, m, l);
        double rsNew = 0.0;
        for (int i = 0; i < n; i++){
            double ri = matB[i] - Ax[i];
            rsNew += ri * ri;
        }

        lastIter = k; lastRes = std::sqrt(rsNew);

        if (std::isnan(rsNew) || std::isinf(rsNew)){
            std::cerr << "GS diverges @ iteration " << k << ", residual: " << rsNew << "\n";
            return false;
        }

        if (std::sqrt(rsNew) < tolNum){ return true; }
    }

    return true;
}
