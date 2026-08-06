#include <iostream>
#include <vector>
#include <cmath>

#include "o02_Mesh.h"
#include "o04_BCG.h"
#include "o09_libArithmetic.h"

bool BCG::newSolve(std::vector<Matrix> matA, std::vector<std::vector<double>>& x, std::vector<double> matB, std::vector<std::vector<int>> ignoreBC){

    // Dimensions
    int n = matB.size(), m = x.size(), l = x[0].size();

    // Control
    std::vector<double> xFlat(n);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < l; j++)
            xFlat[i*l + j] = x[i][j];

    // Residual
    std::vector<double> r = operCombLinVec(matB, operProdMatVec(matA, xFlat, m, l), 1.0, -1.0);

    // Shadow Residual
    std::vector<double> rHat = r;

    // Scalars
    double rhoOld = 1.0, alpha = 1.0, omega = 1.0, rhoNew{}, beta{};

    // Vectors
    std::vector<double> v(n, 0.0), p(n, 0.0), s(n), t(n);

    // Convergence Criterion
    double normB = std::sqrt(operDotProd(matB, matB));
    if (normB == 0.0) normB = 1.0;

    // BiCGSTAB loop
    for (int k = 0; k < maxIter; k++){

        rhoNew = operDotProd(rHat, r);
        if (std::abs(rhoNew) < 1e-300){
            std::cerr << "BiCGSTAB breakdown: rho = 0 at iteration " << k << "\n";
            lastIter = k; lastRes = std::sqrt(operDotProd(r, r)); return false;
        }

        // p = r + beta*(p - omega*v)
        beta = (rhoNew / rhoOld) * (alpha / omega);
        p = operCombLinVec(r, operCombLinVec(p, v, 1.0, -omega), 1.0, beta);

        // v = A·p
        v = operProdMatVec(matA, p, m, l);

        // alpha = rho / (r̂, v)
        double dotRhatV = operDotProd(rHat, v);
        if (std::abs(dotRhatV) < 1e-300){
            std::cerr << "BiCGSTAB breakdown: (rHat, v) = 0 at iteration " << k << "\n";
            lastIter = k; lastRes = std::sqrt(operDotProd(r, r)); return false;
        }
        alpha = rhoNew / dotRhatV;

        // s = r - alpha*v  (intermediate residual)
        s = operCombLinVec(r, v, 1.0, -alpha);

        // Control 
        double normS = std::sqrt(operDotProd(s, s));
        if (normS / normB < tolNum){
            xFlat = operCombLinVec(xFlat, p, 1.0, alpha);
            lastIter = k; lastRes = normS; break;
        }

        // t = A·s
        t = operProdMatVec(matA, s, m, l);

        // omega = (t,s)/(t,t)
        double dotTT = operDotProd(t, t);
        if (std::abs(dotTT) < 1e-300){
            std::cerr << "BiCGSTAB breakdown: (t,t) = 0 at iteration " << k << "\n";
            lastIter = k; lastRes = normS; return false;
        }
        omega = operDotProd(t, s) / dotTT;

        // x = x + alpha*p + omega*s
        xFlat = operCombLinVec(xFlat, operCombLinVec(p, s, alpha, omega), 1.0, 1.0);

        // r = s - omega*t
        r = operCombLinVec(s, t, 1.0, -omega);

        // Convergence / divergence check
        double normR = std::sqrt(operDotProd(r, r));
        lastIter = k; lastRes = normR;

        if (std::isnan(normR) || std::isinf(normR)){
            std::cerr << "BCGS diverges at iteration " << k << ", residual: " << normR << "\n";
            return false;
        }
        if (normR / normB < tolNum){ break; }

        rhoOld = rhoNew;
    }

    // Control
    for (int i = 0; i < m; i++)
        for (int j = 0; j < l; j++)
            x[i][j] = xFlat[i*l + j];

    return true;
}
