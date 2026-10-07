// Imports
#include <complex.h>
/* #include <iostream> // Just to debug constructor rn */

// Self-Imports
#include "o04_SolverSpectral.h"
#include "o03_DiscretizerSpectral.h"

SolverSpectral::SolverSpectral(std::string solverSelection, double tolNumeric, double maxIterations) {
    // Data
    tolNum = tolNumeric; maxIter = maxIterations;

    // Solver Selection
    if (solverSelection == "RK4") { algorithm = solSpectral::RK4; }   
    else if (solverSelection == "LES") { algorithm = solSpectral::LES; }
}

bool SolverSpectral::solveRK4(std::vector<std::complex<double>>& u, double Re, double t, double dt) {

    // Control
    std::vector<std::complex<double>> uTemp{}, STemp{}, SFinal{}; uTemp.resize(u.size()); STemp.resize(u.size()); SFinal.resize(u.size());

    // S1
    STemp = DiscretizerSpectral::calculateRHS(u, Re, t);
    for (size_t k = 0; k < u.size(); k++) { SFinal[k] += dt * STemp[k] / static_cast<double>(6); }

    // S2
    for (size_t k = 0; k < u.size(); k++) { uTemp[k] = u[k] + 0.5 * dt * STemp[k]; }
    STemp = DiscretizerSpectral::calculateRHS(uTemp, Re, t + 0.5 * dt);
    for (size_t k = 0; k < u.size(); k++) { SFinal[k] += dt * 2 * STemp[k] / static_cast<double>(6); }

    // S3
    for (size_t k = 0; k < u.size(); k++) { uTemp[k] = u[k] + 0.5 * dt * STemp[k]; }
    STemp = DiscretizerSpectral::calculateRHS(uTemp, Re, t + 0.5 * dt);
    for (size_t k = 0; k < u.size(); k++) { SFinal[k] += dt * 2 * STemp[k] / static_cast<double>(6); }

    // S4
    for (size_t k = 0; k < u.size(); k++) { uTemp[k] = u[k] + dt * STemp[k]; }
    STemp = DiscretizerSpectral::calculateRHS(uTemp, Re, t + dt);
    for (size_t k = 0; k < u.size(); k++) { SFinal[k] += dt * STemp[k] / static_cast<double>(6); }

    // Update
    for (size_t k = 0; k < u.size(); k++) { u[k] = u[k] + SFinal[k]; }

    return true;
}

/* bool SolverSpectral::solveLES() { */
/*     return false; */
/* } */

bool SolverSpectral::solveModel(std::vector<std::complex<double>>& u, double Re, double t, double dt) {
    switch (algorithm) {
        case solSpectral::RK4: return solveRK4(u, Re, t, dt);
        case solSpectral::LES: return false;
        default: return false;
    }
}
