// Imports
#include <complex.h>
#include <iostream>

// Self-Imports
#include "o04_SolverSpectral.h"
/* #include "o02_MeshSpectral.h" */
#include "o03_DiscretizerSpectral.h"

SolverSpectral::SolverSpectral(std::string solverSelection, double tolNumeric, double maxIterations) {
    // Data
    tolNum = tolNumeric; maxIter = maxIterations;

    std::cout << "Solver type: " << solverSelection << "\n";

    // Solver Selection
    if (solverSelection == "DNS") { algorithm = solSpectral::DNS; }   
    else if (solverSelection == "LES") { algorithm = solSpectral::LES; }
}

template <typename MeshClass> bool SolverSpectral::solveRK4(std::vector<std::complex<double>>& u, double Re, double t, double dt, MeshClass& Msh) {

    // Control
    std::vector<std::complex<double>> uTemp{}, STemp{}, SFinal{}; uTemp.resize(u.size()); STemp.resize(u.size()); SFinal.resize(u.size());

    // S1
    STemp = DiscretizerSpectral::calculateRHS(u, Re, t, algorithm, Msh);
    for (size_t k = 0; k < u.size(); k++) { SFinal[k] += dt * STemp[k] / static_cast<double>(6); }

    // S2
    for (size_t k = 0; k < u.size(); k++) { uTemp[k] = u[k] + 0.5 * dt * STemp[k]; }
    STemp = DiscretizerSpectral::calculateRHS(uTemp, Re, t + 0.5 * dt, algorithm, Msh);
    for (size_t k = 0; k < u.size(); k++) { SFinal[k] += dt * 2 * STemp[k] / static_cast<double>(6); }

    // S3
    for (size_t k = 0; k < u.size(); k++) { uTemp[k] = u[k] + 0.5 * dt * STemp[k]; }
    STemp = DiscretizerSpectral::calculateRHS(uTemp, Re, t + 0.5 * dt, algorithm, Msh);
    for (size_t k = 0; k < u.size(); k++) { SFinal[k] += dt * 2 * STemp[k] / static_cast<double>(6); }

    // S4
    for (size_t k = 0; k < u.size(); k++) { uTemp[k] = u[k] + dt * STemp[k]; }
    STemp = DiscretizerSpectral::calculateRHS(uTemp, Re, t + dt, algorithm, Msh);
    for (size_t k = 0; k < u.size(); k++) { SFinal[k] += dt * STemp[k] / static_cast<double>(6); }

    // Update
    for (size_t k = 0; k < u.size(); k++) { u[k] = u[k] + SFinal[k]; }

    return true;
}

template <typename MeshClass> bool SolverSpectral::solveModel(std::vector<std::complex<double>>& u, double Re, double t, double dt, MeshClass& Msh) {
    switch (algorithm) {
        case solSpectral::DNS: return solveRK4(u, Re, t, dt, Msh);
        case solSpectral::LES: return solveRK4(u, Re, t, dt, Msh);
        default: return false;
    }
}

// Instances
template bool SolverSpectral::solveModel<MeshBurgers>(std::vector<std::complex<double>>&, double, double, double, MeshBurgers&);
template bool SolverSpectral::solveModel<MeshLES>(std::vector<std::complex<double>>&, double, double, double, MeshLES&);
template bool SolverSpectral::solveRK4<MeshBurgers>(std::vector<std::complex<double>>&, double, double, double, MeshBurgers&);
template bool SolverSpectral::solveRK4<MeshLES>(std::vector<std::complex<double>>&, double, double, double, MeshLES&);
