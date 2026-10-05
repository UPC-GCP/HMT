#ifndef SOLVERSPECTRAL_H_
#define SOLVERSPECTRAL_H_

// Self-Imports
#include "o04_Solver.h"

// Enumerator
enum class solSpectral { RK4, LES };

// Class
class SolverSpectral : Solver<1> {
private:

public:
    // Variables
    solSpectral algorithm{};

    // Constructor
    SolverSpectral(std::string solverSelection, double tolNumeric, double maxIterations);

    // Headers
    bool solveModel(std::vector<std::complex<double>>& u, double Re, double t, double dt);
    bool solveRK4(std::vector<std::complex<double>>& u, double Re, double t, double dt);
    bool solveLES(std::vector<std::complex<double>>& u, double Re, double t, double dt);
};

#endif
