#ifndef SOLVERSPECTRAL_H_
#define SOLVERSPECTRAL_H_

// Self-Imports
#include "o04_Solver.h"

// Enumerator
enum class solSpectral { DNS, LES };

// Class
class SolverSpectral : Solver<1> {
private:

public:
    // Variables
    solSpectral algorithm{};

    // Constructor
    SolverSpectral(std::string solverSelection, double tolNumeric, double maxIterations);

    // Headers
    template <typename MeshClass> bool solveRK4(std::vector<std::complex<double>>& u, double Re, double t, double dt, MeshClass& Msh);
    template <typename MeshClass> bool solveLES(std::vector<std::complex<double>>& u, double Re, double t, double dt, MeshClass& Msh);

    // Functions
    template <typename MeshClass> bool solveModel(std::vector<std::complex<double>>& u, double Re, double t, double dt, MeshClass& Msh);
};

#endif
