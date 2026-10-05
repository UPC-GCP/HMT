#ifndef SOLVERFVM_H_
#define SOLVERFVM_H_

// Imports
#include <cstddef>

// Self-Imports
#include "o04_Solver.h"
#include "o02_MeshFVM.h"

// Enumerator
enum class solFVM { CG, BCG, GS, TDMA };

// Class
template <size_t Dim> class SolverFVM : Solver<Dim> {
private:

public:
    // Variables
    solFVM algorithm{};

    // Constructor
    SolverFVM(std::string solverSelection, double tolNumeric, double maxIterations);

    // Headers
    bool solveModel(const std::vector<Matrix<Dim>>& A, std::vector<double>& x, const std::vector<double> B); // Do I need bObs here?
    bool solveCG();
    bool solveBCG();
};

#endif
