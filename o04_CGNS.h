#ifndef CGNS_H_
#define CGNS_H_

#include <string>
#include <vector>

#include "o02_MeshNS.h"
#include "o04_SolverNS.h"

class CG: public Solver{
private:

public:
    // Constructor
    CG(std::string scheme, double maxIterations, double tolNum, double tolTime, std::string fName, std::string fSol) : Solver(scheme, maxIterations, tolNum, tolTime, fName, fSol){};

    // Functions
    bool newSolve(std::vector<Matrix> matA, std::vector<std::vector<double>>& x, std::vector<double> matB) override;
};

#endif
