#ifndef CGSC_H_
#define CGSC_H_

#include <string>
#include <vector>

#include "o02_MeshSC.h"
#include "o04_SolverSC.h"

class CG: public Solver{
private:

public:
    // Constructor
    CG(std::string scheme, double maxIterations, double tolNum, double tolTime, std::string fName, std::string fSol) : Solver(scheme, maxIterations, tolNum, tolTime, fName, fSol){};

    // Functions
    bool newSolve(std::vector<Matrix> matA, std::vector<std::vector<double>>& x, std::vector<double> matB, std::vector<std::vector<bool>> bObs) override;
};

#endif
