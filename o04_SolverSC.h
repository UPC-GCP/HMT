#ifndef SOLVERSC_H_
#define SOLVERSC_H_

#include <string>
#include <vector>

#include "o02_MeshSC.h"

class Solver
{
private:
    
public:
    // Variables
    std::string fileName{};
    double maxIter{}, tolNum{}, tolTemp{}, lastIter{}, lastRes{};

    // Constructor
    Solver(std::string scheme, double maxIterations, double tolNum, double tolTime, std::string fName, std::string fSol);

    // Inheritance Functions
    virtual bool newSolve(std::vector<Matrix> matA, std::vector<std::vector<double>>& x, std::vector<double> matB, std::vector<std::vector<bool>> bObs) = 0;

    // Functions
    double calcErr(std::vector<std::vector<double>> matA, std::vector<std::vector<double>> matB, std::vector<std::vector<bool>> bObs);
};

#endif
