#ifndef SOLVERNS_H_
#define SOLVERNS_H_

#include <string>
#include <vector>
/* #include <fstream> */

/* #include "o01_Material.h" */
#include "o02_MeshNS.h"
/* #include "o03_Discretizer.h" */

// PENDING CLEAN VARIABLES AND FUNCTIONS

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
    virtual bool newSolve(std::vector<Matrix> matA, std::vector<std::vector<double>>& x, std::vector<double> matB, std::vector<std::vector<int>> ignoreBC = {}) = 0;

    // Functions
    double calcErr(std::vector<std::vector<double>> matA, std::vector<std::vector<double>> matB);
};

#endif
