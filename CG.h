#ifndef CG_H_
#define CG_H_

#include <string>
#include <vector>

#include "Material.h"
#include "Mesh.h"
#include "Discretizer.h"
#include "Solver.h"

class CG: public Solver{
private:

public:
    // Constructor
    CG(std::string scheme, double maxIterations, double tolNum, double tolTime, std::string fName, std::string fSol) : Solver(scheme, maxIterations, tolNum, tolTime, fName, fSol){};

    // Functions
    void newSolve(std::vector<Matrix> matA, std::vector<std::vector<double>>& x, std::vector<double> matB, std::vector<std::vector<int>> ignoreBC) override;

};

#endif