#ifndef CG_H_
#define CG_H_

#include <string>
#include <vector>

/* #include "o01_Material.h" */
#include "o02_Mesh.h"
/* #include "o03_Discretizer.h" */
#include "o04_Solver.h"

class CG: public Solver{
private:

public:
    // Constructor
    CG(std::string scheme, double maxIterations, double tolNum, double tolTime, std::string fName, std::string fSol) : Solver(scheme, maxIterations, tolNum, tolTime, fName, fSol){};

    // Functions
    bool newSolve(std::vector<Matrix> matA, std::vector<std::vector<double>>& x, std::vector<double> matB, std::vector<std::vector<int>> ignoreBC) override;
};

#endif
