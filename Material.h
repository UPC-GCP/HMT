#ifndef MATERIAL_H_
#define MATERIAL_H_

#include <vector>
#include <json/json.h>

struct MatPhys{
    double rho, lambda, cp, alpha;
    std::string rhoExpr, lambdaExpr, cpExpr; // Not currently in use. Kept for future variable properties implementation.
};

class Material
{
private:

public:
    // Variables
    double T0{};

    // Vectors
    std::vector<MatPhys> vMat{};

    // Constructor
    Material(Json::Value materials);
    
    // Functions
    void setInitialConditions(double initTemp);
};

#endif