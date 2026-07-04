#ifndef MATERIAL_H_
#define MATERIAL_H_

#include <vector>
#include <json/json.h>

struct MatPhys{
    double rho=0, gamma=0, cp=1, mu=0, beta=0;
};

class Material
{
private:

public:
    // Variables
    double P0{}, T0{}, Phi0{}, g{};

    // Vectors
    std::vector<MatPhys> vMat{};
    std::vector<double> VF0{};

    // Constructor
    Material(Json::Value materials, double g=9.81);
    
    // Functions
    void setInitialConditions(double initPhi);
    void setInitialConditions(double initPhi, Json::Value initVF);
    void setInitialConditions(double initT, double initP, Json::Value initVF);
};

#endif
