#ifndef MATERIAL_H_
#define MATERIAL_H_

#include <vector>
#include <json/json.h>

struct MatPhys{
    double rho=0, gamma=0, cp=1, mu=0;
};

class Material
{
private:

public:
    // Variables
    double Phi0{};

    // Vectors
    std::vector<MatPhys> vMat{};
    std::vector<double> VF0{};

    // Constructor
    Material(Json::Value materials);
    
    // Functions
    void setInitialConditions(double initPhi);
    void setInitialConditions(double initPhi, Json::Value initVF);
};

#endif
