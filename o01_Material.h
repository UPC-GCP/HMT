#ifndef MATERIAL_H_
#define MATERIAL_H_

#include <json/json.h>

struct MatPhys {
    double rho{}, gamma{}, cp{}, mu{}, beta{};
};

class Material {
private:

public:
    // Variables
    bool bPath=false; // Bool for Path/Value
    double P0{}, T0{}, Phi0{}, g{}; // Initial values for solver variables (p, T, phi), general constants (g)
    std::string sP0{}, sT0{}, sPhi0{}; // Path for solver variables

    // Vectors
    std::vector<MatPhys> vMat{}; // Material repository
    std::vector<double> VF0{}; // Initial values for Velocity Field
    std::vector<std::string> sVF0{}; // Path for velocity field

    // Constructor
    Material(Json::Value materials, double g=9.81);
    
    // Functions
    void setInitialConditions(double initPhi, Json::Value initVF); // PHISolver (Value) 
    void setInitialConditions(std::string pathPhi, Json::Value initVF); // PHISolver (Path)
    void setInitialConditions(double initT, double initP, Json::Value initVF); // NSSolver (Value) (Old: DHCSolver)
    void setInitialConditions(std::string pathT, std::string pathP, Json::Value pathV); // NSSolver (Path)
    
    // Deprecated (Pending cleanup)
    void setInitialConditions(double initPhi); // Old: PHISolver
};

#endif
