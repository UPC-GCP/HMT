#ifndef MATERIAL_H_
#define MATERIAL_H_

#include <vector>
#include <json/json.h>

struct MatPhys{
    double rho{}, gamma{}, cp{}, mu{}, beta{};
};

class Material
{
private:

public:
    // Variables
    double P0{}, T0{}, Phi0{}, g{}; // Initial values for solver variables (p, T, phi), general constants (g)

    // Vectors
    std::vector<double> VF0{}; // Initial values for Velocity Field
    std::vector<MatPhys> vMat{}; // Material repository

    // Constructor
    Material(Json::Value materials, double g=9.81);
    
    // Functions
    void setInitialConditions(double initPhi); // MainSolver
    void setInitialConditions(double initPhi, Json::Value initVF); // NS Solver
    void setInitialConditions(double initT, double initP, Json::Value initVF); // DHCSolver
    void setInitialConditions(std::string pathT, std::string pathP, std::string pathV); // Read .csv
    template <typename nType, typename nJson> void setInitialConditions(nType initT, nType initP, nJson initVF);
    // PENDING CHANGES:
    // I want this to be able to read a .csv and use the last instant as the initial conditions of the simulation
    // config.json: Include path to .csv instead of value, code needs to identify the case
    // o01: if (Phi0 is number or path){read .csv; check dimensions for coherence with mesh (should also pass N); store initial value as vector}
};

#endif
