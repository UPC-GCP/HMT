#ifndef DISCRETIZERSPECTRAL_H_
#define DISCRETIZERSPECTRAL_H_

// Self-Imports
#include "o02_MeshSpectral.h"
#include "o03_Discretizer.h"

class DiscretizerSpectral : Discretizer {
private:

public:
    // Headers
    void initializeBurgers(MeshBurgers Msh);
    static std::vector<std::complex<double>> calculateRHS(std::vector<std::complex<double>> const& u, double Re, double t);
    void calculateEnergy(std::vector<double>& E, const std::vector<std::complex<double>>& uHat, double Re);
};

#endif
