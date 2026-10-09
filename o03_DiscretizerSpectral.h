#ifndef DISCRETIZERSPECTRAL_H_
#define DISCRETIZERSPECTRAL_H_

// Self-Imports
#include "o02_MeshSpectral.h"
#include "o03_Discretizer.h"

class DiscretizerSpectral : Discretizer {
private:

public:
    // Headers
    void initializeBurgers(MeshBurgers& Msh);
    void calculateEnergy(std::vector<double>& E, const std::vector<std::complex<double>>& uHat, double Re);
    static std::vector<std::complex<double>> calculateRHS(const std::vector<std::complex<double>>& u, double Re, double t);
    static void calculateConvectiveTerm(std::vector<std::complex<double>>& C, const std::vector<std::complex<double>>& u);
};

#endif
