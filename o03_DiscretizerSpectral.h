#ifndef DISCRETIZERSPECTRAL_H_
#define DISCRETIZERSPECTRAL_H_

// Self-Imports
#include "o01_Material.h"
#include "o02_MeshSpectral.h"
#include "o03_Discretizer.h"
#include "o04_SolverSpectral.h"

class DiscretizerSpectral : Discretizer {
private:

public:
    // Headers
    void evaluateTimeStability(const Material& Mat, const MeshBurgers& Msh, double& timeStep, double C1);
    void initializeBurgers(MeshBurgers& Msh, double Re);
    void initializeBurgers(MeshLES& Msh, double Re);
    void calculateEnergy(std::vector<double>& E, const std::vector<std::complex<double>>& uHat);
    static void calculateConvectiveTerm(std::vector<std::complex<double>>& C, const std::vector<std::complex<double>>& u);
    template <typename MeshClass> static std::vector<std::complex<double>> calculateRHS(const std::vector<std::complex<double>>& u, double Re, double t, solSpectral solVar, MeshClass& Msh); 
    template <typename MeshClass> static void calculateVEff(MeshClass& Msh);
};

#endif
