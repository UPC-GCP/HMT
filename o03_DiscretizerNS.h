#ifndef DISCRETIZERNS_H_
#define DISCRETIZERNS_H_

#include <string>
#include <functional>

#include "o01_Material.h"
#include "o02_MeshNS.h"


class Discretizer
{
private:

public:
    // Variables
    bool bStep=true;
    std::string tempScheme{}, spatScheme{};
    double beta{}, endTime{}, dt{}, epsFind{};
    std::function<double(double)> funcScheme{};

    // Constructor
    Discretizer(std::string temporalScheme, std::string spatialScheme, double endTime, double dt, double epsFind=1e-5);
    
    // Functions

    // There needs to be a general function that runs the entire discretization program and changes depending on the amount of variables involved in this model
    // Will detect p, T, u, v, w and run the proper discretization functions for each case
    // Need to remake my functions to adapt to the templates (this is going to be hell)



    double calcHarmonicMean(double dPF, std::vector<double> lambda, std::vector<double> deltaX);
    void setSchemeParameters(Material& Mat, Mesh& Msh);
    void checkStability(Material Mat, Mesh& Msh);

    void setMomentumCoefficients(Material Mat, Mesh& Msh);
    void setMomentumBoundaries(Material Mat, Mesh& Msh);

    void setPressureCoefficients(Material Mat, Mesh& Msh);
    void setPressureBoundaries(Material Mat, Mesh& Msh);

    void setEnergyCoefficients(Material Mat, Mesh& Msh);
    void setEnergyBoundaries(Material Mat, Mesh& Msh);

    void correctVelocity(Material Mat, Mesh& Msh);

};

#endif
