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
    double calcHarmonicMean(double dPF, std::vector<double> lambda, std::vector<double> deltaX);
    void setSchemeParameters(Material& Mat, Mesh& Msh);
    void checkStability(Material Mat, Mesh& Msh);

    void setMomentumCoefficients(Material Mat, Mesh& Msh);
    void setMomentumBoundaries(Material Mat, Mesh& Msh);

    void setPressureCoefficients(Material Mat, Mesh& Msh);
    void setPressureBoundaries(Material Mat, Mesh& Msh);

    void correctVelocity(Material Mat, Mesh& Msh);

};

#endif
