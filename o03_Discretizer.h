#ifndef DISCRETIZER_H_
#define DISCRETIZER_H_

#include <string>
#include <functional>

#include "o01_Material.h"
#include "o02_Mesh.h"
#include "o09_ExpressionParser.h"


class Discretizer
{
private:

public:
    // Variables
    std::string tempScheme{}, spatScheme{};
    double beta{}, endTime{}, dt{}, epsFind{};
    bool bIgnore=true;
    std::function<double(double)> funcScheme{};

    // Constructor
    Discretizer(std::string temporalScheme, std::string spatialScheme, double endTime, double dt, double epsFind=1e-5);
    
    // Functions
    double calcHarmonicMean(double dPF, std::vector<double> lambda, std::vector<double> deltaX);
    void setSchemeParameters(Material& Mat, Mesh& Msh);
    void setBoundaryConditions(Material& Mat, Mesh& Msh, ExpressionParser& Prs, double t = 0);
    void setCoefficients(Material& Mat, Mesh& Msh);
    void setRHS(Material& Mat, Mesh& Msh);

    void newSetRHS(Material& Mat, Mesh& Msh);
    void newSetCoefficients(Material& Mat, Mesh& Msh);
    void newSetBoundaries(Material& Mat, Mesh& Msh, ExpressionParser& Prs, double t = 0);

};

#endif
