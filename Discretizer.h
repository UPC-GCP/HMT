#ifndef DISCRETIZER_H_
#define DISCRETIZER_H_

#include <string>

#include "Material.h"
#include "Mesh.h"
#include "ExpressionParser.h"

// PENDING CLEAN VARIABLES AND FUNCTIONS

class Discretizer
{
private:

public:
    // Variables
    std::string scheme{};
    double beta{}, endTime{}, dt{}, epsFind{};
    bool bIgnore=true;

    // Constructor
    Discretizer(std::string scheme, double endTime, double dt, double epsFind=1e-5);
    
    // Functions
    double calcHarmonicMean(double dPF, std::vector<double> lambda, std::vector<double> deltaX);
    void setSchemeParameters(Material& Mat, Mesh& Msh);
    void newSetBoundaryConditions(Material& Mat, Mesh& Msh, ExpressionParser& Prs, double t = 0);
    void newSetCoefficients(Material& Mat, Mesh& Msh);
    void newSetRHS(Material& Mat, Mesh& Msh);  
};

#endif