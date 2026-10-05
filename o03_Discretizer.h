#ifndef DISCRETIZER_H_
#define DISCRETIZER_H_

// Self-Imports
#include "o00_Globals.h"
/* #include "o01_Material.h" */
/* #include "o02_MeshDEV.h" */

class Discretizer
{
private:

public:

    // Decide what things stay here:
    // bStep -- What was it, is it general?
    // endTime, dt -- Should just be kept in main
    // beta, tempScheme, spatScheme -- DiscretizerFVM
    // epsFind -- Make ot a global variable that can be accessed from everywhere

    /* // Variables */
    /* bool bStep=true; */ 
    /* std::string tempScheme{}, spatScheme{}; */
    /* double beta{}, endTime{}, dt{}, epsFind{}; */
    /* std::function<double(double)> funcScheme{}; */

    /* // Constructor */
    /* Discretizer(std::string temporalScheme, std::string spatialScheme, double endTime, double dt, double epsFind=1e-5); */
    
    /* // Functions */
    /* double calcHarmonicMean(double dPF, std::vector<double> lambda, std::vector<double> deltaX); */
    /* void setSchemeParameters(Material& Mat, Mesh& Msh); */
    /* void checkStability(Material Mat, Mesh& Msh); */

    /* void setMomentumCoefficients(Material Mat, Mesh& Msh); */
    /* void setMomentumBoundaries(Material Mat, Mesh& Msh); */
    /* void setObstacles(Material Mat, Mesh& Msh); */

    /* void setPressureCoefficients(Material Mat, Mesh& Msh); */
    /* void setPressureBoundaries(Material Mat, Mesh& Msh); */

    /* void setEnergyCoefficients(Material Mat, Mesh& Msh); */
    /* void setEnergyBoundaries(Material Mat, Mesh& Msh); */

    /* void correctVelocity(Material Mat, Mesh& Msh); */

};

#endif
