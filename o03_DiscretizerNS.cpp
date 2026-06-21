// Imports
/* #include "json/value.h" */
#include <csignal>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
/* #include <iterator> */
/* #include <sys/ucontext.h> */
#include <vector>
#include <json/json.h>
#include <cmath>
#include <algorithm>

// Self-Imports
/* #include "exprtk.hpp" */
#include "o01_Material.h"
#include "o02_MeshNS.h"
#include "o03_DiscretizerNS.h"
/* #include "o09_ExpressionParser.h" */


double Discretizer::calcHarmonicMean(double dPF, std::vector<double> lambda, std::vector<double> deltaX) {

    // Denominator
    double A = 0;
    for (int i = 0; i < lambda.size(); i++){
        A += (deltaX[i] / 2) / lambda[i];
    }
    
    return dPF / A;

}


Discretizer::Discretizer(std::string temporalScheme, std::string spatialScheme, double endTime, double dt, double epsFind) {
    
    // Time Parameters
    tempScheme = temporalScheme; spatScheme = spatialScheme;
    this->endTime = endTime; this->dt = dt;
    this->epsFind = epsFind;

}


void checkStability(){

    // Temporal Scheme
    /* if (tempScheme == "explicit"){ */

    /* } else if (tempScheme == "crank-nicolson"){ */

    /* } else if (tempScheme == "implicit"){ */

    /* } else {std::cerr << "";} */

}

void Discretizer::setSchemeParameters(Material& Mat, Mesh& Msh){

    // Temporal Scheme
    if (tempScheme == "explicit") {
        
        // Beta
        beta = 0;

        // Calculate Timestep
        std::vector<double> dtNew(Msh.totNodes, 0); double dtMin{}; int k{};
        for (size_t i = 0; i < Msh.N[0]; i++){
            for (size_t j = 0; j < Msh.N[1]; j++){
                k = i * Msh.N[1] + j;
		        dtNew[k] = 0.5 * pow(Msh.DeltaX[0][i], 2) * pow(Msh.DeltaX[1][j], 2) / ((Mat.vMat[Msh.nMat[i][j]].gamma / (Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp)) * (pow(Msh.DeltaX[0][i], 2) + pow(Msh.DeltaX[1][j], 2)));
            }
        }

        // Update Time-step
        dtMin = *std::min_element(dtNew.begin()+1, dtNew.end()-1);
        if (dtMin < dt) {dt = dtMin;}

    } else if (tempScheme == "crank-nicolson") {

        // Beta
        beta = 0.5;

    } else if (tempScheme == "implicit") {

        // Beta
        beta = 1;

    } else {
        std::cerr << "Error: Invalid temporal discretization scheme " << tempScheme << "\n";
    }

    // Spatial Scheme
    if (spatScheme == "CDS"){
        funcScheme = [](double Pe) {return 1.0 - 0.5 * std::abs(Pe);};
    } else if (spatScheme == "UDS"){
        funcScheme = [](double Pe){return 1;};
    } else if (spatScheme == "Hybrid"){
        funcScheme = [](double Pe){return std::max(0.0, 1 -  0.5 * std::abs(Pe));};
    } else if (spatScheme == "PowerLaw"){
        funcScheme = [](double Pe){return std::max(0.0, std::pow(1 - 0.1 * std::abs(Pe), 5));};
    } else if (spatScheme == "Exponential"){
        funcScheme = [](double Pe){return Pe / (exp(std::abs(Pe)) - 1);};
    } else {
        std::cerr << "Error: Invalid spatial discretization scheme '" << spatScheme << "'\n";
    }

}
