// Imports
/* #include "json/value.h" */
#include <csignal>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
/* #include <sys/ucontext.h> */
/* #include <iterator> */
/* #include <random> */
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


void Discretizer::checkStability(Material Mat, Mesh& Msh){

    // Control
    double dtMin{}, mu{}, rho{}, nu{};

    // uMesh
    for (size_t i = 0; i < Msh.u.N[0]; i++){
        for (size_t j = 0; j < Msh.u.N[1]; j++){

            // Properties 
            if (i == 0){ // Properties from East node
                mu = Mat.vMat[Msh.p.nMat[i][j]].mu; rho = Mat.vMat[Msh.p.nMat[i][j]].rho;
            } else { // Properties from West node
                mu = Mat.vMat[Msh.p.nMat[i-1][j]].mu; rho = Mat.vMat[Msh.p.nMat[i-1][j]].rho;
            } nu = mu / rho;

            // Convective
            dtMin = 0.35 * Msh.u.deltaX[0][i] / std::max(std::abs(Msh.u.Phi[i][j]), epsFind);
            if (dtMin < dt){dt = dtMin;}

            // Diffusive
            dtMin = (0.5 / nu) * std::pow(Msh.u.deltaX[0][i], 2) * std::pow(Msh.u.deltaX[1][j], 2) / (std::pow(Msh.u.deltaX[0][i], 2) + std::pow(Msh.u.deltaX[1][j], 2));
            if (dtMin < dt){dt = dtMin;}

        }
    }

    // vMesh
    for (size_t i = 0; i < Msh.v.N[0]; i++){
        for (size_t j = 0; j < Msh.v.N[1]; j++){

            // Properties
            if (j == 0){ // Properties from North node
                mu = Mat.vMat[Msh.p.nMat[i][j]].mu; rho = Mat.vMat[Msh.p.nMat[i][j]].rho;
            } else { // Properties from South node
                mu = Mat.vMat[Msh.p.nMat[i][j-1]].mu; rho = Mat.vMat[Msh.p.nMat[i][j-1]].rho;
            } nu = mu / rho;

            // Convective
            dtMin = 0.35 * Msh.v.deltaX[1][j] / std::max(std::abs(Msh.v.Phi[i][j]), epsFind);
            if (dtMin < dt){dt = dtMin;}

            // Diffusive
            dtMin = (0.5 / nu) * std::pow(Msh.v.deltaX[0][i], 2) * std::pow(Msh.v.deltaX[1][j], 2) / (std::pow(Msh.v.deltaX[0][i], 2) + std::pow(Msh.v.deltaX[1][j], 2));
            if (dtMin < dt){dt = dtMin;}

        }
    }
    
}

void Discretizer::setSchemeParameters(Material& Mat, Mesh& Msh){

    // Temporal Scheme
    if (tempScheme == "explicit") {
        beta = 0;
    } else if (tempScheme == "crank-nicolson") {
        beta = 0.5;
    } else if (tempScheme == "implicit") {
        beta = 1;
    } else {std::cerr << "Error: Invalid temporal discretization scheme " << tempScheme << "\n";}

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

    // Check Stability
    checkStability(Mat, Msh);

}



void Discretizer::setMomentumCoefficients(){


}


void Discretizer::setMomentumBoundaries(){


}


void Discretizer::setPressureCoefficients(){


}


void Discretizer::setPressureBoundaries(){


}


void Discretizer::correctVelocity(){


}










