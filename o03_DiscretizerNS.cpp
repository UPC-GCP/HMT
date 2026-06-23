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



void Discretizer::setMomentumCoefficients(Material Mat, Mesh& Msh){

    // Control
    int k{};
    double Fw{}, Fe{}, Fs{}, Fn{}, Dw{}, De{}, Ds{}, Dn{}, Pw{}, Pe{}, Ps{}, Pn{}, a0{};
    double mu = Mat.vMat[Msh.p.nMat[0][0]].mu, rho = Mat.vMat[Msh.p.nMat[0][0]].rho; // simplification because of single material

    // u Interior Nodes
    for (size_t i = 1; i < Msh.u.N[0]-1; i++){
        for (size_t j = 1; j < Msh.u.N[1]-1; j++){
            // Control
            k = i * Msh.u.N[1] + j;

            // Coefficients Convection-Diffusion 
            Dw = mu * Msh.u.Sw[i][j] / Msh.u.dX[0][i]; Fw = rho * Msh.u.Sw[i][j] * 0.5 * (Msh.u.Phi[i][j] + Msh.u.Phi[i-1][j]);
            Pw = Fw / Dw; Msh.u.matA[k].aw = Dw * funcScheme(Pw) + std::max(Fw, 0.0);
            De = mu * Msh.u.Se[i][j] / Msh.u.dX[0][i+1]; Fe = rho * Msh.u.Se[i][j] * 0.5 * (Msh.u.Phi[i][j] + Msh.u.Phi[i+1][j]);
            Pe = Fe / De; Msh.u.matA[k].ae = De * funcScheme(Pe) + std::max(-Fe, 0.0);
            Ds = mu * Msh.u.Ss[i][j] / Msh.u.dX[1][j]; Fs = rho * Msh.u.Ss[i][j] * 0.5 * (Msh.v.Phi[i][j] + Msh.v.Phi[i-1][j]);
            Ps = Fs / Ds; Msh.u.matA[k].as = Ds * funcScheme(Ps) + std::max(Fs, 0.0);
            Dn = mu * Msh.u.Sn[i][j] / Msh.u.dX[1][j+1]; Fn = rho * Msh.u.Sn[i][j] * 0.5 * (Msh.v.Phi[i][j+1] + Msh.v.Phi[i-1][j+1]);
            Pn = Fn / Dn; Msh.u.matA[k].an = Dn * funcScheme(Pn) + std::max(-Fn, 0.0);

            // Calculate
            a0 = rho * Msh.u.Vp[i][j] / dt;
            Msh.u.Phi[i][j] = Msh.u.oPhi[i][j] + (1 / a0) * (Msh.u.matA[k].aw * Msh.u.oPhi[i-1][j] + Msh.u.matA[k].ae * Msh.u.oPhi[i+1][j] + Msh.u.matA[k].as * Msh.u.oPhi[i][j-1] + Msh.u.matA[j].an * Msh.u.oPhi[i][j+1] - (Msh.u.matA[k].aw + Msh.u.matA[k].ae + Msh.u.matA[k].as + Msh.u.matA[k].an) * Msh.u.oPhi[i][j]);
        }
    }

    // v Interior Nodes
    for (size_t i = 1 ; i < Msh.v.N[0]-1; i++){
        for (size_t j = 1; j < Msh.v.N[1]-1; j++){
            // Control
            k = i * Msh.v.N[1] + j;

            // Coefficients Convection-Diffusion
            Dw = mu * Msh.v.Sw[i][j] / Msh.v.dX[0][i]; Fw = rho * Msh.v.Sw[i][j] * 0.5 * (Msh.u.Phi[i][j] + Msh.u.Phi[i][j-1]);
            Pw = Fw / Dw; Msh.u.matA[k].aw = Dw * funcScheme(Pw) + std::max(Fw, 0.0);
            De = mu * Msh.v.Se[i][j] / Msh.v.dX[0][i+1]; Fe = rho * Msh.v.Se[i][j] * 0.5 * (Msh.u.Phi[i+1][j] + Msh.u.Phi[i+1][j-1]);
            Pe = Fe / De; Msh.u.matA[k].ae = De * funcScheme(Pe) + std::max(-Fe, 0.0);
            Ds = mu * Msh.v.Ss[i][j] / Msh.v.dX[1][j]; Fs = rho * Msh.v.Ss[i][j] * 0.5 * (Msh.v.Phi[i][j] + Msh.v.Phi[i][j-1]);
            Ps = Fs / Ds; Msh.u.matA[k].as = Ds * funcScheme(Ps) + std::max(Fs, 0.0);
            Dn = mu * Msh.v.Sn[i][j] / Msh.v.dX[1][j+1]; Fn = rho * Msh.v.Sn[i][j] * 0.5 * (Msh.v.Phi[i][j] + Msh.v.Phi[i][j+1]);
            Pn = Fn / Dn; Msh.u.matA[k].an = Dn * funcScheme(Pn) + std::max(-Fn, 0.0);

            // Calculate
            a0 = rho * Msh.u.Vp[i][j] / dt;
            Msh.u.Phi[i][j] = Msh.u.oPhi[i][j] + (1 / a0) * (Msh.u.matA[k].aw * Msh.u.oPhi[i-1][j] + Msh.u.matA[k].ae * Msh.u.oPhi[i+1][j] + Msh.u.matA[k].as * Msh.u.oPhi[i][j-1] + Msh.u.matA[j].an * Msh.u.oPhi[i][j+1] - (Msh.u.matA[k].aw + Msh.u.matA[k].ae + Msh.u.matA[k].as + Msh.u.matA[k].an) * Msh.u.oPhi[i][j]);
        }
    }

}


void Discretizer::setPressureCoefficients(Material Mat, Mesh& Msh){

    // Control
    int k{};
    double Dw{}, De{}, Ds{}, Dn{}, a0{};
    double mu = Mat.vMat[Msh.p.nMat[0][0]].mu, rho = Mat.vMat[Msh.p.nMat[0][0]].rho; // simplification because of single material

    // p Interior Nodes
    for (size_t i = 1; i < Msh.p.N[0]-1; i++){
        for (size_t j = 1; j < Msh.p.N[1]-1; j++){
            // Control
            k = i * Msh.p.N[1] + j;

            // Coefficients A
            Dw = Msh.p.Sw[i][j] / Msh.p.dX[0][i]; De = Msh.p.Se[i][j] / Msh.p.dX[0][i+1];
            Ds = Msh.p.Ss[i][j] / Msh.p.dX[1][j]; Dn = Msh.p.Sn[i][j] / Msh.p.dX[1][j+1];
            Msh.p.matA[k].aw = - Dw; Msh.p.matA[k].ae = - De; Msh.p.matA[k].as = - Ds; Msh.p.matA[k].an = - Dn;
            Msh.p.matA[k].ap = Msh.p.matA[k].aw + Msh.p.matA[k].ae + Msh.p.matA[k].as + Msh.p.matA[k].an;
    
            // Coefficients B
            Msh.p.matB[k] = (rho / dt) * (Msh.u.Phi[i+1][j] * Msh.p.Se[i][j] - Msh.u.Phi[i][j] * Msh.p.Sw[i][j] + Msh.v.Phi[i][j+1] * Msh.p.Sn[i][j] - Msh.v.Phi[i][j] * Msh.p.Ss[i][j]);
        }
    }

}


void Discretizer::setMomentumBoundaries(Material Mat, Mesh& Msh){

    // Control
    int i{}, j{}, k{}; double gammaw{}, gammae{}, gammas{}, gamman{}, a0{};
    double Fw{}, Fe{}, Fs{}, Fn{}, Dw{}, De{}, Ds{}, Dn{}, Pw{}, Pe{}, Ps{}, Pn{};

    // Boundary Node Coefficients
    for (boundVelocity& bC : Msh.boundaryVelocity){



    }


}


void Discretizer::setPressureBoundaries(Material Mat, Mesh& Msh){

    // Control
    Msh.p.tempA.assign(Msh.p.totNodes, {}); Msh.p.tempB.assign(Msh.p.totNodes, 0);
    int i{}, j{}, k{}; double gammaw{}, gammae{}, gammas{}, gamman{}, a0{};
    double Fw{}, Fe{}, Fs{}, Fn{}, Dw{}, De{}, Ds{}, Dn{}, Pw{}, Pe{}, Ps{}, Pn{};
    std::vector<int> Pos0{}, Pos1{}; Pos0.resize(Msh.p.N.size()); Pos1.resize(Msh.p.N.size());

    // Boundary Node Coefficients
    for (boundMain& bC : Msh.boundaryConditions){

        // Positions
        for (size_t i = 0; i < Msh.p.N.size(); i++){
            Pos0[i] = std::lower_bound(Msh.p.Nodes[i].begin(), Msh.p.Nodes[i].end(), bC.x0[i] - epsFind) - Msh.p.Nodes[i].begin();
            Pos1[i] = std::lower_bound(Msh.p.Nodes[i].begin(), Msh.p.Nodes[i].end(), bC.x1[i] - epsFind) - Msh.p.Nodes[i].begin();
        }

        // Boundaries
        if (bC.type == 0){

            // Dirichlet
            if (Pos0[0] == Pos1[0]){

                // xBoundary
                if (bC.side == 0){

                    // West Boundary

                } else if (bC.side == 1){

                    // East Boundary

                } else {std::cerr << "Boundary side not specified correctly.\n";}


            } else if (Pos0[1] == Pos1[1]){

                // yBoundary
                if (bC.side == 0){

                    // South Boundary

                } else if (bC.side == 1){

                    // North Boundary

                } else {std::cerr << "Boundary side not specified correctly.\n";}

            } else {std::cerr << "Boundary range not specified correctly.\n";}

        } else if (bC.type == 1){

            // Neumann
            if (Pos0[0] == Pos1[0]){

                // xBoundary
                if (bC.side == 0){

                    // West Boundary

                } else if (bC.side == 1){

                    // East Boundary

                } else {std::cerr << "Boundary side not specified correctly.\n";}

            } else if (Pos0[1] == Pos1[1]){

                // yBoundary
                if (bC.side == 0){

                    // South Boundary

                } else if (bC.side == 1){

                    // North Boundary

                } else {std::cerr << "Boundary side not specified correctly.\n";}

            } else {std::cerr << "Boundary range not specified correctly.\n";}

        } else {std::cerr << "Boundary type not specified correctly.\n";}

    }

    // Corrections (Change 0,0 to Dirichlet for stability)


}


void Discretizer::correctVelocity(Material Mat, Mesh& Msh){

    // Control
    double mu = Mat.vMat[Msh.p.nMat[0][0]].mu, rho = Mat.vMat[Msh.p.nMat[0][0]].rho; // simplification because of single material

    // u Interior Nodes
    for (size_t i = 1; i < Msh.u.N[0]-1; i++){
        for (size_t j = 1; j < Msh.u.N[1]-1; j++){
            Msh.u.Phi[i][j] -= (dt / rho) * (Msh.p.Phi[i][j] - Msh.p.Phi[i-1][j]) / Msh.u.dX[0][i];
        }
    }

    // v Interior Nodes
    for (size_t i = 1; i < Msh.v.N[0]-1; i++){
        for (size_t j = 1; j < Msh.v.N[1]-1; j++){
            Msh.v.Phi[i][j] -= (dt / rho) * (Msh.p.Phi[i][j] - Msh.p.Phi[i][j-1]) / Msh.v.dX[1][j];    
        }
    }

}


