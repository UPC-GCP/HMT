// Imports
#include <csignal>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <json/json.h>
#include <cmath>
#include <algorithm>

// Self-Imports
#include "o03_DiscretizerSC.h"


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
            dtMin = (0.25 / nu) * std::pow(Msh.u.deltaX[0][i], 2) * std::pow(Msh.u.deltaX[1][j], 2) / (std::pow(Msh.u.deltaX[0][i], 2) + std::pow(Msh.u.deltaX[1][j], 2));
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
            dtMin = (0.25 / nu) * std::pow(Msh.v.deltaX[0][i], 2) * std::pow(Msh.v.deltaX[1][j], 2) / (std::pow(Msh.v.deltaX[0][i], 2) + std::pow(Msh.v.deltaX[1][j], 2));
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
    double Fw{}, Fe{}, Fs{}, Fn{}, Dw{}, De{}, Ds{}, Dn{}, Pw{}, Pe{}, Ps{}, Pn{}, a0{}, Rn{};
    double mu = Mat.vMat[Msh.p.nMat[0][0]].mu, rho = Mat.vMat[Msh.p.nMat[0][0]].rho; // simplification because of single material
    double Cn = bStep ? 1 : 1.5, Co = bStep ? 0 : -0.5;

    // u Interior Nodes
    for (size_t i = 1; i < Msh.u.N[0]-1; i++){
        for (size_t j = 1; j < Msh.u.N[1]-1; j++){
            // Control
            k = i * Msh.u.N[1] + j;

            // Obstacles
            if (Msh.p.bObs[i-1][j] || Msh.p.bObs[i][j]){
                Msh.u.matA[k] = {}; Msh.u.matA[k].ap = 1;
                Msh.u.matB[k] = 0; Msh.u.Phi[i][j] = 0; Msh.u.oR[k] = 0;
                continue;
            }

            // Coefficients Convection-Diffusion 
            Dw = mu * Msh.u.Sw[i][j] / Msh.u.dX[0][i]; Fw = rho * Msh.u.Sw[i][j] * 0.5 * (Msh.u.oPhi[i][j] + Msh.u.oPhi[i-1][j]);
            Pw = Fw / Dw; Msh.u.matA[k].aw = Dw * funcScheme(Pw) + std::max(Fw, 0.0);
            De = mu * Msh.u.Se[i][j] / Msh.u.dX[0][i+1]; Fe = rho * Msh.u.Se[i][j] * 0.5 * (Msh.u.oPhi[i][j] + Msh.u.oPhi[i+1][j]);
            Pe = Fe / De; Msh.u.matA[k].ae = De * funcScheme(Pe) + std::max(-Fe, 0.0);
            Ds = mu * Msh.u.Ss[i][j] / Msh.u.dX[1][j]; Fs = rho * Msh.u.Ss[i][j] * 0.5 * (Msh.v.oPhi[i][j] + Msh.v.oPhi[i-1][j]);
            Ps = Fs / Ds; Msh.u.matA[k].as = Ds * funcScheme(Ps) + std::max(Fs, 0.0);
            Dn = mu * Msh.u.Sn[i][j] / Msh.u.dX[1][j+1]; Fn = rho * Msh.u.Sn[i][j] * 0.5 * (Msh.v.oPhi[i][j+1] + Msh.v.oPhi[i-1][j+1]);
            Pn = Fn / Dn; Msh.u.matA[k].an = Dn * funcScheme(Pn) + std::max(-Fn, 0.0);

            // Calculate
            a0 = rho * Msh.u.Vp[i][j] / dt; Rn = Msh.u.matA[k].aw * Msh.u.oPhi[i-1][j] + Msh.u.matA[k].ae * Msh.u.oPhi[i+1][j] + Msh.u.matA[k].as * Msh.u.oPhi[i][j-1] + Msh.u.matA[k].an * Msh.u.oPhi[i][j+1] - (Msh.u.matA[k].aw + Msh.u.matA[k].ae + Msh.u.matA[k].as + Msh.u.matA[k].an) * Msh.u.oPhi[i][j];
            Msh.u.Phi[i][j] = Msh.u.oPhi[i][j] + (1 / a0) * (Cn * Rn + Co * Msh.u.oR[k]);
            
            // Control
            Msh.u.oR[k] = Rn;
        }
    }

    // v Interior Nodes
    double beta = Mat.vMat[Msh.p.nMat[0][0]].beta, Tv{}, buoyancy{};
    for (size_t i = 1 ; i < Msh.v.N[0]-1; i++){
        for (size_t j = 1; j < Msh.v.N[1]-1; j++){
            // Control
            k = i * Msh.v.N[1] + j;

            // Obstacles
            if (Msh.p.bObs[i][j-1] || Msh.p.bObs[i][j]){
                Msh.v.matA[k] = {}; Msh.v.matA[k].ap = 1;
                Msh.v.matB[k] = 0; Msh.v.Phi[i][j] = 0; Msh.v.oR[k] = 0;
                continue;
            }

            // Coefficients Convection-Diffusion
            Dw = mu * Msh.v.Sw[i][j] / Msh.v.dX[0][i]; Fw = rho * Msh.v.Sw[i][j] * 0.5 * (Msh.u.oPhi[i][j] + Msh.u.oPhi[i][j-1]);
            Pw = Fw / Dw; Msh.v.matA[k].aw = Dw * funcScheme(Pw) + std::max(Fw, 0.0);
            De = mu * Msh.v.Se[i][j] / Msh.v.dX[0][i+1]; Fe = rho * Msh.v.Se[i][j] * 0.5 * (Msh.u.oPhi[i+1][j] + Msh.u.oPhi[i+1][j-1]);
            Pe = Fe / De; Msh.v.matA[k].ae = De * funcScheme(Pe) + std::max(-Fe, 0.0);
            Ds = mu * Msh.v.Ss[i][j] / Msh.v.dX[1][j]; Fs = rho * Msh.v.Ss[i][j] * 0.5 * (Msh.v.oPhi[i][j] + Msh.v.oPhi[i][j-1]);
            Ps = Fs / Ds; Msh.v.matA[k].as = Ds * funcScheme(Ps) + std::max(Fs, 0.0);
            Dn = mu * Msh.v.Sn[i][j] / Msh.v.dX[1][j+1]; Fn = rho * Msh.v.Sn[i][j] * 0.5 * (Msh.v.oPhi[i][j] + Msh.v.oPhi[i][j+1]);
            Pn = Fn / Dn; Msh.v.matA[k].an = Dn * funcScheme(Pn) + std::max(-Fn, 0.0);

            // Buoyancy
            if (!Msh.T.Phi.empty()){
                Tv = 0.5 * (Msh.T.Phi[i][j-1] + Msh.T.Phi[i][j]);
                buoyancy = beta * Mat.g * (Tv - Mat.T0) * dt;
            } else {buoyancy = 0.0;}

            // Calculate
            a0 = rho * Msh.v.Vp[i][j] / dt; Rn = Msh.v.matA[k].aw * Msh.v.oPhi[i-1][j] + Msh.v.matA[k].ae * Msh.v.oPhi[i+1][j] + Msh.v.matA[k].as * Msh.v.oPhi[i][j-1] + Msh.v.matA[k].an * Msh.v.oPhi[i][j+1] - (Msh.v.matA[k].aw + Msh.v.matA[k].ae + Msh.v.matA[k].as + Msh.v.matA[k].an) * Msh.v.oPhi[i][j];
            Msh.v.Phi[i][j] = Msh.v.oPhi[i][j] + (1 / a0) * (Cn * Rn + Co * Msh.v.oR[k]) + buoyancy;

            // Control
            Msh.v.oR[k] = Rn;
        }
    }

    // Control
    bStep = false;

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

            // Obstacles
            if (Msh.p.bObs[i][j]){
                Msh.p.matA[k] = {}; Msh.p.matA[k].ap = 1; Msh.p.matB[k] = 0;
                continue;
            }

            // Coefficients A
            Dw = Msh.p.Sw[i][j] / Msh.p.dX[0][i]; De = Msh.p.Se[i][j] / Msh.p.dX[0][i+1];
            Ds = Msh.p.Ss[i][j] / Msh.p.dX[1][j]; Dn = Msh.p.Sn[i][j] / Msh.p.dX[1][j+1];
            Msh.p.matA[k].aw = Msh.p.bObs[i-1][j] ? 0 : -Dw; Msh.p.matA[k].ae = Msh.p.bObs[i+1][j] ? 0 : -De;
            Msh.p.matA[k].as = Msh.p.bObs[i][j-1] ? 0 : -Ds; Msh.p.matA[k].an = Msh.p.bObs[i][j+1] ? 0 : -Dn;
            Msh.p.matA[k].ap = - Msh.p.matA[k].aw - Msh.p.matA[k].ae - Msh.p.matA[k].as - Msh.p.matA[k].an;
    
            // Coefficients B
            Msh.p.matB[k] = (rho / dt) * (Msh.u.Phi[i+1][j] * Msh.p.Se[i][j] - Msh.u.Phi[i][j] * Msh.p.Sw[i][j] + Msh.v.Phi[i][j+1] * Msh.p.Sn[i][j] - Msh.v.Phi[i][j] * Msh.p.Ss[i][j]);
        }
    }

}


void Discretizer::setMomentumBoundaries(Material Mat, Mesh& Msh){

    // Control
    int i{}, j{}, k{}; double gammaw{}, gammae{}, gammas{}, gamman{}, a0{}, mu = Mat.vMat[Msh.p.nMat[0][0]].mu, rho = Mat.vMat[Msh.p.nMat[0][0]].rho;
    double Fw{}, Fe{}, Fs{}, Fn{}, Dw{}, De{}, Ds{}, Dn{}, Pw{}, Pe{}, Ps{}, Pn{};

    // Boundary Node Coefficients
    for (boundVelocity& bC : Msh.boundaryVelocity){

        // Boundaries
        if (bC.type == 0){

            // Dirichlet
            if (bC.i0[0] == bC.i1[0]){

                // xBoundary
                if (bC.side == 0){
                    
                    // West Boundary
                    i = 0;

                    // uMesh
                    for (size_t j = 0; j < Msh.u.N[1]; j++){
                        // Control
                        k = i * Msh.u.N[1] + j; Msh.u.Phi[i][j] = bC.uVal;
                        Msh.u.matA[k].ap = 1; Msh.u.matB[k] = bC.uVal;
                    }

                    // vMesh
                    for (size_t j = 1; j < Msh.v.N[1]-1; j++){
                        // Control
                        k = i * Msh.v.N[1] + j;

                        // Coefficients Convection-Diffusion
                        Dw = mu * Msh.v.Sw[i][j] / Msh.v.dX[0][i]; Msh.v.matA[k].aw = Dw;

                        // Calculate
                        a0 = rho * Msh.v.Vp[i][j] / dt;
                        Msh.v.matB[k] = (1 / a0) * (Msh.v.matA[k].aw * bC.vVal);
                    }
                
                } else if (bC.side == 1){

                    // East Boundary
                    // uMesh
                    i = Msh.u.N[0]-1; 
                    for (size_t j = 0; j < Msh.u.N[1]; j++){
                        // Control
                        k = i * Msh.u.N[1] + j; Msh.u.Phi[i][j] = bC.uVal;
                        Msh.u.matA[k].ap = 1; Msh.u.matB[k] = bC.uVal;
                    }

                    // vMesh
                    i = Msh.v.N[0]-1;
                    for (size_t j = 1; j < Msh.v.N[1]-1; j++){
                        // Control
                        k = i * Msh.v.N[1] + j;

                        // Coefficients
                        De = mu * Msh.v.Se[i][j] / Msh.v.dX[0][i+1]; Msh.v.matA[k].ae = De;

                        // Calculate
                        a0 = rho * Msh.v.Vp[i][j] / dt;
                        Msh.v.matB[k] = (1 / a0) * (Msh.v.matA[k].ae * bC.vVal);
                    }

                }

            } else if (bC.i0[1] == bC.i1[1]){

                // yBoundary
                if (bC.side == 0){

                    // South Boundary
                    j = 0;

                    // uMesh
                    for (size_t i = 1; i < Msh.u.N[0]-1; i++){
                        // Control
                        k = i * Msh.u.N[1] + j;

                        // Coefficients
                        Ds = mu * Msh.u.Ss[i][j] / Msh.u.dX[1][j]; Msh.u.matA[k].as = Ds;

                        // Calculate
                        a0 = rho * Msh.u.Vp[i][j] / dt;
                        Msh.u.matB[k] = (1 / a0) * (Msh.u.matA[k].as * bC.uVal);
                    }

                    // vMesh
                    for (size_t i = 0; i < Msh.v.N[0]; i++){
                        k = i * Msh.v.N[1] + j; Msh.v.Phi[i][j] = bC.vVal;
                        Msh.v.matA[k].ap = 1; Msh.v.matB[k] = bC.vVal;
                    }

                } else if (bC.side == 1){

                    // North Boundary
                    // uMesh
                    j = Msh.u.N[1]-1;
                    for (size_t i = 1; i < Msh.u.N[0]-1; i++){
                        // Control
                        k = i * Msh.u.N[1] + j;

                        // Coefficients
                        Dn = mu * Msh.u.Sn[i][j] / Msh.u.dX[1][j+1]; Msh.u.matA[k].an = Dn;

                        // Calculate
                        a0 = rho * Msh.u.Vp[i][j] / dt;
                        Msh.u.matB[k] = (1 / a0) * (Msh.u.matA[k].an * bC.uVal);
                    }

                    // vMesh
                    j = Msh.v.N[1]-1;
                    for (size_t i = 0; i < Msh.v.N[0]; i++){
                        k = i * Msh.v.N[1] + j; Msh.v.Phi[i][j] = bC.vVal;
                        Msh.v.matA[k].ap = 1; Msh.v.matB[k] = bC.vVal;
                    }

                }

            }

        } else if (bC.type == 1){

            // Neumann
            if (bC.i0[0] == bC.i1[0]){

                // xBoundary
                if (bC.side == 0){

                    // West Boundary
                    // uMesh

                    // vMesh

                } else if (bC.side == 1){

                    // East Boundary
                    // uMesh
                    i = Msh.u.N[0]-1;
                    for (size_t j = 0; j < Msh.u.N[1]; j++){
                        Msh.u.Phi[i][j] = Msh.u.oPhi[i-1][j];
                    }
                    // vMesh
                    i = Msh.v.N[0]-1;
                    for (size_t j = 1; j < Msh.v.N[1]-1; j++){
                        k = i * Msh.v.N[1] + j;
                        Msh.v.matA[k].ae = 0; Msh.v.matB[k] = bC.vVal;
                    }

                } else {std::cerr << "Boundary side not specified correctly.\n";}

            } else if (bC.i0[1] == bC.i1[1]){

                // yBoundary
                if (bC.side == 0){

                    // South Boundary
                    // uMesh

                    // vMesh

                } else if (bC.side == 1){

                    // North Boundary
                    // uMesh

                    // vMesh

                } else {std::cerr << "Boundary side not specified correctly.\n";}

            } else {std::cerr << "Boundary range not specified correctly.\n";}

        } else {std::cerr << "Boundary type not recognized.\n";}

    }

    ///// Interior Node Coefficients
    
    /// West Boundary (Edges vMesh)
    i = 0; 
    for (size_t j = 1; j < Msh.v.N[1]-1; j++){
        // Control
        k = i * Msh.v.N[1] + j;

        // Coefficients Convection-Diffusion
        De = mu * Msh.v.Se[i][j] / Msh.v.dX[0][i+1]; Fe = rho * Msh.v.Se[i][j] * 0.5 * (Msh.u.oPhi[i+1][j] + Msh.u.oPhi[i+1][j-1]);
        Pe = Fe / De; Msh.v.matA[k].ae = De * funcScheme(Pe) + std::max(-Fe, 0.0);
        Ds = mu * Msh.v.Ss[i][j] / Msh.v.dX[1][j]; Fs = rho * Msh.v.Ss[i][j] * 0.5 * (Msh.v.oPhi[i][j] + Msh.v.oPhi[i][j-1]);
        Ps = Fs / Ds; Msh.v.matA[k].as = Ds * funcScheme(Ps) + std::max(Fs, 0.0);
        Dn = mu * Msh.v.Sn[i][j] / Msh.v.dX[1][j+1]; Fn = rho * Msh.v.Sn[i][j] * 0.5 * (Msh.v.oPhi[i][j] + Msh.v.oPhi[i][j+1]);
        Pn = Fn / Dn; Msh.v.matA[k].an = Dn * funcScheme(Pn) + std::max(-Fn, 0.0);

        // Calculate
        a0 = rho * Msh.v.Vp[i][j] / dt;
        Msh.v.Phi[i][j] = Msh.v.oPhi[i][j] + (1 / a0) * (Msh.v.matA[k].ae * Msh.v.oPhi[i+1][j] + Msh.v.matA[k].as * Msh.v.oPhi[i][j-1] + Msh.v.matA[k].an * Msh.v.oPhi[i][j+1] - (Msh.v.matA[k].aw + Msh.v.matA[k].ae + Msh.v.matA[k].as + Msh.v.matA[k].an) * Msh.v.oPhi[i][j]) + Msh.v.matB[k];
    }

    /// East Boundary (Edges vMesh)
    i = Msh.v.N[0]-1;
    for (size_t j = 1; j < Msh.v.N[1]-1; j++){
        // Control
        k = i * Msh.v.N[1] + j;

        // Coefficients Convection-Diffusion
        Dw = mu * Msh.v.Sw[i][j] / Msh.v.dX[0][i]; Fw = rho * Msh.v.Sw[i][j] * 0.5 * (Msh.u.oPhi[i][j] + Msh.u.oPhi[i][j-1]);
        Pw = Fw / Dw; Msh.v.matA[k].aw = Dw * funcScheme(Pw) + std::max(Fw, 0.0);
        Ds = mu * Msh.v.Ss[i][j] / Msh.v.dX[1][j]; Fs = rho * Msh.v.Ss[i][j] * 0.5 * (Msh.v.oPhi[i][j] + Msh.v.oPhi[i][j-1]);
        Ps = Fs / Ds; Msh.v.matA[k].as = Ds * funcScheme(Ps) + std::max(Fs, 0.0);
        Dn = mu * Msh.v.Sn[i][j] / Msh.v.dX[1][j+1]; Fn = rho * Msh.v.Sn[i][j] * 0.5 * (Msh.v.oPhi[i][j] + Msh.v.oPhi[i][j+1]);
        Pn = Fn / Dn; Msh.v.matA[k].an = Dn * funcScheme(Pn) + std::max(-Fn, 0.0);

        // Calculate
        a0 = rho * Msh.v.Vp[i][j] / dt;
        Msh.v.Phi[i][j] = Msh.v.oPhi[i][j] + (1 / a0) * (Msh.v.matA[k].aw * Msh.v.oPhi[i-1][j] + Msh.v.matA[k].as * Msh.v.oPhi[i][j-1] + Msh.v.matA[k].an * Msh.v.oPhi[i][j+1] - (Msh.v.matA[k].aw + Msh.v.matA[k].ae + Msh.v.matA[k].as + Msh.v.matA[k].an) * Msh.v.oPhi[i][j]) + Msh.v.matB[k];
    }

    /// South Boundary (Edges uMesh)
    j = 0;
    for (size_t i = 1; i < Msh.u.N[0]-1; i++){
        // Control
        k = i * Msh.u.N[1] + j;

        // Coefficients Convection-Diffusion 
        Dw = mu * Msh.u.Sw[i][j] / Msh.u.dX[0][i]; Fw = rho * Msh.u.Sw[i][j] * 0.5 * (Msh.u.oPhi[i][j] + Msh.u.oPhi[i-1][j]);
        Pw = Fw / Dw; Msh.u.matA[k].aw = Dw * funcScheme(Pw) + std::max(Fw, 0.0);
        De = mu * Msh.u.Se[i][j] / Msh.u.dX[0][i+1]; Fe = rho * Msh.u.Se[i][j] * 0.5 * (Msh.u.oPhi[i][j] + Msh.u.oPhi[i+1][j]);
        Pe = Fe / De; Msh.u.matA[k].ae = De * funcScheme(Pe) + std::max(-Fe, 0.0);
        Dn = mu * Msh.u.Sn[i][j] / Msh.u.dX[1][j+1]; Fn = rho * Msh.u.Sn[i][j] * 0.5 * (Msh.v.oPhi[i][j+1] + Msh.v.oPhi[i-1][j+1]);
        Pn = Fn / Dn; Msh.u.matA[k].an = Dn * funcScheme(Pn) + std::max(-Fn, 0.0);

        // Calculate
        a0 = rho * Msh.u.Vp[i][j] / dt;
        Msh.u.Phi[i][j] = Msh.u.oPhi[i][j] + (1 / a0) * (Msh.u.matA[k].aw * Msh.u.oPhi[i-1][j] + Msh.u.matA[k].ae * Msh.u.oPhi[i+1][j] + Msh.u.matA[k].an * Msh.u.oPhi[i][j+1] - (Msh.u.matA[k].aw + Msh.u.matA[k].ae + Msh.u.matA[k].as + Msh.u.matA[k].an) * Msh.u.oPhi[i][j]) + Msh.u.matB[k];
    }

    /// North Boundary (Edges uMesh)
    j = Msh.u.N[1]-1;
    for (size_t i = 1; i < Msh.u.N[0]-1; i++){
        // Control
        k = i * Msh.u.N[1] + j;

        // Coefficients Convection-Diffusion 
        Dw = mu * Msh.u.Sw[i][j] / Msh.u.dX[0][i]; Fw = rho * Msh.u.Sw[i][j] * 0.5 * (Msh.u.oPhi[i][j] + Msh.u.oPhi[i-1][j]);
        Pw = Fw / Dw; Msh.u.matA[k].aw = Dw * funcScheme(Pw) + std::max(Fw, 0.0);
        De = mu * Msh.u.Se[i][j] / Msh.u.dX[0][i+1]; Fe = rho * Msh.u.Se[i][j] * 0.5 * (Msh.u.oPhi[i][j] + Msh.u.oPhi[i+1][j]);
        Pe = Fe / De; Msh.u.matA[k].ae = De * funcScheme(Pe) + std::max(-Fe, 0.0);
        Ds = mu * Msh.u.Ss[i][j] / Msh.u.dX[1][j]; Fs = rho * Msh.u.Ss[i][j] * 0.5 * (Msh.v.oPhi[i][j] + Msh.v.oPhi[i-1][j]);
        Ps = Fs / Ds; Msh.u.matA[k].as = Ds * funcScheme(Ps) + std::max(Fs, 0.0);

        // Calculate
        a0 = rho * Msh.u.Vp[i][j] / dt;
        Msh.u.Phi[i][j] = Msh.u.oPhi[i][j] + (1 / a0) * (Msh.u.matA[k].aw * Msh.u.oPhi[i-1][j] + Msh.u.matA[k].ae * Msh.u.oPhi[i+1][j] + Msh.u.matA[k].as * Msh.u.oPhi[i][j-1] - (Msh.u.matA[k].aw + Msh.u.matA[k].ae + Msh.u.matA[k].as + Msh.u.matA[k].an) * Msh.u.oPhi[i][j]) + Msh.u.matB[k];
    }

}


void Discretizer::setPressureBoundaries(Material Mat, Mesh& Msh){

    // Control
    Msh.p.tempA.assign(Msh.p.totNodes, {}); Msh.p.tempB.assign(Msh.p.totNodes, 0);
    int i{}, j{}, k{}; double gammaw{}, gammae{}, gammas{}, gamman{}, a0{}, rho{};
    double Fw{}, Fe{}, Fs{}, Fn{}, Dw{}, De{}, Ds{}, Dn{}, Pw{}, Pe{}, Ps{}, Pn{};

    // Boundary Node Coefficients
    for (boundMain& bC : Msh.boundaryPressure){

        // Boundaries
        if (bC.type == 0){

            // Dirichlet
            if (bC.i0[0] == bC.i1[0]){
                // xBoundary
                if (bC.side == 0){ // West Boundary
                } else if (bC.side == 1){ // East Boundary
                } else {std::cerr << "Boundary side not specified correctly.\n";}
            } else if (bC.i0[1] == bC.i1[1]){
                // yBoundary
                if (bC.side == 0){ // South Boundary
                } else if (bC.side == 1){ // North Boundary
                } else {std::cerr << "Boundary side not specified correctly.\n";}
            } else {std::cerr << "Boundary range not specified correctly.\n";}

        } else if (bC.type == 1){

            // Neumann
            if (bC.i0[0] == bC.i1[0]){

                // xBoundary
                if (bC.side == 0){
                    
                    // West Boundary
                    i = bC.i0[0];
                    for (size_t j = bC.i0[1]; j < bC.i1[1]; j++){
                        // Control
                        k = i * Msh.p.N[1] + j;
                        
                        // Coefficients A 
                        Msh.p.matA[k].aw = bC.value;
                    }

                } else if (bC.side == 1){
                    
                    // East Boundary
                    i = bC.i0[0]-1;
                    for (size_t j = bC.i0[1]; j < bC.i1[1]; j++){
                        // Control
                        k = i * Msh.p.N[1] + j;

                        // Coefficients
                        Msh.p.matA[k].ae = bC.value;
                    }

                } else {std::cerr << "Boundary side not specified correctly.\n";}

            } else if (bC.i0[1] == bC.i1[1]){

                // yBoundary
                if (bC.side == 0){

                    // South Boundary
                    j = bC.i0[1];
                    for (size_t i = bC.i0[0]; i < bC.i1[0]; i++){
                        // Control
                        k = i * Msh.p.N[1] + j;

                        // Coefficients
                        Msh.p.matA[k].as = bC.value;
                    }

                } else if (bC.side == 1){

                    // North Boundary
                    j = bC.i0[1]-1;
                    for (size_t i = bC.i0[0]; i < bC.i1[0]; i++){
                        // Control
                        k = i * Msh.p.N[1] + j;

                        // Coefficients
                        Msh.p.matA[k].an = bC.value;
                    }

                } else {std::cerr << "Boundary side not specified correctly.\n";}

            } else {std::cerr << "Boundary range not specified correctly.\n";}

        } else {std::cerr << "Boundary type not specified correctly.\n";}

    }


    ///// Internal Node Coefficients /////

    // West Boundary (Edges & Corners)
    i = 0;
    for (int j = 0; j < Msh.p.N[1]; j++){
        // Index
        k = i * Msh.p.N[1] + j; rho = Mat.vMat[Msh.p.nMat[i][j]].rho;

        // Corners (~as/~an)
        if (j != 0){Ds = Msh.p.Ss[i][j] / Msh.p.dX[1][j]; Msh.p.matA[k].as = -Ds;} 
        if (j != Msh.p.N[1]-1){Dn = Msh.p.Sn[i][j] / Msh.p.dX[1][j+1]; Msh.p.matA[k].an = -Dn;}
        
        // Coefficients A
        De = Msh.p.Se[i][j] / Msh.p.dX[0][i+1]; Msh.p.matA[k].ae = -De;
        Msh.p.matA[k].ap = - Msh.p.matA[k].aw - Msh.p.matA[k].ae - Msh.p.matA[k].as - Msh.p.matA[k].an;

        // Coefficients B
        Msh.p.matB[k] = (rho / dt) * (Msh.u.Phi[i+1][j] * Msh.p.Se[i][j] + Msh.v.Phi[i][j+1] * Msh.p.Sn[i][j] - Msh.v.Phi[i][j] * Msh.p.Ss[i][j]);
    }

    // East Boundary (Edges & Corners)
    i = Msh.p.N[0]-1;
    for (int j = 0; j < Msh.p.N[1]; j++){
        // Index
        k = i * Msh.p.N[1] + j; rho = Mat.vMat[Msh.p.nMat[i][j]].rho;

        // Corners (~as/~an)
        if (j != 0){Ds = Msh.p.Ss[i][j] / Msh.p.dX[1][j]; Msh.p.matA[k].as = -Ds;}
        if (j != Msh.p.N[1]-1){Dn = Msh.p.Sn[i][j] / Msh.p.dX[1][j+1]; Msh.p.matA[k].an = -Dn;}

        // Coefficients A
        Dw = Msh.p.Sw[i][j] / Msh.p.dX[0][i]; Msh.p.matA[k].aw = -Dw;
        Msh.p.matA[k].ap = - Msh.p.matA[k].aw - Msh.p.matA[k].ae - Msh.p.matA[k].as - Msh.p.matA[k].an;

        // Coefficients B
        Msh.p.matB[k] = (rho / dt) * (- Msh.u.Phi[i][j] * Msh.p.Sw[i][j] + Msh.v.Phi[i][j+1] * Msh.p.Sn[i][j] - Msh.v.Phi[i][j] * Msh.p.Ss[i][j]);
    }
    
    // South Boundary (Edges)
    j = 0;
    for (int i = 1; i < Msh.p.N[0]-1; i++){
        // Index
        k = i * Msh.p.N[1] + j; rho = Mat.vMat[Msh.p.nMat[i][j]].rho;

        // Coefficients A
        Dw = Msh.p.Sw[i][j] / Msh.p.dX[0][i]; Msh.p.matA[k].aw = -Dw;
        De = Msh.p.Se[i][j] / Msh.p.dX[0][i+1]; Msh.p.matA[k].ae = -De;
        Dn = Msh.p.Sn[i][j] / Msh.p.dX[1][j+1]; Msh.p.matA[k].an = -Dn;
        Msh.p.matA[k].ap = - Msh.p.matA[k].aw - Msh.p.matA[k].ae - Msh.p.matA[k].as - Msh.p.matA[k].an;

        // Coefficients B
        Msh.p.matB[k] = (rho / dt) * (Msh.u.Phi[i+1][j] * Msh.p.Se[i][j] - Msh.u.Phi[i][j] * Msh.p.Sw[i][j] + Msh.v.Phi[i][j+1] * Msh.p.Sn[i][j]);
    }

    // North Boundary (Edges)
    j = Msh.p.N[1]-1;
    for (int i = 1; i < Msh.p.N[0]-1; i++){
        // Index
        k = i * Msh.p.N[1] + j; rho = Mat.vMat[Msh.p.nMat[i][j]].rho;

        // Coefficients A
        Dw = Msh.p.Sw[i][j] / Msh.p.dX[0][i]; Msh.p.matA[k].aw = -Dw;
        De = Msh.p.Se[i][j] / Msh.p.dX[0][i+1]; Msh.p.matA[k].ae = -De;
        Ds = Msh.p.Ss[i][j] / Msh.p.dX[1][j]; Msh.p.matA[k].as = -Ds;
        Msh.p.matA[k].ap = - Msh.p.matA[k].aw - Msh.p.matA[k].ae - Msh.p.matA[k].as - Msh.p.matA[k].an;

        // Coefficients B
        Msh.p.matB[k] = (rho / dt) * (Msh.u.Phi[i+1][j] * Msh.p.Se[i][j] - Msh.u.Phi[i][j] * Msh.p.Sw[i][j] - Msh.v.Phi[i][j] * Msh.p.Ss[i][j]);
    }

    /* // Corrections */
    /* Msh.p.matA[0].aw = 0; Msh.p.matA[0].ae = 0; Msh.p.matA[0].as = 0; Msh.p.matA[0].an = 0; */
    /* Msh.p.matA[0].ap = 1; Msh.p.matB[0] = 0; */    


}


void Discretizer::correctVelocity(Material Mat, Mesh& Msh){

    // Control
    double mu = Mat.vMat[Msh.p.nMat[0][0]].mu, rho = Mat.vMat[Msh.p.nMat[0][0]].rho; // simplification because of single material

    // u Interior Nodes
    for (size_t i = 1; i < Msh.u.N[0]-1; i++){
        for (size_t j = 1; j < Msh.u.N[1]-1; j++){
            if (Msh.p.bObs[i-1][j] || Msh.p.bObs[i][j]){continue;}
            Msh.u.Phi[i][j] += (dt / rho) * (Msh.p.Phi[i][j] - Msh.p.Phi[i-1][j]) / Msh.p.dX[0][i];
        }
    }

    // v Interior Nodes
    for (size_t i = 1; i < Msh.v.N[0]-1; i++){
        for (size_t j = 1; j < Msh.v.N[1]-1; j++){
            if (Msh.p.bObs[i][j-1] || Msh.p.bObs[i][j]){continue;}
            Msh.v.Phi[i][j] += (dt / rho) * (Msh.p.Phi[i][j] - Msh.p.Phi[i][j-1]) / Msh.p.dX[1][j];
        }
    }

}


/* void Discretizer::setEnergyBoundaries(Material Mat, Mesh& Msh){ */


/*     // Control */
/*     Msh.T.tempA.assign(Msh.T.totNodes, {}); Msh.T.tempB.assign(Msh.T.totNodes, 0); */
/*     int i{}, j{}, k{}; double gammaw{}, gammae{}, gammas{}, gamman{}, a0{}, rho{}, cp{}, gamma{}; */
/*     double Fw{}, Fe{}, Fs{}, Fn{}, Dw{}, De{}, Ds{}, Dn{}, Pw{}, Pe{}, Ps{}, Pn{}; */
/*     std::vector<int> Pos0{}, Pos1{}; Pos0.resize(Msh.T.N.size()); Pos1.resize(Msh.T.N.size()); */

/*     ////////// Boundary Node Coefficients ////////// */
/*     for (boundMain& bC : Msh.boundaryEnergy){ */

/*         // Positions */
/*         for (size_t i = 0; i < Msh.T.N.size(); i++){ */
/*             Pos0[i] = std::lower_bound(Msh.T.Nodes[i].begin(), Msh.T.Nodes[i].end(), bC.x0[i] - epsFind) - Msh.T.Nodes[i].begin(); */
/*             Pos1[i] = std::lower_bound(Msh.T.Nodes[i].begin(), Msh.T.Nodes[i].end(), bC.x1[i] - epsFind) - Msh.T.Nodes[i].begin(); */
/*         } */

/*         // Boundaries */
/*         if (bC.type == 0){ */

/*             // Dirichlet */
/*             if (Pos0[0] == Pos1[0]){ */

/*                 // xBoundary */
/*                 if (bC.side == 0){ */

/*                     // West Boundary */
/*                     i = Pos0[0]; */
/*                     for (size_t j = Pos0[1]; j < Pos1[1]; j++){ */
/*                         k = i * Msh.T.N[1] + j; */
/*                         gammaw = Mat.vMat[Msh.T.nMat[i][j]].gamma; rho = Mat.vMat[Msh.T.nMat[i][j]].rho; cp = Mat.vMat[Msh.T.nMat[i][j]].cp; */

/*                         // Coefficients Conv-Diff */
/*                         Dw = gammaw * Msh.T.Sw[i][j] / Msh.T.dX[0][i]; Fw = rho * cp * Msh.T.Sw[i][j] * Msh.u.Phi[i][j]; */
/*                         Pw = Fw / Dw; Msh.T.tempA[k].aw = Dw * funcScheme(Pw) + std::max(Fw, 0.0); */

/*                         // Coefficients B */
/*                         Msh.T.tempB[k] += beta * Msh.T.tempA[k].aw * bC.Phi[j] + (1 - beta) * Msh.T.tempA[k].aw * bC.oPhi[j]; */
/*                     } */

/*                 } else if (bC.side == 1){ */

/*                     // East Boundary */
/*                     i = Pos0[0]-1; */
/*                     for (size_t j = Pos0[1]; j < Pos1[1]; j++){ */
/*                         k = i * Msh.T.N[1] + j; */
/*                         gammae = Mat.vMat[Msh.T.nMat[i][j]].gamma; rho = Mat.vMat[Msh.T.nMat[i][j]].rho; cp = Mat.vMat[Msh.T.nMat[i][j]].cp; */

/*                         // Coefficients Conv-Diff */
/*                         De = gammae * Msh.T.Se[i][j] / Msh.T.dX[0][i+1]; Fe = rho * cp * Msh.T.Se[i][j] * Msh.u.Phi[i+1][j]; */
/*                         Pe = Fe / De; Msh.T.tempA[k].ae = De * funcScheme(Pe) + std::max(-Fe, 0.0); */

/*                         // Coefficients B */
/*                         Msh.T.tempB[k] += beta * Msh.T.tempA[k].ae * bC.Phi[j] + (1 - beta) * Msh.T.tempA[k].ae * bC.oPhi[j]; */
/*                     } */

/*                 } else {std::cerr << "Boundary side not specified correctly.\n";} */

/*             } else if (Pos0[1] == Pos1[1]){ */

/*                 // yBoundary */
/*                 if (bC.side == 0){ */

/*                     // South Boundary */
/*                     j = Pos0[1]; */
/*                     for (size_t i = Pos0[0]; i < Pos1[0]; i++){ */
/*                         k = i * Msh.T.N[1] + j; */
/*                         gammas = Mat.vMat[Msh.T.nMat[i][j]].gamma; rho = Mat.vMat[Msh.T.nMat[i][j]].rho; cp = Mat.vMat[Msh.T.nMat[i][j]].cp; */

/*                         // Coefficients Conv-Diff */
/*                         Ds = gammas * Msh.T.Ss[i][j] / Msh.T.dX[1][j]; Fs = rho * cp * Msh.T.Ss[i][j] * Msh.v.Phi[i][j]; */
/*                         Ps = Fs / Ds; Msh.T.tempA[k].as = Ds * funcScheme(Ps) + std::max(Fs, 0.0); */

/*                         // Coefficients B */
/*                         Msh.T.tempB[k] += beta * Msh.T.tempA[k].as * bC.Phi[i] + (1 - beta) * Msh.T.tempA[k].as * bC.oPhi[i]; */
/*                     } */

/*                 } else if (bC.side == 1){ */

/*                     // North Boundary */
/*                     j = Pos0[1]-1; */
/*                     for (size_t i = Pos0[0]; i < Pos1[0]; i++){ */
/*                         k = i * Msh.T.N[1] + j; */
/*                         gamman = Mat.vMat[Msh.T.nMat[i][j]].gamma; rho = Mat.vMat[Msh.T.nMat[i][j]].rho; cp = Mat.vMat[Msh.T.nMat[i][j]].cp; */

/*                         // Coefficients Conv-Diff */
/*                         Dn = gamman * Msh.T.Sn[i][j] / Msh.T.dX[1][j+1]; Fn = rho * cp * Msh.T.Sn[i][j] * Msh.v.Phi[i][j+1]; */
/*                         Pn = Fn / Dn; Msh.T.tempA[k].an = Dn * funcScheme(Pn) + std::max(-Fn, 0.0); */

/*                         // Coefficients B */
/*                         Msh.T.tempB[k] += beta * Msh.T.tempA[k].an * bC.Phi[i] + (1 - beta) * Msh.T.tempA[k].an * bC.oPhi[i]; */
/*                     } */

/*                 } else {std::cerr << "Boundary side not specified correctly.\n";} */

/*             } else {std::cerr << "Boundary range not specified correctly.\n";} */

/*         } else if (bC.type == 1){ */

/*             // Neumann */
/*             if (Pos0[0] == Pos1[0]){ */

/*                 // xBoundary */
/*                 if (bC.side == 0){ */

/*                     // West Boundary */
/*                     i = Pos0[0]; */
/*                     for (size_t j = Pos0[1]; j < Pos1[1]; j++){ */
/*                         k = i * Msh.T.N[1] + j; rho = Mat.vMat[Msh.T.nMat[i][j]].rho; cp = Mat.vMat[Msh.T.nMat[i][j]].cp; */

/*                         // Coefficients A (no diffusion through Neumann face, upwind convection only) */
/*                         Fw = rho * cp * Msh.T.Sw[i][j] * Msh.u.Phi[i][j]; */
/*                         Msh.T.tempA[k].aw = std::max(Fw, 0.0); */

/*                         // Coefficients B */
/*                         Msh.T.tempB[k] += bC.value * Msh.T.Sw[i][j]; */
/*                     } */

/*                 } else if (bC.side == 1){ */

/*                     // East Boundary */
/*                     i = Pos0[0]-1; */
/*                     for (size_t j = Pos0[1]; j < Pos1[1]; j++){ */
/*                         k = i * Msh.T.N[1] + j; rho = Mat.vMat[Msh.T.nMat[i][j]].rho; cp = Mat.vMat[Msh.T.nMat[i][j]].cp; */

/*                         // Coefficients A */
/*                         Fe = rho * cp * Msh.T.Se[i][j] * Msh.u.Phi[i+1][j]; */
/*                         Msh.T.tempA[k].ae = std::max(-Fe, 0.0); */

/*                         // Coefficients B */
/*                         Msh.T.tempB[k] += bC.value * Msh.T.Se[i][j]; */
/*                     } */

/*                 } else {std::cerr << "Boundary side not specified correctly.\n";} */

/*             } else if (Pos0[1] == Pos1[1]){ */

/*                 // yBoundary */
/*                 if (bC.side == 0){ */

/*                     // South Boundary */
/*                     j = Pos0[1]; */
/*                     for (size_t i = Pos0[0]; i < Pos1[0]; i++){ */
/*                         k = i * Msh.T.N[1] + j; rho = Mat.vMat[Msh.T.nMat[i][j]].rho; cp = Mat.vMat[Msh.T.nMat[i][j]].cp; */

/*                         // Coefficients A */
/*                         Fs = rho * cp * Msh.T.Ss[i][j] * Msh.v.Phi[i][j]; */
/*                         Msh.T.tempA[k].as = std::max(Fs, 0.0); */

/*                         // Coefficients B */
/*                         Msh.T.tempB[k] += bC.value * Msh.T.Ss[i][j]; */
/*                     } */

/*                 } else if (bC.side == 1){ */

/*                     // North Boundary */
/*                     j = Pos0[1]-1; */
/*                     for (size_t i = Pos0[0]; i < Pos1[0]; i++){ */
/*                         k = i * Msh.T.N[1] + j; rho = Mat.vMat[Msh.T.nMat[i][j]].rho; cp = Mat.vMat[Msh.T.nMat[i][j]].cp; */

/*                         // Coefficients A */
/*                         Fn = rho * cp * Msh.T.Sn[i][j] * Msh.v.Phi[i][j+1]; */
/*                         Msh.T.tempA[k].an = std::max(-Fn, 0.0); */

/*                         // Coefficients B */
/*                         Msh.T.tempB[k] += bC.value * Msh.T.Sn[i][j]; */
/*                     } */

/*                 } else {std::cerr << "Boundary side not specified correctly.\n";} */

/*             } else {std::cerr << "Boundary range not specified correctly.\n";} */

/*         } else {std::cerr << "Boundary type not specified correctly.\n";} */

/*     } */


/*     ////////// Internal Node Coefficients ////////// */

/*     ///// West Boundary (Edges & Corners) */
/*     i = 0; */
/*     for (int j = 0; j < Msh.T.N[1]; j++){ */

/*         // Index */
/*         k = i * Msh.T.N[1] + j; rho = Mat.vMat[Msh.T.nMat[i][j]].rho; gamma = Mat.vMat[Msh.T.nMat[i][j]].gamma; cp = Mat.vMat[Msh.T.nMat[i][j]].cp; */

/*         // SW Corner (~as) */
/*         if (j != 0){ */
/*             Ds = gamma * Msh.T.Ss[i][j] / Msh.T.dX[1][j]; Fs = rho * cp * Msh.T.Ss[i][j] * Msh.v.Phi[i][j]; */
/*             Ps = Fs / Ds; Msh.T.tempA[k].as = Ds * funcScheme(Ps) + std::max(Fs, 0.0); */
/*             Msh.T.tempB[k] += (1 - beta) * Msh.T.oPhi[i][j-1] * Msh.T.tempA[k].as; */
/*         } */

/*         // NW Corner (~an) */
/*         if (j != Msh.T.N[1]-1){ */
/*             Dn = gamma * Msh.T.Sn[i][j] / Msh.T.dX[1][j+1]; Fn = rho * cp * Msh.T.Sn[i][j] * Msh.v.Phi[i][j+1]; */
/*             Pn = Fn / Dn; Msh.T.tempA[k].an = Dn * funcScheme(Pn) + std::max(-Fn, 0.0); */
/*             Msh.T.tempB[k] += (1 - beta) * Msh.T.oPhi[i][j+1] * Msh.T.tempA[k].an; */
/*         } */

/*         // W Edge (ae) */
/*         De = gamma * Msh.T.Se[i][j] / Msh.T.dX[0][i+1]; Fe = rho * cp * Msh.T.Se[i][j] * Msh.u.Phi[i+1][j]; */
/*         Pe = Fe / De; Msh.T.tempA[k].ae = De * funcScheme(Pe) + std::max(-Fe, 0.0); */

/*         // Coefficients A */
/*         a0 = rho * cp * Msh.T.Vp[i][j] / dt; */
/*         Msh.T.matA[k].aw = -beta * Msh.T.tempA[k].aw; Msh.T.matA[k].ae = -beta * Msh.T.tempA[k].ae; */
/*         Msh.T.matA[k].as = -beta * Msh.T.tempA[k].as; Msh.T.matA[k].an = -beta * Msh.T.tempA[k].an; */
/*         Msh.T.matA[k].ap = a0 - Msh.T.matA[k].aw - Msh.T.matA[k].ae - Msh.T.matA[k].as - Msh.T.matA[k].an; */

/*         // Coefficients B */
/*         Msh.T.matB[k] = Msh.T.sPhi[i][j] * Msh.T.Vp[i][j] + a0 * Msh.T.oPhi[i][j] + (1 - beta) * (Msh.T.oPhi[i+1][j] * Msh.T.tempA[k].ae - Msh.T.oPhi[i][j] * (Msh.T.tempA[k].aw + Msh.T.tempA[k].ae + Msh.T.tempA[k].as + Msh.T.tempA[k].an)) + Msh.T.tempB[k]; */
/*     } */

/*     ///// East Boundary (Edges & Corners) */
/*     i = Msh.T.N[0]-1; */
/*     for (int j = 0; j < Msh.T.N[1]; j++){ */

/*         // Index */
/*         k = i * Msh.T.N[1] + j; rho = Mat.vMat[Msh.T.nMat[i][j]].rho; gamma = Mat.vMat[Msh.T.nMat[i][j]].gamma; cp = Mat.vMat[Msh.T.nMat[i][j]].cp; */

/*         // SE Corner (~as) */
/*         if (j != 0){ */
/*             Ds = gamma * Msh.T.Ss[i][j] / Msh.T.dX[1][j]; Fs = rho * cp * Msh.T.Ss[i][j] * Msh.v.Phi[i][j]; */
/*             Ps = Fs / Ds; Msh.T.tempA[k].as = Ds * funcScheme(Ps) + std::max(Fs, 0.0); */
/*             Msh.T.tempB[k] += (1 - beta) * Msh.T.oPhi[i][j-1] * Msh.T.tempA[k].as; */
/*         } */

/*         // NE Corner (~an) */
/*         if (j != Msh.T.N[1]-1){ */
/*             Dn = gamma * Msh.T.Sn[i][j] / Msh.T.dX[1][j+1]; Fn = rho * cp * Msh.T.Sn[i][j] * Msh.v.Phi[i][j+1]; */
/*             Pn = Fn / Dn; Msh.T.tempA[k].an = Dn * funcScheme(Pn) + std::max(-Fn, 0.0); */
/*             Msh.T.tempB[k] += (1 - beta) * Msh.T.oPhi[i][j+1] * Msh.T.tempA[k].an; */
/*         } */

/*         // E Edge (aw) */
/*         Dw = gamma * Msh.T.Sw[i][j] / Msh.T.dX[0][i]; Fw = rho * cp * Msh.T.Sw[i][j] * Msh.u.Phi[i][j]; */
/*         Pw = Fw / Dw; Msh.T.tempA[k].aw = Dw * funcScheme(Pw) + std::max(Fw, 0.0); */

/*         // Coefficients A */
/*         a0 = rho * cp * Msh.T.Vp[i][j] / dt; */
/*         Msh.T.matA[k].aw = -beta * Msh.T.tempA[k].aw; Msh.T.matA[k].ae = -beta * Msh.T.tempA[k].ae; */
/*         Msh.T.matA[k].as = -beta * Msh.T.tempA[k].as; Msh.T.matA[k].an = -beta * Msh.T.tempA[k].an; */
/*         Msh.T.matA[k].ap = a0 - Msh.T.matA[k].aw - Msh.T.matA[k].ae - Msh.T.matA[k].as - Msh.T.matA[k].an; */

/*         // Coefficients B */
/*         Msh.T.matB[k] = Msh.T.sPhi[i][j] * Msh.T.Vp[i][j] + a0 * Msh.T.oPhi[i][j] + (1 - beta) * (Msh.T.oPhi[i-1][j] * Msh.T.tempA[k].aw - Msh.T.oPhi[i][j] * (Msh.T.tempA[k].aw + Msh.T.tempA[k].ae + Msh.T.tempA[k].as + Msh.T.tempA[k].an)) + Msh.T.tempB[k]; */
/*     } */

/*     ///// South Boundary (Edges) (aw, ae, an, ap) */
/*     j = 0; */
/*     for (int i = 1; i < Msh.T.N[0]-1; i++){ */

/*         // Index */
/*         k = i * Msh.T.N[1] + j; rho = Mat.vMat[Msh.T.nMat[i][j]].rho; gamma = Mat.vMat[Msh.T.nMat[i][j]].gamma; cp = Mat.vMat[Msh.T.nMat[i][j]].cp; */

/*         // Coefficients Conv-Diff */
/*         Dw = gamma * Msh.T.Sw[i][j] / Msh.T.dX[0][i]; Fw = rho * cp * Msh.T.Sw[i][j] * Msh.u.Phi[i][j]; */
/*         Pw = Fw / Dw; Msh.T.tempA[k].aw = Dw * funcScheme(Pw) + std::max(Fw, 0.0); */
/*         De = gamma * Msh.T.Se[i][j] / Msh.T.dX[0][i+1]; Fe = rho * cp * Msh.T.Se[i][j] * Msh.u.Phi[i+1][j]; */
/*         Pe = Fe / De; Msh.T.tempA[k].ae = De * funcScheme(Pe) + std::max(-Fe, 0.0); */
/*         Dn = gamma * Msh.T.Sn[i][j] / Msh.T.dX[1][j+1]; Fn = rho * cp * Msh.T.Sn[i][j] * Msh.v.Phi[i][j+1]; */
/*         Pn = Fn / Dn; Msh.T.tempA[k].an = Dn * funcScheme(Pn) + std::max(-Fn, 0.0); */

/*         // Coefficients A */
/*         a0 = rho * cp * Msh.T.Vp[i][j] / dt; */
/*         Msh.T.matA[k].aw = -beta * Msh.T.tempA[k].aw; Msh.T.matA[k].ae = -beta * Msh.T.tempA[k].ae; */
/*         Msh.T.matA[k].as = -beta * Msh.T.tempA[k].as; Msh.T.matA[k].an = -beta * Msh.T.tempA[k].an; */
/*         Msh.T.matA[k].ap = a0 - Msh.T.matA[k].aw - Msh.T.matA[k].ae - Msh.T.matA[k].as - Msh.T.matA[k].an; */

/*         // Coefficients B */
/*         Msh.T.matB[k] = Msh.T.sPhi[i][j] * Msh.T.Vp[i][j] + a0 * Msh.T.oPhi[i][j] + (1 - beta) * (Msh.T.oPhi[i-1][j] * Msh.T.tempA[k].aw + Msh.T.oPhi[i+1][j] * Msh.T.tempA[k].ae + Msh.T.oPhi[i][j+1] * Msh.T.tempA[k].an - Msh.T.oPhi[i][j] * (Msh.T.tempA[k].aw + Msh.T.tempA[k].ae + Msh.T.tempA[k].as + Msh.T.tempA[k].an)) + Msh.T.tempB[k]; */
/*     } */

/*     ///// North Boundary (Edges) (aw, ae, as, ap) */
/*     j = Msh.T.N[1]-1; */
/*     for (int i = 1; i < Msh.T.N[0]-1; i++){ */

/*         // Index */
/*         k = i * Msh.T.N[1] + j; rho = Mat.vMat[Msh.T.nMat[i][j]].rho; gamma = Mat.vMat[Msh.T.nMat[i][j]].gamma; cp = Mat.vMat[Msh.T.nMat[i][j]].cp; */

/*         // Coefficients Conv-Diff */
/*         Dw = gamma * Msh.T.Sw[i][j] / Msh.T.dX[0][i]; Fw = rho * cp * Msh.T.Sw[i][j] * Msh.u.Phi[i][j]; */
/*         Pw = Fw / Dw; Msh.T.tempA[k].aw = Dw * funcScheme(Pw) + std::max(Fw, 0.0); */
/*         De = gamma * Msh.T.Se[i][j] / Msh.T.dX[0][i+1]; Fe = rho * cp * Msh.T.Se[i][j] * Msh.u.Phi[i+1][j]; */
/*         Pe = Fe / De; Msh.T.tempA[k].ae = De * funcScheme(Pe) + std::max(-Fe, 0.0); */
/*         Ds = gamma * Msh.T.Ss[i][j] / Msh.T.dX[1][j]; Fs = rho * cp * Msh.T.Ss[i][j] * Msh.v.Phi[i][j]; */
/*         Ps = Fs / Ds; Msh.T.tempA[k].as = Ds * funcScheme(Ps) + std::max(Fs, 0.0); */

/*         // Coefficients A */
/*         a0 = rho * cp * Msh.T.Vp[i][j] / dt; */
/*         Msh.T.matA[k].aw = -beta * Msh.T.tempA[k].aw; Msh.T.matA[k].ae = -beta * Msh.T.tempA[k].ae; */
/*         Msh.T.matA[k].as = -beta * Msh.T.tempA[k].as; Msh.T.matA[k].an = -beta * Msh.T.tempA[k].an; */
/*         Msh.T.matA[k].ap = a0 - Msh.T.matA[k].aw - Msh.T.matA[k].ae - Msh.T.matA[k].as - Msh.T.matA[k].an; */

/*         // Coefficients B */
/*         Msh.T.matB[k] = Msh.T.sPhi[i][j] * Msh.T.Vp[i][j] + a0 * Msh.T.oPhi[i][j] + (1 - beta) * (Msh.T.oPhi[i-1][j] * Msh.T.tempA[k].aw + Msh.T.oPhi[i+1][j] * Msh.T.tempA[k].ae + Msh.T.oPhi[i][j-1] * Msh.T.tempA[k].as - Msh.T.oPhi[i][j] * (Msh.T.tempA[k].aw + Msh.T.tempA[k].ae + Msh.T.tempA[k].as + Msh.T.tempA[k].an)) + Msh.T.tempB[k]; */
/*     } */

/* } */


/* void Discretizer::setEnergyCoefficients(Material Mat, Mesh& Msh){ */

/*     // Control */
/*     int k{}; */
/*     double Dw{}, De{}, Ds{}, Dn{}, a0{}, Fw{}, Fe{}, Fs{}, Fn{}, Pw{}, Pe{}, Ps{}, Pn{}; */
/*     double gamma = Mat.vMat[Msh.T.nMat[0][0]].gamma, rho = Mat.vMat[Msh.T.nMat[0][0]].rho, cp = Mat.vMat[Msh.T.nMat[0][0]].cp; // simplification because of single material */

/*     // p Interior Nodes */
/*     for (size_t i = 1; i < Msh.T.N[0]-1; i++){ */
/*         for (size_t j = 1; j < Msh.T.N[1]-1; j++){ */
/*             // Control */
/*             k = i * Msh.T.N[1] + j; a0 = rho * cp * Msh.T.Vp[i][j] / dt; */

/*             // Coefficients Convection-Diffusion */
/*             Dw = gamma * Msh.T.Sw[i][j] / Msh.T.dX[0][i]; Fw = rho * cp * Msh.T.Sw[i][j] * Msh.u.Phi[i][j]; */
/*             Pw = Fw / Dw; Msh.T.tempA[k].aw = Dw * funcScheme(Pw) + std::max(Fw, 0.0); */
/*             De = gamma * Msh.T.Se[i][j] / Msh.T.dX[0][i+1]; Fe = rho * cp * Msh.T.Se[i][j] * Msh.u.Phi[i+1][j]; */
/*             Pe = Fe / De; Msh.T.tempA[k].ae = De * funcScheme(Pe) + std::max(-Fe, 0.0); */
/*             Ds = gamma * Msh.T.Ss[i][j] / Msh.T.dX[1][j];  Fs = rho * cp * Msh.T.Ss[i][j] * Msh.v.Phi[i][j]; */
/*             Ps = Fs / Ds; Msh.T.tempA[k].as = Ds * funcScheme(Ps) + std::max(Fs, 0.0); */
/*             Dn = gamma * Msh.T.Sn[i][j] / Msh.T.dX[1][j+1]; Fn = rho * cp * Msh.T.Sn[i][j] * Msh.v.Phi[i][j+1]; */
/*             Pn = Fn / Dn; Msh.T.tempA[k].an = Dn * funcScheme(Pn) + std::max(-Fn, 0.0); */

/*             // Coefficients A */
/*             Msh.T.matA[k].aw = - beta * Msh.T.tempA[k].aw; Msh.T.matA[k].ae = - beta * Msh.T.tempA[k].ae; */
/*             Msh.T.matA[k].as = - beta * Msh.T.tempA[k].as; Msh.T.matA[k].an = - beta * Msh.T.tempA[k].an; */
/*             Msh.T.matA[k].ap = a0 - Msh.T.matA[k].aw - Msh.T.matA[k].ae - Msh.T.matA[k].as - Msh.T.matA[k].an; */
    
/*             // Coefficients B */
/*             Msh.T.matB[k] = Msh.T.sPhi[i][j] * Msh.T.Vp[i][j] + a0 * Msh.T.oPhi[i][j] + (1 - beta) * (Msh.T.oPhi[i-1][j] * Msh.T.tempA[k].aw + Msh.T.oPhi[i+1][j] * Msh.T.tempA[k].ae + Msh.T.oPhi[i][j-1] * Msh.T.tempA[k].as + Msh.T.oPhi[i][j+1] * Msh.T.tempA[k].an - Msh.T.oPhi[i][j] * (Msh.T.tempA[k].ae + Msh.T.tempA[k].aw + Msh.T.tempA[k].an + Msh.T.tempA[k].as)); */
/*         } */
/*     } */

/* } */
