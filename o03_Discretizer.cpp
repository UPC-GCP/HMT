// Imports
/* #include "json/value.h" */
#include <csignal>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
/* #include <iterator> */
#include <sys/ucontext.h>
#include <vector>
#include <json/json.h>
#include <cmath>
#include <algorithm>

// Self-Imports
/* #include "exprtk.hpp" */
#include "o01_Material.h"
#include "o02_Mesh.h"
#include "o03_Discretizer.h"
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
        funcScheme = [](double Pe) {return 1.0 - 0.5 * Pe;};
    } else if (spatScheme == "UDS"){
        funcScheme = [](double Pe){return 1;};
    } else if (spatScheme == "Hybrid"){
        funcScheme = [](double Pe){return std::max(0.0, 1 -  0.5 * Pe);};
    } else if (spatScheme == "PowerLaw"){
        funcScheme = [](double Pe){return std::max(0.0, std::pow(1 - 0.1 * Pe, 5));};
    } else if (spatScheme == "Exponential"){
        funcScheme = [](double Pe){return Pe / (exp(Pe) - 1);};
    } else {
        std::cerr << "Error: Invalid spatial discretization scheme '" << spatScheme << "'\n";
    }

}


void Discretizer::newSetBoundaries(Material& Mat, Mesh& Msh, ExpressionParser& Prs, double t){

    // Control
    Msh.tempB.resize(Msh.totNodes, 0);
    int i{}, j{}, k{}; double gammaw, gammae, gammas{}, gamman{};
    double Fw{}, Fe{}, Fs{}, Fn{}, Dw{}, De{}, Ds{}, Dn{}, Pw{}, Pe{}, Ps{}, Pn{};
    std::vector<int> Pos0{}, Pos1{}; Pos0.resize(Msh.N.size()); Pos1.resize(Msh.N.size());

    ////////// Boundary Node Coefficients //////////
    for (Boundary bC : Msh.boundaryConditions){

        // Positions
        for (size_t i = 0; i < Msh.N.size(); i++){
            Pos0[i] = std::lower_bound(Msh.Nodes[i].begin(), Msh.Nodes[i].end(), bC.x0[i] - epsFind) - Msh.Nodes[i].begin();
            Pos1[i] = std::lower_bound(Msh.Nodes[i].begin(), Msh.Nodes[i].end(), bC.x1[i] - epsFind) - Msh.Nodes[i].begin();
        }

        // Boundaries (Non-nD)
        if (bC.type == 0){

            // Update Value 
            if (bC.bUpdate && bC.iEq == 0){bC.value = Prs.evaluateTime(bC.iExpr, t);}

            // Dirichlet
            if (Pos0[0] == Pos1[0]){

                // xBoundary
                if (bC.side == 0){

                    // West Boundary
                    i = Pos0[0];

                    for (size_t j = Pos0[1]; j < Pos1[1]; j++){
                        // Control
                        k = i * Msh.N[1] + j;
                        if (bC.bUpdate && bC.iEq == 1){bC.value = Prs.evaluateCoordinates(bC.iExpr, Msh.Faces[0][i], Msh.Nodes[1][j]);}
                        
                        // Gamma
                        gammaw = Mat.vMat[Msh.nMat[i][j]].gamma;

                        // Coefficients Conv-Diff
                        Dw = gammaw * Msh.Sw[i][j] / Msh.dX[0][i]; Fw = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sw[i][j] * Msh.vConv[0][i][j].Vw;
                        Pw = Fw / Dw; Msh.tempA[k].aw = Dw * funcScheme(std::abs(Pw)) + std::max(-Fw, 0.0);

                        // Coefficients B
                        Msh.tempB[k] += - (1 - beta) * (Msh.tempA[k].aw * bC.value);
                    }

                } else if (bC.side == 1){
                    
                    // East Boundary
                    i = Pos0[0] - 1;

                    for (size_t j = Pos0[1]; j < Pos1[1]; j++){
                        // Control
                        k = i * Msh.N[1] + j;
                        if (bC.bUpdate && bC.iEq == 1){bC.value = Prs.evaluateCoordinates(bC.iExpr, Msh.Faces[0][i+1], Msh.Nodes[1][j]);}

                        // Gamma
                        gammae = Mat.vMat[Msh.nMat[i][j]].gamma;

                        // Coefficients Conv-Diff
                        De = gammae * Msh.Se[i][j] / Msh.dX[0][i+1]; Fe = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Se[i][j] * Msh.vConv[0][i][j].Ve;
                        Pe = Fe / De; Msh.tempA[k].ae = De * funcScheme(std::abs(Pe)) + std::max(Fe, 0.0);

                        // Coefficients B
                        Msh.tempB[k] += - (1 - beta) * (Msh.tempA[k].ae * bC.value);
                    }

                } else {std::cerr << "Boundary side not specified correctly.\n";}

            } else if (Pos0[1] == Pos1[1]){
    
                // yBoundary
                if (bC.side == 0){

                    // South Boundary
                    j = Pos0[1];

                    for (size_t i = Pos0[0]; i < Pos1[0]; i++){
                        // Control
                        k = i * Msh.N[1] + j;
                        if (bC.bUpdate && bC.iEq == 1){bC.value = Prs.evaluateCoordinates(bC.iExpr, Msh.Nodes[0][i], Msh.Faces[1][j]);}

                        // Gamma
                        gammas = Mat.vMat[Msh.nMat[i][j]].gamma;

                        // Coefficients Conv-Diff
                        Ds = gammas * Msh.Ss[i][j] / Msh.dX[1][j]; Fs = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Ss[i][j] * Msh.vConv[1][i][j].Vs;
                        Ps = Fs / Ds; Msh.tempA[k].as = Ds * funcScheme(std::abs(Ps)) + std::max(-Fs, 0.0);

                        // Coefficients B
                        Msh.tempB[k] += - (1 - beta) * (Msh.tempA[k].as * bC.value);
                    }

                } else if (bC.side == 1){

                    // North Boundary
                    j = Pos0[1] - 1;

                    for (size_t i = Pos0[0]; i < Pos1[0]; i++){
                        // North Boundary
                        k = i * Msh.N[1] + j;
                        if (bC.bUpdate && bC.iEq == 1){bC.value = Prs.evaluateCoordinates(bC.iExpr, Msh.Nodes[0][i], Msh.Faces[1][j+1]);}

                        // Gamma
                        gamman = Mat.vMat[Msh.nMat[i][j]].gamma;

                        // Coefficients Conv-Diff
                        Dn = gamman * Msh.Sn[i][j] / Msh.dX[1][j+1]; Fn = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sn[i][j] * Msh.vConv[1][i][j].Vn;
                        Pn = Fn / Dn; Msh.tempA[k].an = Dn * funcScheme(std::abs(Pn)) + std::max(Fn, 0.0);

                        // Coefficients B
                        Msh.tempB[k] += - (1 - beta) * (Msh.tempA[k].an * bC.value);
                    }

                } else {std::cerr << "Boundary range not specified correctly.\n";}

            }
            
        } else if (bC.type == 1){

            // Neumann
            if (Pos0[0] == Pos1[0]){

                // Control
                double phiBC{}, aTemp{};

                // xBoundary
                if (bC.side == 0){

                    // West Boundary
                    i = Pos0[0]; 

                    for (size_t j = Pos0[1]; j < Pos1[1]; j++){
                        // Control
                        k = i * Msh.N[1] + j; gammaw = Mat.vMat[Msh.nMat[i][j]].gamma;

                        // Coefficients A
                        Msh.tempA[k].aw = 0;
                        Dw = gammaw * Msh.Sw[i][j] / Msh.dX[0][i]; Fw = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sw[i][j] * Msh.vConv[0][i][j].Vw;
                        Pw = Fw / Dw; aTemp = Dw * funcScheme(std::abs(Pw)) + std::max(-Fw, 0.0);

                        // Coefficients B
                        phiBC = (bC.value + Msh.vPhi[i][j] * gammaw / Msh.dX[0][i]) / (gammaw / Msh.dX[0][i]);
                        Msh.tempB[k] += bC.value * Msh.Sw[i][j] - (1 - beta) * (phiBC * aTemp - Msh.vPhi[i][j] * aTemp);
                    }

                } else if (bC.side == 1){
                    
                    // East Boundary
                    i = Pos0[0]-1; 

                    for (size_t j = Pos0[1]; j < Pos1[1]; j++){
                        // Control
                        k = i * Msh.N[1] + j;

                        // Coefficients A
                        Msh.tempA[k].ae = 0;
                        De = gammae * Msh.Se[i][j] / Msh.dX[0][i+1]; Fe = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Se[i][j] * Msh.vConv[0][i][j].Ve;
                        Pe = Fe / De; aTemp = De * funcScheme(std::abs(Pe)) + std::max(Fe, 0.0);

                        // Coefficients B
                        phiBC = (bC.value - Msh.vPhi[i][j] * gammae / Msh.dX[0][i+1]) / (gammae / Msh.dX[0][i+1]); // CHECK SIGN, NOT SURE (PENDING)
                        Msh.tempB[k] += bC.value * Msh.Se[i][j] - (1 - beta) * (phiBC * aTemp - Msh.vPhi[i][j] * aTemp);
                    }

                } else {std::cerr << "Boundary side not specified correctly.\n";}
                
            } else if (Pos0[1] == Pos1[1]){

                // Control
                double phiBC{}, aTemp{};

                // yBoundary
                if (bC.side == 0){
                    
                    // South Boundary
                    j = Pos0[1];

                    for (size_t i = Pos0[0]; i < Pos1[0]; i++){
                        // Control
                        k = i * Msh.N[1] + j;

                        // Coefficients A
                        Msh.tempA[k].as = 0;
                        Ds = gammas * Msh.Ss[i][j] / Msh.dX[1][j]; Fs = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Ss[i][j] * Msh.vConv[1][i][j].Vs;
                        Ps = Fs / Ds; aTemp = Ds * funcScheme(std::abs(Ps)) + std::max(-Fs, 0.0);

                        // Coefficients B
                        phiBC = (bC.value + Msh.vPhi[i][j] * gammas / Msh.dX[1][j]) / (gammas / Msh.dX[1][j]);
                        Msh.tempB[k] += bC.value * Msh.Ss[i][j] - (1 - beta) * (phiBC * aTemp - Msh.vPhi[i][j] * aTemp);
                    }
                    
                } else if (bC.side == 1){

                    // North Boundary
                    j = Pos0[1]-1;

                    for (size_t i = Pos0[0]; i < Pos1[0]; i++){
                        // Control
                        k = i * Msh.N[1] + j;

                        // Coefficients A
                        Msh.tempA[k].an = 0;
                        Dn = gamman * Msh.Sn[i][j] / Msh.dX[1][j+1]; Fn = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sn[i][j] * Msh.vConv[1][i][j].Vn;
                        Pn = Fn / Dn; aTemp = Dn * funcScheme(std::abs(Pn)) + std::max(Fn, 0.0);

                        // Coefficients B
                        phiBC = (bC.value - Msh.vPhi[i][j] * gamman / Msh.dX[1][j+1]) / (gamman / Msh.dX[1][j+1]);
                        Msh.tempB[k] += bC.value * Msh.Sn[i][j] - (1 - beta) * (phiBC * aTemp - Msh.vPhi[i][j] * aTemp);
                    }

                } else {std::cerr << "Boundary side not specified correctly.\n";}
                
            } else {std::cerr << "Boundary range not specified correctly.\n";}


        } else if (bC.type == 2){ // PENDING THIS SHOULD BE IMPLEMENTED AS PURE DIFFUSION

            // Robin (Hybrid)
            if (Pos0[0] == Pos1[0]){

                // xBoundary

            } else if (Pos0[1] == Pos1[1]){

                // yBoundary

            } else {std::cerr << "Boundary range not specified correctly.\n";}

        } else {std::cerr << "Boundary type not specified correctly.\n";}

    }

    
    ////////// Internal Node Coefficients //////////
    
    ///// West Boundary (Edges & Corners)
    i = 0; 
    for (int j = 0; j < Msh.N[1]; j++){

        // Index
        k = i * Msh.N[1] + j;

        // SW Corner (~as)
        if (j != 0){
            gammas = calcHarmonicMean(Msh.dX[1][j], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i][j-1]].gamma}, {Msh.DeltaX[1][j], Msh.DeltaX[1][j-1]});
            Ds = gammas * Msh.Ss[i][j] / Msh.dX[1][j]; Fs = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Ss[i][j] * Msh.vConv[1][i][j].Vs;
            Ps = Fs / Ds; Msh.tempA[k].as = Ds * funcScheme(std::abs(Ps)) + std::max(-Fs, 0.0);
            Msh.tempB[k] += - (1 - beta) * Msh.vPhi[i][j-1] * Msh.tempA[k].as;
        }

        // NW Corner (~an)
        if (j != Msh.N[1]-1){
            gamman = calcHarmonicMean(Msh.dX[1][j+1], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i][j+1]].gamma}, {Msh.DeltaX[1][j], Msh.DeltaX[1][j+1]});
            Dn = gamman * Msh.Sn[i][j] / Msh.dX[1][j+1]; Fn = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sn[i][j] * Msh.vConv[1][i][j].Vn;
            Pn = Fn / Dn; Msh.tempA[k].an = Dn * funcScheme(std::abs(Pn)) + std::max(Fn, 0.0);
            Msh.tempB[k] += - (1 - beta) * Msh.vPhi[i][j+1] * Msh.tempA[k].an;
        } 

        // W Edge (ae)
        gammae = calcHarmonicMean(Msh.dX[0][i+1], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i+1][j]].gamma}, {Msh.DeltaX[0][i], Msh.DeltaX[0][i+1]});
        De = gammae * Msh.Se[i][j] / Msh.dX[0][i+1]; Fe = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Se[i][j] * Msh.vConv[0][i][j].Ve;
        Pe = Fe / De; Msh.tempA[k].ae = De * funcScheme(std::abs(Pe)) + std::max(Fe, 0.0);
       
        // Coefficients A
        Msh.matA[k].aw = - beta * Msh.tempA[k].aw; Msh.matA[k].ae = - beta * Msh.tempA[k].ae;
        Msh.matA[k].as = - beta * Msh.tempA[k].as; Msh.matA[k].an = - beta * Msh.tempA[k].an; 
        Msh.matA[k].ap = Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp * Msh.Vp[i][j] / dt - Msh.matA[k].ae - Msh.matA[k].aw - Msh.matA[k].an - Msh.matA[k].as;
            
        // Coefficients B 
        // Msh.tempB[k] = Msh.vPhi[i-1][j] * Msh.tempA[k].aw + Msh.vPhi[i][j-1] * Msh.tempA[k].as + Msh.vPhi[i][j+1] * Msh.tempA[k].an
        // No Msh.vPhi[i-1][j] so Msh.tempB[k] needs to be calculated at boundary
        // No Msh.vPhi[i][j+1] || Msh.vPhi[i][j-1] for corners so Msh.tempB[k] needs to be added with each term
        
        Msh.bp[k] = Msh.sPhi[i][j] * Msh.Vp[i][j] + Mat.vMat[Msh.nMat[i][j]].rho * Msh.Vp[i][j] * Msh.vPhi[i][j] / dt - (1 - beta) * (Msh.vPhi[i+1][j] * Msh.tempA[k].ae - Msh.vPhi[i][j] * (Msh.tempA[k].ae + Msh.tempA[k].aw + Msh.tempA[k].an + Msh.tempA[k].as)) + Msh.tempB[k];

    }

    ///// East Boundary (Edges & Corners)
    i = Msh.N[0]-1; 
    for (int j = 0; j < Msh.N[1]; j++){
      
        // Index
        k = i * Msh.N[1] + j;

        // SE Corner (~as)
        if (j != 0){
            gammas = calcHarmonicMean(Msh.dX[1][j], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i][j-1]].gamma}, {Msh.DeltaX[1][j], Msh.DeltaX[1][j-1]});
            Ds = gammas * Msh.Ss[i][j] / Msh.dX[1][j]; Fs = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Ss[i][j] * Msh.vConv[1][i][j].Vs;
            Ps = Fs / Ds; Msh.tempA[k].as = Ds * funcScheme(std::abs(Ps)) + std::max(-Fs, 0.0);
            Msh.tempB[k] += - (1 - beta) * Msh.vPhi[i][j-1] * Msh.tempA[k].as;
        }

        // NE Corner (~an)
        if (j != Msh.N[1]-1){
            gamman = calcHarmonicMean(Msh.dX[1][j+1], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i][j+1]].gamma}, {Msh.DeltaX[1][j], Msh.DeltaX[1][j+1]});
            Dn = gamman * Msh.Sn[i][j] / Msh.dX[1][j+1]; Fn = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sn[i][j] * Msh.vConv[1][i][j].Vn;
            Pn = Fn / Dn; Msh.tempA[k].an = Dn * funcScheme(std::abs(Pn)) + std::max(Fn, 0.0);
            Msh.tempB[k] += - (1 - beta) * Msh.vPhi[i][j+1] * Msh.tempA[k].an;
        } 

        // East Edge (aw)
        Dw = gammaw * Msh.Sw[i][j] / Msh.dX[0][i]; Fw = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sw[i][j] * Msh.vConv[0][i][j].Vw;
        Pw = Fw / Dw; Msh.tempA[k].aw = Dw * funcScheme(std::abs(Pw)) + std::max(-Fw, 0.0);

        // Coefficients A
        Msh.matA[k].aw = - beta * Msh.tempA[k].aw; Msh.matA[k].ae = - beta * Msh.tempA[k].ae;
        Msh.matA[k].as = - beta * Msh.tempA[k].as; Msh.matA[k].an = - beta * Msh.tempA[k].an; 
        Msh.matA[k].ap = Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp * Msh.Vp[i][j] / dt - Msh.matA[k].ae - Msh.matA[k].aw - Msh.matA[k].an - Msh.matA[k].as;
            
        // Coefficients B
        // Msh.tempB[k] = Msh.vPhi[i+1][j] * Msh.tempA[k].ae + Msh.vPhi[i][j-1] * Msh.tempA[k].as + Msh.vPhi[i][j+1] * Msh.tempA[k].an
        // No Msh.vPhi[i+1][j] so Msh.tempB[k] needs to be calculated at boundary
        // No Msh.vPhi[i][j+1] || Msh.vPhi[i][j-1] for corners so Msh.tempB[k] needs to be added with each term
        Msh.bp[k] = Msh.sPhi[i][j] * Msh.Vp[i][j] + Mat.vMat[Msh.nMat[i][j]].rho * Msh.Vp[i][j] * Msh.vPhi[i][j] / dt - (1 - beta) * (Msh.vPhi[i-1][j] * Msh.tempA[k].aw - Msh.vPhi[i][j] * (Msh.tempA[k].ae + Msh.tempA[k].aw + Msh.tempA[k].an + Msh.tempA[k].as)) + Msh.tempB[k];
  
    }

    ///// South Boundary (Edges) (aw, ae, an, ap)
    j = 0;
    for (int i = 1; i < Msh.N[0]-1; i++){

        // Index
        k = i * Msh.N[1] + j;

        // Harmonic Mean
        gammaw = calcHarmonicMean(Msh.dX[0][i], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i-1][j]].gamma}, {Msh.DeltaX[0][i], Msh.DeltaX[0][i-1]});
        gammae = calcHarmonicMean(Msh.dX[0][i+1], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i+1][j]].gamma}, {Msh.DeltaX[0][i], Msh.DeltaX[0][i+1]});
        gamman = calcHarmonicMean(Msh.dX[1][j+1], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i][j+1]].gamma}, {Msh.DeltaX[1][j], Msh.DeltaX[1][j+1]});

        // Coefficients Conv-Diff
        Dw = gammaw * Msh.Sw[i][j] / Msh.dX[0][i]; Fw = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sw[i][j] * Msh.vConv[0][i][j].Vw;
        Pw = Fw / Dw; Msh.tempA[k].aw = Dw * funcScheme(std::abs(Pw)) + std::max(-Fw, 0.0);
        De = gammae * Msh.Se[i][j] / Msh.dX[0][i+1]; Fe = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Se[i][j] * Msh.vConv[0][i][j].Ve;
        Pe = Fe / De; Msh.tempA[k].ae = De * funcScheme(std::abs(Pe)) + std::max(Fe, 0.0);
        Dn = gamman * Msh.Sn[i][j] / Msh.dX[1][j+1]; Fn = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sn[i][j] * Msh.vConv[1][i][j].Vn;
        Pn = Fn / Dn; Msh.tempA[k].an = Dn * funcScheme(std::abs(Pn)) + std::max(Fn, 0.0);
        
        // Coefficients A - tempA keeps a_k values without beta
        Msh.matA[k].aw = - beta * Msh.tempA[k].aw; Msh.matA[k].ae = - beta * Msh.tempA[k].ae;
        Msh.matA[k].as = - beta * Msh.tempA[k].as; Msh.matA[k].an = - beta * Msh.tempA[k].an; 
        Msh.matA[k].ap = Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp * Msh.Vp[i][j] / dt - Msh.matA[k].ae - Msh.matA[k].aw - Msh.matA[k].an - Msh.matA[k].as;
            
        // Coefficients B - tempB keeps boundary values until needed PENDIENTE CAMBIAR A SUM TEMPB
        // Msh.tempB[k] = Msh.vPhi[i][j-1] * Msh.tempA[k].as
        // No Msh.vPhi[i][j-1] so Msh.tempB[k] needs to be calculated at boundary
        Msh.bp[k] = Msh.sPhi[i][j] * Msh.Vp[i][j] + Mat.vMat[Msh.nMat[i][j]].rho * Msh.Vp[i][j] * Msh.vPhi[i][j] / dt - (1 - beta) * (Msh.vPhi[i+1][j] * Msh.tempA[k].ae + Msh.vPhi[i-1][j] * Msh.tempA[k].aw + Msh.vPhi[i][j+1] * Msh.tempA[k].an - Msh.vPhi[i][j] * (Msh.tempA[k].ae + Msh.tempA[k].aw + Msh.tempA[k].an + Msh.tempA[k].as)) + Msh.tempB[k];

    }

    ///// North Boundary (Edges) (aw, ae, as, ap)
    j = Msh.N[1]-1;
    for (int i = 1; i < Msh.N[0]-1; i++){
        
        // Index
        k = i * Msh.N[1] + j;

        // Harmonic Mean
        gammaw = calcHarmonicMean(Msh.dX[0][i], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i-1][j]].gamma}, {Msh.DeltaX[0][i], Msh.DeltaX[0][i-1]});
        gammae = calcHarmonicMean(Msh.dX[0][i+1], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i+1][j]].gamma}, {Msh.DeltaX[0][i], Msh.DeltaX[0][i+1]});
        gammas = calcHarmonicMean(Msh.dX[1][j], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i][j-1]].gamma},{Msh.DeltaX[1][j], Msh.DeltaX[1][j-1]});

        // Coefficients Conv-Diff
        Dw = gammaw * Msh.Sw[i][j] / Msh.dX[0][i]; Fw = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sw[i][j] * Msh.vConv[0][i][j].Vw;
        Pw = Fw / Dw; Msh.tempA[k].aw = Dw * funcScheme(std::abs(Pw)) + std::max(-Fw, 0.0);
        De = gammae * Msh.Se[i][j] / Msh.dX[0][i+1]; Fe = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Se[i][j] * Msh.vConv[0][i][j].Ve;
        Pe = Fe / De; Msh.tempA[k].ae = De * funcScheme(std::abs(Pe)) + std::max(Fe, 0.0);
        Ds = gammas * Msh.Ss[i][j] / Msh.dX[1][j]; Fs = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Ss[i][j] * Msh.vConv[1][i][j].Vs;
        Ps = Fs / Ds; Msh.tempA[k].as = Ds * funcScheme(std::abs(Ps)) + std::max(-Fs, 0.0);
        
        // Coefficients A
        Msh.matA[k].aw = - beta * Msh.tempA[k].aw; Msh.matA[k].ae = - beta * Msh.tempA[k].ae;
        Msh.matA[k].as = - beta * Msh.tempA[k].as; Msh.matA[k].an = - beta * Msh.tempA[k].an; 
        Msh.matA[k].ap = Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp * Msh.Vp[i][j] / dt - Msh.matA[k].ae - Msh.matA[k].aw - Msh.matA[k].an - Msh.matA[k].as;
            
        // Coefficients B
        // Msh.tempB[k] = Msh.vPhi[i][j+1] * Msh.tempA[k].an
        // No Msh.vPhi[i][j+1] so Msh.tempB[k] needs to be calculated at boundary
        Msh.bp[k] = Msh.sPhi[i][j] * Msh.Vp[i][j] + Mat.vMat[Msh.nMat[i][j]].rho * Msh.Vp[i][j] * Msh.vPhi[i][j] / dt - (1 - beta) * (Msh.vPhi[i+1][j] * Msh.tempA[k].ae + Msh.vPhi[i-1][j] * Msh.tempA[k].aw + Msh.vPhi[i][j-1] * Msh.tempA[k].as - Msh.vPhi[i][j] * (Msh.tempA[k].ae + Msh.tempA[k].aw + Msh.tempA[k].an + Msh.tempA[k].as)) + Msh.tempB[k];

    }

}


void Discretizer::newSetCoefficients(Material& Mat, Mesh& Msh){

    // Control
    int k{}; double gammaw, gammae, gammas{}, gamman{};
    double Fw{}, Fe{}, Fs{}, Fn{}, Dw{}, De{}, Ds{}, Dn{}, Pw{}, Pe{}, Ps{}, Pn{};

    // Interior Nodes Loop
    for (size_t i = 1; i < Msh.N[0]-1; i++){
        for (size_t j = 1; j < Msh.N[1]-1; j++){
        
            // Harmonic Mean
            gammaw = calcHarmonicMean(Msh.dX[0][i], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i-1][j]].gamma}, {Msh.DeltaX[0][i], Msh.DeltaX[0][i-1]});
            gammae = calcHarmonicMean(Msh.dX[0][i+1], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i+1][j]].gamma}, {Msh.DeltaX[0][i], Msh.DeltaX[0][i+1]});
            gammas = calcHarmonicMean(Msh.dX[1][j], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i][j-1]].gamma}, {Msh.DeltaX[1][j], Msh.DeltaX[1][j-1]});
            gamman = calcHarmonicMean(Msh.dX[1][j+1], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i][j+1]].gamma}, {Msh.DeltaX[1][j], Msh.DeltaX[1][j+1]});

            // Index
            k = i * Msh.N[1] + j;

            // Coefficients Conv-Diff
            Dw = gammaw * Msh.Sw[i][j] / Msh.dX[0][i]; Fw = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sw[i][j] * Msh.vConv[0][i][j].Vw;
            Pw = Fw / Dw; Msh.tempA[k].aw = Dw * funcScheme(std::abs(Pw)) + std::max(-Fw, 0.0);
            De = gammae * Msh.Se[i][j] / Msh.dX[0][i+1]; Fe = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Se[i][j] * Msh.vConv[0][i][j].Ve;
            Pe = Fe / De; Msh.tempA[k].ae = De * funcScheme(std::abs(Pe)) + std::max(Fe, 0.0);
            Ds = gammas * Msh.Ss[i][j] / Msh.dX[1][j]; Fs = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Ss[i][j] * Msh.vConv[1][i][j].Vs;
            Ps = Fs / Ds; Msh.tempA[k].as = Ds * funcScheme(std::abs(Ps)) + std::max(-Fs, 0.0);
            Dn = gamman * Msh.Sn[i][j] / Msh.dX[1][j+1]; Fn = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sn[i][j] * Msh.vConv[1][i][j].Vn;
            Pn = Fn / Dn; Msh.tempA[k].an = Dn * funcScheme(std::abs(Pn)) + std::max(Fn, 0.0);

            // Coefficients A
            Msh.matA[k].aw = - beta * Msh.tempA[k].aw; Msh.matA[k].ae = - beta * Msh.tempA[k].ae;
            Msh.matA[k].as = - beta * Msh.tempA[k].as; Msh.matA[k].an = - beta * Msh.tempA[k].an; 
            Msh.matA[k].ap = Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp * Msh.Vp[i][j] / dt - Msh.matA[k].ae - Msh.matA[k].aw - Msh.matA[k].an - Msh.matA[k].as;
            
            // Coefficients B 
            Msh.bp[k] = Msh.sPhi[i][j] * Msh.Vp[i][j] + Mat.vMat[Msh.nMat[i][j]].rho * Msh.Vp[i][j] * Msh.vPhi[i][j] / dt - (1 - beta) * (Msh.vPhi[i+1][j] * Msh.tempA[k].ae + Msh.vPhi[i-1][j] * Msh.tempA[k].aw + Msh.vPhi[i][j+1] * Msh.tempA[k].an + Msh.vPhi[i][j-1] * Msh.tempA[k].as - Msh.vPhi[i][j] * (Msh.tempA[k].ae + Msh.tempA[k].aw + Msh.tempA[k].an + Msh.tempA[k].as));

        }
    }
    
}

// CURRENT CHANGES IN DISCRETIZER 
// WRITING EVERYTHING HERE SO I DON'T WASTE TIME TRYING THE SAME FIXES AND TO REMEMBER EVERYTHING IF I NEED TO REVERT ANYTHING
// 1. + (1 - beta) ---> - (1 - beta)





/* void Discretizer::newSetBoundaries(Material& Mat, Mesh& Msh, ExpressionParser& Prs, double t){ */

/*     // Control */
/*     double gamm{}, Dp{}, Fp{}, Pp{}, sumConv{}, sumDiff{}; int k{}; */
/*     std::vector<int> Pos0{}, Pos1{}; Pos0.resize(Msh.N.size()); Pos1.resize(Msh.N.size()); */

/*     // Boundary Conditions */
/*     for (Boundary bC : Msh.boundaryConditions){ */

/*         // Positions (nD) */
/*         for (size_t i = 0; i < Msh.N.size(); i++){ */
/*             Pos0[i] = std::lower_bound(Msh.Nodes[i].begin(), Msh.Nodes[i].end(), bC.x0[i] - epsFind) - Msh.Nodes[i].begin(); */
/*             Pos1[i] = std::lower_bound(Msh.Nodes[i].begin(), Msh.Nodes[i].end(), bC.x1[i] - epsFind) - Msh.Nodes[i].begin(); */
/*         } */

/*         std::cout << "\n"; */
/*         std::cout << "BC - x0: " << bC.x0[0] << " " << bC.x0[1] << "\n"; */
/*         std::cout << "Positions: " << Msh.Nodes[0][Pos0[0]] << " " << Msh.Nodes[1][Pos0[1]] << "\n"; */
/*         std::cout << "BC - x1: " << bC.x1[0] << " " << bC.x1[1] << "\n"; */
/*         std::cout << "Positions: " << Msh.Nodes[0][Pos1[0]] << " " << Msh.Nodes[1][Pos1[1]] << "\n"; */
/*         std::cout << bC.type << "\n"; */

/*         // Boundaries (Non-nD) */
/*         if (bC.type == 0){ */


/*                 /1* else if (bC.iEq == 1){bC.value == Prs.evaluateCoordinates(bC.iExpr, Msh.Nodes[0][Pos0[i]], Msh.Nodes[1][Pos0[1]])} *1/ */

/*             // Update Value */
/*             if (bC.bUpdate && bC.iEq == 0){bC.value = Prs.evaluateTime(bC.iExpr, t);} */

/*             // Dirichlet */
/*             if (Pos0[0] == Pos1[0]){ */

/*                 // xBoundary */
/*                 std::cout << "xBoundary\n"; */
/*                 for (int i = Pos0[1]; i < Pos1[1]; i++){ */
/*                     // Ignore Corners */
/*                     if (i == 0 || i == Msh.N[1]-1){continue;} */

/*                     // Value */
/*                     std::cout << "Update Values\n"; */
/*                     std::cout << "x: " << Pos0[0] << " " << Msh.Nodes[0][Pos0[0]] << "\n"; */
/*                     std::cout << "y: " << i << " " << Msh.Nodes[1][i] << "\n"; */
/*                     if (bC.bUpdate && bC.iEq == 1){bC.value = Prs.evaluateCoordinates(bC.iExpr, Msh.Nodes[0][Pos0[0]], Msh.Nodes[1][i]);} */
/*                     Msh.vPhi[Pos0[0]][i] = bC.value; */
/*                     std::cout << "value: " << bC.value << "\n"; */

/*                     // Coefficients */
/*                     k = Pos0[0] * Msh.N[1] + i; */
/*                     Msh.matA[k].ap = 1; */
/*                     Msh.bp[k] = bC.value; */

/*                     /1* // Control *1/ */
/*                     /1* if (bIgnore){Msh.nIgnore.push_back({Pos0[0], i});} *1/ */
/*                 } */

/*             } else if (Pos0[1] == Pos1[1]){ */
                
/*                 // yBoundary */
/*                 std::cout << "\nyBoundary\n"; */
/*                 for (int i = Pos0[0]; i < Pos1[0]; i++){ */
/*                     // Ignore Corners */
/*                     if (i == 0 || i == Msh.N[0]-1){continue;} */
                    
/*                     // Value */
/*                     std::cout << "Update Values\n"; */
/*                     std::cout << "x: " << i << " " << Msh.Nodes[0][i] << "\n"; */
/*                     std::cout << "y: " << Pos0[1] << " " << Msh.Nodes[1][Pos0[1]] << "\n"; */
/*                     if (bC.bUpdate && bC.iEq == 1){bC.value = Prs.evaluateCoordinates(bC.iExpr, Msh.Nodes[0][i], Msh.Nodes[1][Pos0[1]]);} */
/*                     Msh.vPhi[i][Pos0[1]] = bC.value; */
/*                     std::cout << "value: " << bC.value << "\n"; */

/*                     // Coefficients */
/*                     k = i * Msh.N[1] + Pos0[1]; */
/*                     Msh.matA[k].ap = 1; */
/*                     Msh.bp[k] = bC.value; */

/*                     /1* // Control *1/ */
/*                     /1* if (bIgnore){Msh.nIgnore.push_back({i, Pos0[1]});} *1/ */
/*                 } */

/*             } */
            
/*         } else if (bC.type == 1){ // NEED TO CHECK IF RANGES ARE OKAY AND THEN IF VALUES MAKE SENSE FOR EACH EXPRESSION */
/*             // ALSO CHECK IF */ 

/*             // Neumann */
/*             if (Pos0[0] == Pos1[0]){ */
                
/*                 // xBoundary */
/*                 if (bC.side == 0){ */

/*                     // West Boundary // ae, ap, bp */
/*                     for (int i = Pos0[1]; i < Pos1[1]; i++){ */
/*                         // Ignore Corners */
/*                         if (i == 0 || i == Msh.N[1]-1){continue;} */

/*                         // Thermal Conductivity */
/*                         gamm = Mat.vMat[Msh.nMat[Pos0[0]][i]].gamma; */

/*                         // Value */
/*                         /1* Msh.vPhi[Pos0[0]][i] = bC.value * Msh.nd[0][Pos0[0]] / gamm + Msh.vPhi[Pos0[0]+1][i]; *1/ */

/*                         // Index */
/*                         k = Pos0[0] * Msh.N[1] + i; */

/*                         // Coefficients Conv-Diff */
/*                         Dp = gamm * Msh.Se[Pos0[0]][i] / Msh.nd[0][Pos0[0]]; Fp = Mat.vMat[Msh.nMat[Pos0[0]][i]].rho * Msh.Se[Pos0[0]][i] * Msh.vField.vConv[0][Pos0[0]+1][i]; Pp = Fp / Dp; */

/*                         // Coefficients A */
/*                         Msh.matA[k].ae = - beta * (Dp * funcScheme(std::abs(Pp)) + std::max(Fp, 0.0)); */
/*                         Msh.matA[k].ap = - Msh.matA[k].ae; */
                        
/*                         // Coefficients B */
/*                         Msh.bp[k] = bC.value + (1 - beta) * (Msh.vPhi[Pos0[0]+1][i] * Msh.matA[k].ae / beta + Msh.vPhi[Pos0[0]][i] * Msh.matA[k].ap / beta); */
                        
/*                         // OLD (WRONG) */
/*                         /1* sumConv = Mat.vMat[Msh.nMat[Pos0[0]][i]].rho * (Msh.vField.vConv[0][Pos0[0]+1][i] * Msh.vPhi[Pos0[0]+1][i] * Msh.Se[Pos0[0]][i] + Msh.vField.vConv[0][Pos0[0]+1][i] * Msh.Se[Pos0[0]][i] * Msh.vPhi[Pos0[0]][i]); *1/ */
/*                         /1* sumConv = Mat.vMat[Msh.nMat[Pos0[0]][i]].rho * Msh.vField.vConv[0][Pos0[0]+1][i] * Msh.Se[Pos0[0]][i] * Msh.vPhi[Pos0[0]+1][i]; *1/ */
/*                         /1* sumDiff = gamm * (Msh.Se[Pos0[0]][i] * Msh.vPhi[Pos0[0]+1][i] / Msh.nd[0][Pos0[0]] - Msh.Se[Pos0[0]][i] * Msh.vPhi[Pos0[0]][i]) / Msh.nd[0][Pos0[0]]; *1/ */
/*                         /1* Msh.bp[k] = bC.value * Msh.Se[Pos0[0]][i] + (1 - beta) * (sumConv + sumDiff); // Not entirely sure if this should have q_N * S or just q_N *1/ */

/*   			            /1* // Control (Testing) *1/ */
/* 	                    /1* if (bIgnore){Msh.nIgnore.push_back({Pos0[0], i});} *1/ */

/*                     } */

/*                 } else if (bC.side == 1){ */

/*                     // East Boundary // aw, ap, bp */
/*                     for (int i = Pos0[1]; i < Pos1[1]; i++){ */
/*                         // Ignore Corners */
/*                         if (i == 0 || i == Msh.N[1]-1){continue;} */
                        
/*                         // Thermal Conductivity */
/*                         gamm = Mat.vMat[Msh.nMat[Pos0[0]][i]].gamma; */

/*                         // Value */
/*                         /1* Msh.vPhi[Pos0[0]][i] = bC.value * Msh.nd[0][Pos0[0]-1] / gamm + Msh.vPhi[Pos0[0]-1][i]; *1/ */

/*                         // Index */
/*                         k = Pos0[0] * Msh.N[1] + i; */

/*                         // Coefficients Conv-Diff */
/*                         Dp = gamm * Msh.Sw[Pos0[0]][i] / Msh.nd[0][Pos0[0]-1]; Fp = Mat.vMat[Msh.nMat[Pos0[0]][i]].rho * Msh.Sw[Pos0[0]][i] * Msh.vField.vConv[0][Pos0[0]-1][i]; Pp = Fp / Dp; */

/*                         // Coefficients A */
/*                         Msh.matA[k].aw = - beta * (Dp * funcScheme(std::abs(Pp)) + std::max(-Fp, 0.0)); */
/*                         Msh.matA[k].ap = - Msh.matA[k].aw; */ 

/*                         // Coefficients B */
/*                         Msh.bp[k] = bC.value + (1 - beta) * (Msh.vPhi[Pos0[0]-1][i] * Msh.matA[k].aw / beta + Msh.vPhi[Pos0[0]][i] * Msh.matA[k].ap / beta); */

/*                         // OLD (WRONG) */
/*                         /1* sumConv = Mat.vMat[Msh.nMat[Pos0[0]][i]].rho * Msh.vField.vConv[0][Pos0[0]-1][i] * Msh.Sw[Pos0[0]][i] * Msh.vPhi[Pos0[0]][i]; *1/ */
/*                         /1* sumDiff = gamm * (Msh.Sw[Pos0[0]][i] * Msh.vPhi[Pos0[0]-1][i] / Msh.nd[0][Pos0[0]-1] - Msh.Sw[Pos0[0]][i] * Msh.vPhi[Pos0[0]][i] / Msh.nd[0][Pos0[0]-1]); *1/ */
/*                         /1* Msh.bp[k] = bC.value * Msh.Sw[Pos0[0]][i] + (1 - beta) * (sumConv + sumDiff); *1/ */



/*                     // CHANGE COEFFICIENT CALCULATION - INCLUDE DIRECT CALCULATION FIRST AND THEN INCLUDE BETA / (1 - BETA) */
/*                     // OVERHAUL EVERYTHING */
/*                     // INTERIOR NODES ONLY */
/*                     // AHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH */

/*           			/1* // Control (Testing) *1/ */
/* 	                /1* if (bIgnore){Msh.nIgnore.push_back({Pos0[0], i});} *1/ */

/*                     } */

/*                 } else {std::cerr << "Boundary side not specified correcly.\n";} */

/*             } else if (Pos0[1] == Pos1[1]){ */

/*                 // yBoundary */
/*                 std::cout << "\nyBoundary\n"; */
/*                 if (bC.side == 0){ */

/*                     // South Boundary // an, ap, bp */
/*                     for (int i = Pos0[0]; i < Pos1[0]; i++){ */
/*                         // Ignore Corners */
/*                         if (i == 0 || i == Msh.N[0]-1){continue;} */

/*                         // Thermal Conductivity */
/*                         gamm = Mat.vMat[Msh.nMat[i][Pos0[1]]].gamma; */

/*                         // Value */
/*                         /1* Msh.vPhi[i][Pos0[1]] = bC.value * Msh.nd[1][Pos0[1]] / gamm + Msh.vPhi[i][Pos0[1]+1]; *1/ */
/*                         std::cout << "x: " << Pos0[0] << " " << Msh.Nodes[0][Pos0[0]] << "\n"; */
/*                         std::cout << "y: " << i << " " << Msh.Nodes[1][i] << "\n"; */
/*                         /1* std::cout << "value: " << bC.value << "\n"; *1/ */


/*                         // Index */
/*                         k = i * Msh.N[1] + Pos0[1]; */
                        
/*                         // Coefficients Conv-Diff */
/*                         Dp = gamm * Msh.Sn[i][Pos0[1]] / Msh.nd[1][Pos0[1]]; Fp = Mat.vMat[Msh.nMat[i][Pos0[1]]].rho * Msh.Sn[i][Pos0[1]] * Msh.vField.vConv[1][i][Pos0[1]+1]; Pp = Fp / Dp; */
                        
/*                         // Coefficients A */
/*                         Msh.matA[k].an = - beta * (Dp * funcScheme(std::abs(Pp)) + std::max(Fp, 0.0)); */
/*                         Msh.matA[k].ap = - Msh.matA[k].an; */

/*                         // Coefficients B */
/*                         Msh.bp[k] = bC.value + (1 - beta) * (Msh.vPhi[i][Pos0[1]+1] * Msh.matA[k].an / beta + Msh.vPhi[Pos0[0]][i] * Msh.matA[k].ap / beta); */
/*                         std::cout << "ap: " << Msh.matA[k].ap << "\n"; */
/*                         std::cout << "an: " << Msh.matA[k].an << "\n"; */
/*                         std::cout << "bp: " << Msh.bp[k] << "\n"; */
/*                         std::cout << "Others : " << Msh.matA[k].as << " " << Msh.matA[k].ae << " " << Msh.matA[k].aw << "\n"; */


/*                         // OLD (WRONG) */
/*                         /1* sumConv = Mat.vMat[Msh.nMat[i][Pos0[1]]].rho * Msh.vField.vConv[1][i][Pos0[1]+1] * Msh.Sn[i][Pos0[1]] * Msh.vPhi[i][Pos0[1]+1]; *1/ */
/*                         /1* sumDiff = gamm * (Msh.Sn[i][Pos0[1]] * Msh.vPhi[i][Pos0[1]+1] / Msh.nd[1][Pos0[1]] - Msh.Sn[i][Pos0[1]] * Msh.vPhi[i][Pos0[1]] / Msh.nd[1][Pos0[1]]); *1/ */
/*                         /1* Msh.bp[k] = bC.value * Msh.Sn[i][Pos0[1]] + (1 - beta) * (sumConv + sumDiff); // VALIDATE IF AREAS HERE HAVE A VALUE OR IF THEY'RE JUST 0 - IF 0 CHANGE *1/ */

/*   			        /1* // Control (Testing) *1/ */
/* 	                /1* if (bIgnore){Msh.nIgnore.push_back({i, Pos0[1]});} *1/ */

/*                     } */

/*                 } else if (bC.side == 1){ */

/*                     // North Boundary // as, ap, bp */
/*                     for (int i = Pos0[0]; i < Pos1[0]; i++){ */
/*                         // Ignore Corners */
/*                         if (i == 0 || i == Msh.N[0]-1){continue;} */
                        
/*                         // Thermal Conductivity */
/*                         gamm = Mat.vMat[Msh.nMat[i][Pos0[1]]].gamma; */

/*                         // Value */
/*                         Msh.vPhi[i][Pos0[1]] = bC.value * Msh.nd[1][Pos0[1]-1] / gamm + Msh.vPhi[i][Pos0[1]-1]; */

/*                         // Index */
/*                         k = i * Msh.N[1] + Pos0[1]; */

/*                         // Coefficients Conv-Diff */
/*                         Dp = gamm * Msh.Ss[i][Pos0[1]] / Msh.nd[1][Pos0[1]-1]; Fp = Mat.vMat[Msh.nMat[i][Pos0[1]]].rho * Msh.Ss[i][Pos0[1]] * Msh.vField.vConv[1][i][Pos0[1]-1]; Pp = Fp / Dp; */

/*                         // Coefficients A */
/*                         Msh.matA[k].as = - beta * (Dp * funcScheme(std::abs(Pp)) + std::max(-Fp, 0.0)); */
/*                         Msh.matA[k].ap = - Msh.matA[k].as; */

/*                         // Coefficients B */
/*                         Msh.bp[k] = bC.value + (1 - beta) * (Msh.vPhi[i][Pos0[1]-1] * Msh.matA[k].as / beta + Msh.vPhi[Pos0[0]][i] * Msh.matA[k].ap / beta); */

/*                         // OLD (WRONG) */
/*                         /1* sumConv = Mat.vMat[Msh.nMat[i][Pos0[1]]].rho * Msh.vField.vConv[1][i][Pos0[1]-1] * Msh.Ss[i][Pos0[1]] * Msh.vPhi[i][Pos0[1]]; *1/ */
/*                         /1* sumDiff = gamm * (Msh.Ss[i][Pos0[1]] * Msh.vPhi[i][Pos0[1]-1] / Msh.nd[1][Pos0[1]-1] - Msh.Ss[i][Pos0[1]] * Msh.vPhi[i][Pos0[1]] / Msh.nd[1][Pos0[1]-1]); *1/ */
/*                         /1* Msh.bp[k] = bC.value * Msh.Ss[i][Pos0[1]] + (1 - beta) * (sumConv + sumDiff); *1/ */

/*   	        		/1* // Control (Testing) *1/ */
/* 	                /1* if (bIgnore){Msh.nIgnore.push_back({i, Pos0[1]});} *1/ */

/*                     } */

/*                 } else {std::cerr << "Boundary side not specified correcly.\n";} */

/*             } */

/*         } else if (bC.type == 2){ */

/*         } else {std::cerr << "Boundary side not specified correctly.\n";} */

/*     } */


/*     // Corners */
/*     Msh.vPhi[0][0] = 0.5 * (Msh.vPhi[1][0] + Msh.vPhi[0][1]); */
/*     Msh.vPhi[0][Msh.N[1]-1] = 0.5 * (Msh.vPhi[1][Msh.N[1]-1] + Msh.vPhi[0][Msh.N[1]-2]); */
/*     Msh.vPhi[Msh.N[0]-1][0] = 0.5 * (Msh.vPhi[Msh.N[0]-2][0] + Msh.vPhi[Msh.N[0]-1][1]); */
/*     Msh.vPhi[Msh.N[0]-1][Msh.N[1]-1] = 0.5 * (Msh.vPhi[Msh.N[0]-2][Msh.N[1]-1] + Msh.vPhi[Msh.N[0]-1][Msh.N[1]-2]); */

/*     // Corners */
/*     int l{}, m{}, n{}; */

/*     l = 0; m = 0; n = l * Msh.N[1] + m; // SW */
/*     Msh.vPhi[l][m] = 0.5 * (Msh.vPhi[l+1][m] + Msh.vPhi[l][m+1]); */
/*     Msh.matA[n].ap = 2; Msh.matA[n].ae = -1; Msh.matA[n].an = -1; */

/*     l = 0; m = Msh.N[1]-1; n = l * Msh.N[1] + m; // NW */
/*     Msh.vPhi[l][m] = 0.5 * (Msh.vPhi[l+1][m] + Msh.vPhi[l][m-1]); */
/*     Msh.matA[n].ap = 2; Msh.matA[n].ae = -1; Msh.matA[n].as = -1; */

/*     l = Msh.N[0]-1; m = 0; n = l * Msh.N[1] + m; // SE */
/*     Msh.vPhi[l][m] = 0.5 * (Msh.vPhi[l-1][m] + Msh.vPhi[l][m+1]); */
/*     Msh.matA[n].ap = 2; Msh.matA[n].aw = -1; Msh.matA[n].an = -1; */

/*     l = Msh.N[0]-1; m = Msh.N[1]-1; n = l * Msh.N[1] + m; // NE */
/*     Msh.vPhi[l][m] = 0.5 * (Msh.vPhi[l-1][m] + Msh.vPhi[l][m-1]); */
/*     Msh.matA[n].ap = 2; Msh.matA[n].aw = -1; Msh.matA[n].as = -1; */

/*     /1* // Control *1/ */
/*     /1* if (bIgnore){Msh.nIgnore.push_back({0, 0}); Msh.nIgnore.push_back({0, Msh.N[1]-1}); Msh.nIgnore.push_back({Msh.N[0]-1, 0}); Msh.nIgnore.push_back({Msh.N[0]-1, Msh.N[1]-1});} *1/ */
/*     /1* bIgnore = false; *1/ */


/* } */






/* void Discretizer::newSetRHS(Material& Mat, Mesh& Msh){ */

/*     // Control */
/*     double gammaw{}, gammae{}, gammas{}, gamman{}, sumConv{}, sumDiff{}; int k; */
/*     double Fw{}, Fe{}, Fs{}, Fn{}, Dw{}, De{}, Ds{}, Dn{}, Pw{}, Pe{}, Ps{}, Pn{}; */

/*     // Interior Nodes (Non-nD) */
/*     for (size_t i = 1; i < Msh.N[0]-1; i++){ */
/*         for (size_t j = 1; j < Msh.N[1]-1; j++){ */
            
/*             // Harmonic Mean */
/*             gammaw = calcHarmonicMean(Msh.nd[0][i-1], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i-1][j]].gamma}, {Msh.ndelta[0][i], Msh.ndelta[0][i-1]}); */
/*             gammae = calcHarmonicMean(Msh.nd[0][i], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i+1][j]].gamma}, {Msh.ndelta[0][i], Msh.ndelta[0][i+1]}); */
/*             gammas = calcHarmonicMean(Msh.nd[1][j-1], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i][j-1]].gamma}, {Msh.ndelta[1][j], Msh.ndelta[1][j-1]}); */
/*             gamman = calcHarmonicMean(Msh.nd[1][j], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i][j+1]].gamma}, {Msh.ndelta[1][j], Msh.ndelta[1][j+1]}); */

/*             // Index */
/*             k = i * Msh.N[1] + j; */

/*             // Coefficients Conv-Diff */
/*             De = - gammae * Msh.Se[i][j] / Msh.nd[0][i]; Fe = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Se[i][j] * Msh.vField.vConv[0][i+1][j]; Pe = Fe / De; */
/*             Dw = - gammaw * Msh.Sw[i][j] / Msh.nd[0][i-1]; Fw = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sw[i][j] * Msh.vField.vConv[0][i-1][j]; Pw = Fw / Dw; */
/*             Dn = - gamman * Msh.Sn[i][j] / Msh.nd[1][j]; Fn = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sn[i][j] * Msh.vField.vConv[1][i][j+1]; Pn = Fn / Dn; */
/*             Ds = - gammas * Msh.Ss[i][j] / Msh.nd[1][j-1]; Fs = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Ss[i][j] * Msh.vField.vConv[1][i][j-1]; Ps = Fs / Ds; */

/*             // Coefficients B */ 
/*             Msh.bp[k] = Msh.sPhi[i][j] * Msh.Vp[i][j] + Mat.vMat[Msh.nMat[i][j]].rho * Msh.Vp[i][j] * Msh.vPhi[i][j] / dt + (1 - beta) * (Msh.vPhi[i+1][j] * Msh.matA[k].ae / beta + Msh.vPhi[i-1][j] * Msh.matA[k].aw / beta + Msh.vPhi[i][j+1] * Msh.matA[k].an / beta + Msh.vPhi[i][j-1] * Msh.matA[k].as / beta + Msh.vPhi[i][j] * Msh.matA[k].ap / beta); */

/*             /1* sumConv = Mat.vMat[Msh.nMat[i][j]].rho * (Msh.vField.vConv[0][i+1][j] * Msh.vPhi[i+1][j] * Msh.Se[i][j] - Msh.vField.vConv[0][i-1][j] * Msh.vPhi[i-1][j] * Msh.Sw[i][j] + Msh.vField.vConv[1][i][j+1] * Msh.vPhi[i][j+1] * Msh.Sn[i][j] - Msh.vField.vConv[1][i][j-1] * Msh.vPhi[i][j-1] * Msh.Ss[i][j] + (Msh.vField.vConv[0][i+1][j] * Msh.Se[i][j] - Msh.vField.vConv[0][i-1][j] * Msh.Sw[i][j] + Msh.vField.vConv[1][i][j+1] * Msh.Sn[i][j] - Msh.vField.vConv[1][i][j-1] * Msh.Ss[i][j]) * Msh.vPhi[i][j]); *1/ */
/*             /1* sumDiff = (gammaw * Msh.Sw[i][j] * Msh.vPhi[i-1][j] / Msh.nd[0][i-1] + gammae * Msh.Se[i][j] * Msh.vPhi[i+1][j] / Msh.nd[0][i] + gammas * Msh.Ss[i][j] * Msh.vPhi[i][j-1] / Msh.nd[1][j-1] + gamman * Msh.Sn[i][j] * Msh.vPhi[i][j+1] / Msh.nd[1][j] - (gammaw * Msh.Sw[i][j] / Msh.nd[0][i-1] + gammae * Msh.Se[i][j] / Msh.nd[0][i] + gammas * Msh.Ss[i][j] / Msh.nd[1][j-1] + gamman * Msh.Sn[i][j] / Msh.nd[1][j]) * Msh.vPhi[i][j]); *1/ */
/*             /1* Msh.bp[k] = Msh.sPhi[i][j] * Msh.Vp[i][j] + Mat.vMat[Msh.nMat[i][j]].rho * Msh.Vp[i][j] * Msh.vPhi[i][j] / dt + (1 - beta) * (sumConv + sumDiff); *1/ */
/*             /1* Msh.bp[k] = Msh.sPhi[i][j] * Msh.Vp[i][j] + Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp * Msh.Vp[i][j] * Msh.vPhi[i][j] / dt + (1 - beta) * (sumConv + sumDiff); *1/ */

/*         } */
/*     } */

/* } */
