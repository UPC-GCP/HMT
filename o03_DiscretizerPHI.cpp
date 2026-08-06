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


void Discretizer::newSetBoundaries(Material& Mat, Mesh& Msh, ExpressionParser& Prs, double t){

    // Control
    Msh.tempA.assign(Msh.totNodes, {}); Msh.temp2A.assign(Msh.totNodes, {}); Msh.tempB.assign(Msh.totNodes, 0);
    int i{}, j{}, k{}; double gammaw{}, gammae{}, gammas{}, gamman{}, a0{}, tempStore{};
    double Fw{}, Fe{}, Fs{}, Fn{}, Dw{}, De{}, Ds{}, Dn{}, Pw{}, Pe{}, Ps{}, Pn{};
    std::vector<int> Pos0{}, Pos1{}; Pos0.resize(Msh.N.size()); Pos1.resize(Msh.N.size());

    ////////// Boundary Node Coefficients //////////
    for (Boundary& bC : Msh.boundaryConditions){

        // Positions
        for (size_t i = 0; i < Msh.N.size(); i++){
            Pos0[i] = std::lower_bound(Msh.Nodes[i].begin(), Msh.Nodes[i].end(), bC.x0[i] - epsFind) - Msh.Nodes[i].begin();
            Pos1[i] = std::lower_bound(Msh.Nodes[i].begin(), Msh.Nodes[i].end(), bC.x1[i] - epsFind) - Msh.Nodes[i].begin();
        }

        // Boundaries (Non-nD)
        if (bC.type == 0){
            
            // Update Value 
            tempStore = bC.value;
            if (bC.bUpdate && bC.iEq == 0){bC.value = Prs.evaluateTime(bC.iExpr, t); bC.Phi.assign(bC.Phi.size(), bC.value); bC.oPhi.assign(bC.oPhi.size(), tempStore);}

            // Dirichlet
            if (Pos0[0] == Pos1[0]){

                // xBoundary
                if (bC.side == 0){

                    // West Boundary
                    i = Pos0[0];
                    for (size_t j = Pos0[1]; j < Pos1[1]; j++){
                        // Control
                        k = i * Msh.N[1] + j;
                        if (bC.bUpdate && bC.iEq == 1){bC.value = Prs.evaluateCoordinates(bC.iExpr, Msh.Faces[0][i], Msh.Nodes[1][j]); tempStore = bC.value; bC.Phi[j] = bC.value; bC.oPhi[j] = bC.value;}
                        
                        // Gamma
                        gammaw = Mat.vMat[Msh.nMat[i][j]].gamma;

                        // Coefficients Conv-Diff
                        Dw = gammaw * Msh.Sw[i][j] / Msh.dX[0][i]; Fw = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sw[i][j] * Msh.vConv[0][i][j].Vw;
                        Pw = Fw / Dw; Msh.tempA[k].aw = Dw * funcScheme(Pw) + std::max(Fw, 0.0); Msh.temp2A[k].aw = Dw * funcScheme(Pw) + std::max(Fw, 0.0);

                        // Coefficients B
                        Msh.tempB[k] += beta * Msh.tempA[k].aw * bC.Phi[j] + (1 - beta) * Msh.tempA[k].aw * bC.oPhi[j];
                    }

                } else if (bC.side == 1){
                    
                    // East Boundary
                    i = Pos0[0] - 1;
                    for (size_t j = Pos0[1]; j < Pos1[1]; j++){
                        // Control
                        k = i * Msh.N[1] + j;
                        if (bC.bUpdate && bC.iEq == 1){bC.value = Prs.evaluateCoordinates(bC.iExpr, Msh.Faces[0][i+1], Msh.Nodes[1][j]); tempStore = bC.value; bC.Phi[j] = bC.value; bC.oPhi[j] = bC.value;}

                        // Gamma
                        gammae = Mat.vMat[Msh.nMat[i][j]].gamma;

                        // Coefficients Conv-Diff
                        De = gammae * Msh.Se[i][j] / Msh.dX[0][i+1]; Fe = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Se[i][j] * Msh.vConv[0][i][j].Ve;
                        Pe = Fe / De; Msh.tempA[k].ae = De * funcScheme(Pe) + std::max(-Fe, 0.0); Msh.temp2A[k].ae = De * funcScheme(Pe) + std::max(-Fe, 0.0);

                        // Coefficients B
                        Msh.tempB[k] += beta * Msh.tempA[k].ae * bC.Phi[j] + (1 - beta) * Msh.tempA[k].ae * bC.oPhi[j];
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
                        if (bC.bUpdate && bC.iEq == 1){bC.value = Prs.evaluateCoordinates(bC.iExpr, Msh.Nodes[0][i], Msh.Faces[1][j]); tempStore = bC.value; bC.Phi[i] = bC.value; bC.oPhi[i] = bC.value;}

                        // Gamma
                        gammas = Mat.vMat[Msh.nMat[i][j]].gamma;

                        // Coefficients Conv-Diff
                        Ds = gammas * Msh.Ss[i][j] / Msh.dX[1][j]; Fs = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Ss[i][j] * Msh.vConv[1][i][j].Vs;
                        Ps = Fs / Ds; Msh.tempA[k].as = Ds * funcScheme(Ps) + std::max(Fs, 0.0); Msh.temp2A[k].as = Ds * funcScheme(Ps) + std::max(Fs, 0.0);

                        // Coefficients B
                        Msh.tempB[k] += beta * Msh.tempA[k].as * bC.Phi[i] + (1 - beta) * Msh.tempA[k].as * bC.oPhi[i];
                    }

                } else if (bC.side == 1){

                    // North Boundary
                    j = Pos0[1] - 1;
                    for (size_t i = Pos0[0]; i < Pos1[0]; i++){
                        // North Boundary
                        k = i * Msh.N[1] + j;
                        if (bC.bUpdate && bC.iEq == 1){bC.value = Prs.evaluateCoordinates(bC.iExpr, Msh.Nodes[0][i], Msh.Faces[1][j+1]); tempStore = bC.value; bC.Phi[i] = bC.value; bC.oPhi[i] = bC.value;}

                        // Gamma
                        gamman = Mat.vMat[Msh.nMat[i][j]].gamma;

                        // Coefficients Conv-Diff
                        Dn = gamman * Msh.Sn[i][j] / Msh.dX[1][j+1]; Fn = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sn[i][j] * Msh.vConv[1][i][j].Vn;
                        Pn = Fn / Dn; Msh.tempA[k].an = Dn * funcScheme(Pn) + std::max(-Fn, 0.0); Msh.temp2A[k].an = Dn * funcScheme(Pn) + std::max(-Fn, 0.0);

                        // Coefficients B
                        Msh.tempB[k] += beta * Msh.tempA[k].an * bC.Phi[i] + (1 - beta) * Msh.tempA[k].an * bC.oPhi[i];
                    }

                } else {std::cerr << "Boundary range not specified correctly.\n";}

            }
            
        } else if (bC.type == 1){

            // Neumann
            if (Pos0[0] == Pos1[0]){

                // xBoundary
                if (bC.side == 0){

                    // West Boundary
                    i = Pos0[0]; 
                    for (size_t j = Pos0[1]; j < Pos1[1]; j++){
                        // Control
                        k = i * Msh.N[1] + j; gammaw = Mat.vMat[Msh.nMat[i][j]].gamma;

                        // Coefficients A
                        Dw = gammaw * Msh.Sw[i][j] / Msh.dX[0][i]; Fw = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sw[i][j] * Msh.vConv[0][i][j].Vw;
                        Pw = Fw / Dw; Msh.tempA[k].aw = std::max(Fw, 0.0); Msh.temp2A[k].aw = std::max(Fw, 0.0);

                        // Coefficients B
                        Msh.bcPhi[i][j] = (bC.value + Msh.vPhi[i][j] * gammaw / Msh.dX[0][i]) / (gammaw / Msh.dX[0][i]);
                        bC.oPhi[j] = bC.Phi[j]; bC.Phi[j] = Msh.bcPhi[i][j];
                        Msh.tempB[k] += bC.value * Msh.Sw[i][j] + beta * Msh.tempA[k].aw * bC.Phi[j] + (1 - beta) * Msh.tempA[k].aw * bC.oPhi[j];
                    }

                } else if (bC.side == 1){
                    
                    // East Boundary
                    i = Pos0[0]-1; 
                    for (size_t j = Pos0[1]; j < Pos1[1]; j++){
                        // Control
                        k = i * Msh.N[1] + j; gammae = Mat.vMat[Msh.nMat[i][j]].gamma;

                        // Coefficients A
                        Msh.tempA[k].ae = 0;
                        De = gammae * Msh.Se[i][j] / Msh.dX[0][i+1]; Fe = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Se[i][j] * Msh.vConv[0][i][j].Ve;
                        Pe = Fe / De; Msh.tempA[k].ae = std::max(-Fe, 0.0); Msh.temp2A[k].ae = std::max(-Fe, 0.0);

                        // Coefficients B
                        Msh.bcPhi[i][j] = (bC.value + Msh.vPhi[i][j] * gammae / Msh.dX[0][i+1]) / (gammae / Msh.dX[0][i+1]);
                        bC.oPhi[j] = bC.Phi[j]; bC.Phi[j] = Msh.bcPhi[i][j];
                        Msh.tempB[k] += bC.value * Msh.Se[i][j] + beta * Msh.tempA[k].ae * bC.Phi[j] + (1 - beta) * Msh.tempA[k].ae * bC.oPhi[j];
                    }

                } else {std::cerr << "Boundary side not specified correctly.\n";}
                
            } else if (Pos0[1] == Pos1[1]){

                // yBoundary
                if (bC.side == 0){
                    // South Boundary
                    j = Pos0[1];

                    for (size_t i = Pos0[0]; i < Pos1[0]; i++){
                        // Control
                        k = i * Msh.N[1] + j; gammas = Mat.vMat[Msh.nMat[i][j]].gamma;

                        // Coefficients A
                        Msh.tempA[k].as = 0;
                        Ds = gammas * Msh.Ss[i][j] / Msh.dX[1][j]; Fs = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Ss[i][j] * Msh.vConv[1][i][j].Vs;
                        Ps = Fs / Ds; Msh.tempA[k].as = std::max(Fs, 0.0); Msh.temp2A[k].as = std::max(Fs, 0.0);

                        // Coefficients B
                        Msh.bcPhi[i][j] = (bC.value + Msh.vPhi[i][j] * gammas / Msh.dX[1][j]) / (gammas / Msh.dX[1][j]);
                        bC.oPhi[i] = bC.Phi[i]; bC.Phi[i] = Msh.bcPhi[i][j];
                        Msh.tempB[k] += bC.value * Msh.Ss[i][j] + beta * Msh.tempA[k].as * bC.Phi[i] + (1 - beta) * Msh.tempA[k].as * bC.oPhi[i];
                    }
                    
                } else if (bC.side == 1){

                    // North Boundary
                    j = Pos0[1]-1;
                    for (size_t i = Pos0[0]; i < Pos1[0]; i++){
                        // Control
                        k = i * Msh.N[1] + j; gamman = Mat.vMat[Msh.nMat[i][j]].gamma;

                        // Coefficients A
                        Msh.tempA[k].an = 0;
                        Dn = gamman * Msh.Sn[i][j] / Msh.dX[1][j+1]; Fn = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sn[i][j] * Msh.vConv[1][i][j].Vn;
                        Pn = Fn / Dn; Msh.tempA[k].an = std::max(-Fn, 0.0); Msh.temp2A[k].an = std::max(-Fn, 0.0);

                        // Coefficients B
                        Msh.bcPhi[i][j] = (bC.value + Msh.vPhi[i][j] * gamman / Msh.dX[1][j+1]) / (gamman / Msh.dX[1][j+1]);
                        bC.oPhi[i] = bC.Phi[i]; bC.Phi[i] = Msh.bcPhi[i][j];
                        Msh.tempB[k] += bC.value * Msh.Sn[i][j] + beta * Msh.tempA[k].an * bC.Phi[i] + (1 - beta) * Msh.tempA[k].an * bC.oPhi[i];
                    }

                } else {std::cerr << "Boundary side not specified correctly.\n";}
                
            } else {std::cerr << "Boundary range not specified correctly.\n";}


        } else if (bC.type == 2){

            // Robin (Hybrid)
            if (Pos0[0] == Pos1[0]){

                // xBoundary
                if (bC.side == 0){

                    // West Boundary
                    i = Pos0[0];
                    for (size_t j = Pos0[1]; j < Pos1[1]; j++){
                        // Control
                        k = i * Msh.N[1] + j; gammaw = Mat.vMat[Msh.nMat[i][j]].gamma;

                        // Coefficients A
                        Msh.tempA[k].aw = Msh.Sw[i][j] * (bC.alpha * gammaw / Msh.dX[0][i]) / (bC.alpha + gammaw / Msh.dX[0][i]);

                        // Coefficients B
                        Msh.tempB[k] += Msh.tempA[k].aw * bC.value + (1 - beta) * (Msh.tempA[k].aw * (bC.value - Msh.vPhi[i][j]));
                    }

                } else if (bC.side == 1){

                    // East Boundary
                    i = Pos0[0] - 1;
                    for (size_t j = Pos0[1]; j < Pos1[1]; j++){
                        // Control
                        k = i * Msh.N[1] + j; gammae = Mat.vMat[Msh.nMat[i][j]].gamma;

                        // Coefficients A
                        Msh.tempA[k].ae = Msh.Se[i][j] * (bC.alpha * gammae / Msh.dX[0][i+1]) / (bC.alpha + gammae / Msh.dX[0][i+1]);

                        // Coefficients B
                        Msh.tempB[k] += Msh.tempA[k].ae * bC.value + (1 - beta) * (Msh.tempA[k].ae * (bC.value - Msh.vPhi[i][j]));
                    }

                } else {std::cerr << "Boundary side not specified correctly.\n";}

            } else if (Pos0[1] == Pos1[1]){

                // yBoundary
                if (bC.side == 0){

                    // South Boundary
                    j = Pos0[1];
                    for (size_t i = Pos0[0]; i < Pos1[0]; i++){
                        // Control
                        k = i * Msh.N[1] + j; gammas = Mat.vMat[Msh.nMat[i][j]].gamma;

                        // Coefficients A
                        Msh.tempA[k].as = Msh.Ss[i][j] * (bC.alpha * gammas / Msh.dX[1][j]) / (bC.alpha + gammas / Msh.dX[1][j]);

                        // Coefficients B
                        Msh.tempB[k] += Msh.tempA[k].as * bC.value + (1 - beta) * (Msh.tempA[k].as * (bC.value - Msh.vPhi[i][j]));
                    }

                } else if (bC.side == 1){

                    // North Boundary
                    j = Pos0[1] - 1;
                    for (size_t i = Pos0[0]; i < Pos1[0]; i++){
                        // Control
                        k = i * Msh.N[1] + j; gamman = Mat.vMat[Msh.nMat[i][j]].gamma;

                        // Coefficients A
                        Msh.tempA[k].an = Msh.Sn[i][j] * (bC.alpha * gamman / Msh.dX[1][j+1]) / (bC.alpha + gamman / Msh.dX[1][j+1]);

                        // Coefficients B
                        Msh.tempB[k] += Msh.tempA[k].an * bC.value + (1 - beta) * (Msh.tempA[k].an * (bC.value - Msh.vPhi[i][j]));
                    }

                } else {std::cerr << "Boundary side not specified correctly.\n";}

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
            Ps = Fs / Ds; Msh.tempA[k].as = Ds * funcScheme(Ps) + std::max(Fs, 0.0);
            Msh.tempB[k] += (1 - beta) * Msh.vPhi[i][j-1] * Msh.tempA[k].as;
        } 

        // NW Corner (~an)
        if (j != Msh.N[1]-1){
            gamman = calcHarmonicMean(Msh.dX[1][j+1], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i][j+1]].gamma}, {Msh.DeltaX[1][j], Msh.DeltaX[1][j+1]});
            Dn = gamman * Msh.Sn[i][j] / Msh.dX[1][j+1]; Fn = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sn[i][j] * Msh.vConv[1][i][j].Vn;
            Pn = Fn / Dn; Msh.tempA[k].an = Dn * funcScheme(Pn) + std::max(-Fn, 0.0);
            Msh.tempB[k] += (1 - beta) * Msh.vPhi[i][j+1] * Msh.tempA[k].an;
        }

        // W Edge (ae)
        gammae = calcHarmonicMean(Msh.dX[0][i+1], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i+1][j]].gamma}, {Msh.DeltaX[0][i], Msh.DeltaX[0][i+1]});
        De = gammae * Msh.Se[i][j] / Msh.dX[0][i+1]; Fe = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Se[i][j] * Msh.vConv[0][i][j].Ve;
        Pe = Fe / De; Msh.tempA[k].ae = De * funcScheme(Pe) + std::max(-Fe, 0.0);
       
        // Coefficients A
        a0 = Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp * Msh.Vp[i][j] / dt;
        Msh.matA[k].aw = - beta * Msh.tempA[k].aw; Msh.matA[k].ae = - beta * Msh.tempA[k].ae;
        Msh.matA[k].as = - beta * Msh.tempA[k].as; Msh.matA[k].an = - beta * Msh.tempA[k].an; 
        Msh.matA[k].ap = a0 - Msh.matA[k].ae - Msh.matA[k].aw - Msh.matA[k].an - Msh.matA[k].as;
            
        // Coefficients B
        Msh.bp[k] = Msh.sPhi[i][j] * Msh.Vp[i][j] + a0 * Msh.vPhi[i][j] + (1 - beta) * (Msh.vPhi[i+1][j] * Msh.tempA[k].ae - Msh.vPhi[i][j] * (Msh.tempA[k].aw + Msh.tempA[k].ae + Msh.tempA[k].as + Msh.tempA[k].an)) + Msh.tempB[k];

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
            Ps = Fs / Ds; Msh.tempA[k].as = Ds * funcScheme(Ps) + std::max(Fs, 0.0);
            Msh.tempB[k] += (1 - beta) * Msh.vPhi[i][j-1] * Msh.tempA[k].as;
        }

        // NE Corner (~an)
        if (j != Msh.N[1]-1){
            gamman = calcHarmonicMean(Msh.dX[1][j+1], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i][j+1]].gamma}, {Msh.DeltaX[1][j], Msh.DeltaX[1][j+1]});
            Dn = gamman * Msh.Sn[i][j] / Msh.dX[1][j+1]; Fn = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sn[i][j] * Msh.vConv[1][i][j].Vn;
            Pn = Fn / Dn; Msh.tempA[k].an = Dn * funcScheme(Pn) + std::max(-Fn, 0.0);
            Msh.tempB[k] += (1 - beta) * Msh.vPhi[i][j+1] * Msh.tempA[k].an;
        }

        // East Edge (aw)
        gammaw = calcHarmonicMean(Msh.dX[0][i], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i-1][j]].gamma}, {Msh.DeltaX[0][i], Msh.DeltaX[0][i-1]});
        Dw = gammaw * Msh.Sw[i][j] / Msh.dX[0][i]; Fw = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sw[i][j] * Msh.vConv[0][i][j].Vw;
        Pw = Fw / Dw; Msh.tempA[k].aw = Dw * funcScheme(Pw) + std::max(Fw, 0.0);

        // Coefficients A
        a0 = Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp * Msh.Vp[i][j] / dt;
        Msh.matA[k].aw = - beta * Msh.tempA[k].aw; Msh.matA[k].ae = - beta * Msh.tempA[k].ae;
        Msh.matA[k].as = - beta * Msh.tempA[k].as; Msh.matA[k].an = - beta * Msh.tempA[k].an; 
        Msh.matA[k].ap = a0 - Msh.matA[k].ae - Msh.matA[k].aw - Msh.matA[k].an - Msh.matA[k].as;
            
        // Coefficients B
        Msh.bp[k] = Msh.sPhi[i][j] * Msh.Vp[i][j] + a0 * Msh.vPhi[i][j] + (1 - beta) * (Msh.vPhi[i-1][j] * Msh.tempA[k].aw - Msh.vPhi[i][j] * (Msh.tempA[k].aw + Msh.tempA[k].ae + Msh.tempA[k].as + Msh.tempA[k].an)) + Msh.tempB[k];
  
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
        Pw = Fw / Dw; Msh.tempA[k].aw = Dw * funcScheme(Pw) + std::max(Fw, 0.0);
        De = gammae * Msh.Se[i][j] / Msh.dX[0][i+1]; Fe = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Se[i][j] * Msh.vConv[0][i][j].Ve;
        Pe = Fe / De; Msh.tempA[k].ae = De * funcScheme(Pe) + std::max(-Fe, 0.0);
        Dn = gamman * Msh.Sn[i][j] / Msh.dX[1][j+1]; Fn = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sn[i][j] * Msh.vConv[1][i][j].Vn;
        Pn = Fn / Dn; Msh.tempA[k].an = Dn * funcScheme(Pn) + std::max(-Fn, 0.0);
        
        // Coefficients A
        a0 = Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp * Msh.Vp[i][j] / dt;
        Msh.matA[k].aw = - beta * Msh.tempA[k].aw; Msh.matA[k].ae = - beta * Msh.tempA[k].ae;
        Msh.matA[k].as = - beta * Msh.tempA[k].as; Msh.matA[k].an = - beta * Msh.tempA[k].an; 
        Msh.matA[k].ap = a0 - Msh.matA[k].ae - Msh.matA[k].aw - Msh.matA[k].an - Msh.matA[k].as;
            
        // Coefficients B
        Msh.bp[k] = Msh.sPhi[i][j] * Msh.Vp[i][j] + a0 * Msh.vPhi[i][j] + (1 - beta) * (Msh.vPhi[i-1][j] * Msh.tempA[k].aw + Msh.vPhi[i+1][j] * Msh.tempA[k].ae + Msh.vPhi[i][j+1] * Msh.tempA[k].an - Msh.vPhi[i][j] * (Msh.tempA[k].aw + Msh.tempA[k].ae + Msh.tempA[k].as + Msh.tempA[k].an)) + Msh.tempB[k];

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
        Pw = Fw / Dw; Msh.tempA[k].aw = Dw * funcScheme(Pw) + std::max(Fw, 0.0);
        De = gammae * Msh.Se[i][j] / Msh.dX[0][i+1]; Fe = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Se[i][j] * Msh.vConv[0][i][j].Ve;
        Pe = Fe / De; Msh.tempA[k].ae = De * funcScheme(Pe) + std::max(-Fe, 0.0);
        Ds = gammas * Msh.Ss[i][j] / Msh.dX[1][j]; Fs = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Ss[i][j] * Msh.vConv[1][i][j].Vs;
        Ps = Fs / Ds; Msh.tempA[k].as = Ds * funcScheme(Ps) + std::max(Fs, 0.0);
        
        // Coefficients A
        a0 = Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp * Msh.Vp[i][j] / dt;
        Msh.matA[k].aw = - beta * Msh.tempA[k].aw; Msh.matA[k].ae = - beta * Msh.tempA[k].ae;
        Msh.matA[k].as = - beta * Msh.tempA[k].as; Msh.matA[k].an = - beta * Msh.tempA[k].an; 
        Msh.matA[k].ap = a0 - Msh.matA[k].ae - Msh.matA[k].aw - Msh.matA[k].an - Msh.matA[k].as;
            
        // Coefficients B
        Msh.bp[k] = Msh.sPhi[i][j] * Msh.Vp[i][j] + a0 * Msh.vPhi[i][j] + (1 - beta) * (Msh.vPhi[i-1][j] * Msh.tempA[k].aw + Msh.vPhi[i+1][j] * Msh.tempA[k].ae + Msh.vPhi[i][j-1] * Msh.tempA[k].as - Msh.vPhi[i][j] * (Msh.tempA[k].aw + Msh.tempA[k].ae + Msh.tempA[k].as + Msh.tempA[k].an)) + Msh.tempB[k];
    
    }

}


void Discretizer::newSetCoefficients(Material& Mat, Mesh& Msh){

    // Control
    int k{}; double gammaw, gammae, gammas{}, gamman{}, a0{};
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
            Pw = Fw / Dw; Msh.tempA[k].aw = Dw * funcScheme(Pw) + std::max(Fw, 0.0);
            De = gammae * Msh.Se[i][j] / Msh.dX[0][i+1]; Fe = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Se[i][j] * Msh.vConv[0][i][j].Ve;
            Pe = Fe / De; Msh.tempA[k].ae = De * funcScheme(Pe) + std::max(-Fe, 0.0);
            Ds = gammas * Msh.Ss[i][j] / Msh.dX[1][j]; Fs = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Ss[i][j] * Msh.vConv[1][i][j].Vs;
            Ps = Fs / Ds; Msh.tempA[k].as = Ds * funcScheme(Ps) + std::max(Fs, 0.0);
            Dn = gamman * Msh.Sn[i][j] / Msh.dX[1][j+1]; Fn = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sn[i][j] * Msh.vConv[1][i][j].Vn;
            Pn = Fn / Dn; Msh.tempA[k].an = Dn * funcScheme(Pn) + std::max(-Fn, 0.0);

            // Coefficients A
            a0 = Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp * Msh.Vp[i][j] / dt;
            Msh.matA[k].aw = - beta * Msh.tempA[k].aw; Msh.matA[k].ae = - beta * Msh.tempA[k].ae;
            Msh.matA[k].as = - beta * Msh.tempA[k].as; Msh.matA[k].an = - beta * Msh.tempA[k].an; 
            Msh.matA[k].ap = a0 - Msh.matA[k].ae - Msh.matA[k].aw - Msh.matA[k].an - Msh.matA[k].as;
            
            // Coefficients B 
            Msh.bp[k] = Msh.sPhi[i][j] * Msh.Vp[i][j] + a0 * Msh.vPhi[i][j] + (1 - beta) * (Msh.vPhi[i-1][j] * Msh.tempA[k].aw + Msh.vPhi[i+1][j] * Msh.tempA[k].ae + Msh.vPhi[i][j-1] * Msh.tempA[k].as + Msh.vPhi[i][j+1] * Msh.tempA[k].an - Msh.vPhi[i][j] * (Msh.tempA[k].ae + Msh.tempA[k].aw + Msh.tempA[k].an + Msh.tempA[k].as));

        }
    }
    
}
