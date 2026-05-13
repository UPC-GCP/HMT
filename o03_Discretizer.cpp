// Imports
/* #include "json/value.h" */
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
/* #include "exprtk.hpp" */
#include "o01_Material.h"
#include "o02_Mesh.h"
#include "o03_Discretizer.h"
#include "o09_ExpressionParser.h"


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
	int interiorNodes = 1; for (int NVal : Msh.N){interiorNodes *= (NVal-2);}
        std::vector<double> dtNew(interiorNodes, 0); double dtMin{}; int k{};
        for (size_t i = 1; i < Msh.N[0]-1; i++){
            for (size_t j = 1; j < Msh.N[1]-1; j++){
                k = (i-1) * (Msh.N[1]-2) + (j-1);
		dtNew[k] = 0.5 * pow(Msh.ndelta[0][i], 2) * pow(Msh.ndelta[1][j], 2) / ((Mat.vMat[Msh.nMat[i][j]].gamma / (Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp)) * (pow(Msh.ndelta[0][i], 2) + pow(Msh.ndelta[1][j], 2)));
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
        std::cerr << "Error: Invalid spatial discretization scheme " << spatScheme << "\n";
    }

}


void Discretizer::setBoundaryConditions(Material& Mat, Mesh& Msh, ExpressionParser& Prs, double t){ // DIFFUSION, NO CONVECTION

    // Boundary Conditions
    std::vector<int> Pos0, Pos1; Pos0.resize(Msh.N.size()); Pos1.resize(Msh.N.size()); double lamb; int j;
    for (Boundary bC : Msh.boundaryConditions){

        // Positions (nD)
        for (int i = 0; i < Msh.N.size(); i++){	
            Pos0[i] = std::lower_bound(Msh.Nodes[i].begin(), Msh.Nodes[i].end(), bC.x0[i] - epsFind) - Msh.Nodes[i].begin();
            Pos1[i] = std::lower_bound(Msh.Nodes[i].begin(), Msh.Nodes[i].end(), bC.x1[i] - epsFind) - Msh.Nodes[i].begin();
        }

        // Boundaries (Non-nD)
        if (bC.type == 0){
            
            // Update Value
            if (bC.bUpdate){
                bC.value = Prs.evaluateTime(bC.iExpr, t);
            }

            // Dirichlet
            if (Pos0[0] == Pos1[0]){

                // xBoundary
                for (int i = Pos0[1]; i <= Pos1[1]; i++){
                    // Ignore Corners
                    if (i == 0 || i == Msh.N[1]-1){continue;}

                    // Value
                    Msh.vPhi[Pos0[0]][i] = bC.value;

                    // Coefficients
                    j = Pos0[0] * Msh.N[1] + i;
                    Msh.matA[j].ap = 1;
                    Msh.bp[j] = bC.value;

                    /* // Control */
                    /* if (bIgnore){Msh.nIgnore.push_back({Pos0[0], i});} */
                }

            } else if (Pos0[1] == Pos1[1]){
                
                // yBoundary
                for (int i = Pos0[0]; i <= Pos1[0]; i++){
                    // Ignore Corners
                    if (i == 0 || i == Msh.N[0]-1){continue;}
                    
                    // Value
                    Msh.vPhi[i][Pos0[1]] = bC.value;

                    // Coefficients
                    j = i * Msh.N[1] + Pos0[1];
                    Msh.matA[j].ap = 1;
                    Msh.bp[j] = bC.value;

                    /* // Control */
                    /* if (bIgnore){Msh.nIgnore.push_back({i, Pos0[1]});} */
                }

            }

        } else if (bC.type == 1){

            // Neumann
            if (Pos0[0] == Pos1[0]){
                
                // xBoundary
                if (bC.side == 0){

                    // West Boundary // ae, ap, bp
                    for (int i = Pos0[1]; i <= Pos1[1]; i++){
                        // Ignore Corners
                        if (i == 0 || i == Msh.N[1]-1){continue;}

                        // Thermal Conductivity
                        lamb = Mat.vMat[Msh.nMat[Pos0[0]][i]].gamma;

                        // Value
                        Msh.vPhi[Pos0[0]][i] = bC.value * Msh.nd[0][Pos0[0]] / lamb + Msh.vPhi[Pos0[0]+1][i];

                        // Coefficients
                        j = Pos0[0] * Msh.N[1] + i;
                        Msh.matA[j].ap = 1;
                        Msh.bp[j] = Msh.vPhi[Pos0[0]][i];

  			/* // Control (Testing) */
	                /* if (bIgnore){Msh.nIgnore.push_back({Pos0[0], i});} */

                    }

                } else if (bC.side == 1){

                    // East Boundary // aw, ap, bp
                    for (int i = Pos0[1]; i <= Pos1[1]; i++){
                        // Ignore Corners
                        if (i == 0 || i == Msh.N[1]-1){continue;}
                        
                        // Thermal Conductivity
                        lamb = Mat.vMat[Msh.nMat[Pos0[0]][i]].gamma;

                        // Value
                        Msh.vPhi[Pos0[0]][i] = bC.value * Msh.nd[0][Pos0[0]-1] / lamb + Msh.vPhi[Pos0[0]-1][i];

                        // Coefficients
                        j = Pos0[0] * Msh.N[1] + i;
                        Msh.matA[j].ap = 1;
                        Msh.bp[j] = Msh.vPhi[Pos0[0]][i];

  			/* // Control (Testing) */
	                /* if (bIgnore){Msh.nIgnore.push_back({Pos0[0], i});} */

                    }

                } else {std::cerr << "Boundary side not specified correcly.\n";}

            } else if (Pos0[1] == Pos1[1]){

                // yBoundary
                if (bC.side == 0){

                    // South Boundary // an, ap, bp
                    for (int i = Pos0[0]; i <= Pos1[0]; i++){
                        // Ignore Corners
                        if (i == 0 || i == Msh.N[0]-1){continue;}

                        // Thermal Conductivity
                        lamb = Mat.vMat[Msh.nMat[i][Pos0[1]]].gamma;

                        // Value
                        Msh.vPhi[i][Pos0[1]] = bC.value * Msh.nd[1][Pos0[1]] / lamb + Msh.vPhi[i][Pos0[1]+1];

                        // Coefficients
                        j = i * Msh.N[1] + Pos0[1];
                        Msh.matA[j].ap = 1;
                        Msh.bp[j] = Msh.vPhi[i][Pos0[1]];

  			/* // Control (Testing) */
	                /* if (bIgnore){Msh.nIgnore.push_back({i, Pos0[1]});} */

                    }

                } else if (bC.side == 1){

                    // North Boundary // as, ap, bp
                    for (int i = Pos0[0]; i <= Pos1[0]; i++){
                        // Ignore Corners
                        if (i == 0 || i == Msh.N[0]-1){continue;}
                        
                        // Thermal Conductivity
                        lamb = Mat.vMat[Msh.nMat[i][Pos0[1]]].gamma;

                        // Value
                        Msh.vPhi[i][Pos0[1]] = bC.value * Msh.nd[1][Pos0[1]-1] / lamb + Msh.vPhi[i][Pos0[1]-1];

                        // Coefficients
                        j = i * Msh.N[1] + Pos0[1];
                        Msh.matA[j].ap = 1;
                        Msh.bp[j] = Msh.vPhi[i][Pos0[1]];

  			/* // Control (Testing) */
	                /* if (bIgnore){Msh.nIgnore.push_back({i, Pos0[1]});} */

                    }

                } else {std::cerr << "Boundary side not specified correcly.\n";}

            }

        } else if (bC.type == 2){
            
            // Convection
            if (Pos0[0] == Pos1[0]){
                
                // xBoundary
                if (bC.side == 0){

                    // West Boundary // ae, ap, bp
                    for (int i = Pos0[1]; i <= Pos1[1]; i++){
                        // Ignore Corners
                        if (i == 0 || i == Msh.N[1]-1){continue;}

                        // Thermal Conductivity
                        lamb = Mat.vMat[Msh.nMat[Pos0[0]][i]].gamma;

                        // Value
                        Msh.vPhi[Pos0[0]][i] = (bC.alpha * bC.value + lamb * Msh.vPhi[Pos0[0]+1][i] / Msh.nd[0][Pos0[0]]) / (lamb/Msh.nd[0][Pos0[0]] + bC.alpha);

                        // Coefficients
                        j = Pos0[0] * Msh.N[1] + i;
                        Msh.matA[j].ap = 1;
                        Msh.bp[j] = Msh.vPhi[Pos0[0]][i];
			
			/* // Control (Testing) */
	                /* if (bIgnore){Msh.nIgnore.push_back({Pos0[0], i});} */

                    }

                } else if (bC.side == 1){

                    // East Boundary // aw, ap, bp
                    for (int i = Pos0[1]; i <= Pos1[1]; i++){
                        // Ignore Corners
                        if (i == 0 || i == Msh.N[1]-1){continue;}
                        
                        // Thermal Conductivity
                        lamb = Mat.vMat[Msh.nMat[Pos0[0]][i]].gamma;

                        // Value
                        Msh.vPhi[Pos0[0]][i] = (bC.alpha * bC.value + lamb * Msh.vPhi[Pos0[0]-1][i] / Msh.nd[0][Pos0[0]-1]) / (lamb/Msh.nd[0][Pos0[0]-1] + bC.alpha);

                        // Coefficients
                        j = Pos0[0] * Msh.N[1] + i;
                        Msh.matA[j].ap = 1;
                        Msh.bp[j] = Msh.vPhi[Pos0[0]][i];
			
			/* // Control (Testing) */
	                /* if (bIgnore){Msh.nIgnore.push_back({Pos0[0], i});} */

                    }

                } else {std::cerr << "Boundary side not specified correcly.\n";}

            } else if (Pos0[1] == Pos1[1]){
                
                // yBoundary
                if (bC.side == 0){

                    // South Boundary // an, ap, bp
                    for (int i = Pos0[0]; i <= Pos1[0]; i++){
                        // Ignore Corners
                        if (i == 0 || i == Msh.N[0]-1){continue;}

                        // Thermal Conductivity
                        lamb = Mat.vMat[Msh.nMat[i][Pos0[1]]].gamma;

                        // Value
                        Msh.vPhi[i][Pos0[1]] = (bC.alpha * bC.value + lamb * Msh.vPhi[i][Pos0[1]+1] / Msh.nd[1][Pos0[1]]) / (lamb/Msh.nd[1][Pos0[1]] + bC.alpha);

                        // Coefficients
                        j = i * Msh.N[1] + Pos0[1];
                        Msh.matA[j].ap = 1;
                        Msh.bp[j] = Msh.vPhi[i][Pos0[1]];

  			/* // Control (Testing) */
	                /* if (bIgnore){Msh.nIgnore.push_back({i, Pos0[1]});} */

                    }

                } else if (bC.side == 1){

                    // North Boundary // as, ap, bp
                    for (int i = Pos0[0]; i <= Pos1[0]; i++){
                        // Ignore Corners
                        if (i == 0 || i == Msh.N[0]-1){continue;}
                        
                        // Thermal Conductivity
                        lamb = Mat.vMat[Msh.nMat[i][Pos0[1]]].gamma;

                        // Value
                        Msh.vPhi[i][Pos0[1]] = (bC.alpha * bC.value + lamb * Msh.vPhi[i][Pos0[1]-1] / Msh.nd[1][Pos0[1]-1]) / (lamb/Msh.nd[1][Pos0[1]-1] + bC.alpha);

                        // Coefficients
                        j = i * Msh.N[1] + Pos0[1];
                        Msh.matA[j].ap = 1;
                        Msh.bp[j] = Msh.vPhi[i][Pos0[1]];
  			
			/* // Control (Testing) */
	                /* if (bIgnore){Msh.nIgnore.push_back({i, Pos0[1]});} */
                        
                    }

                } else {std::cerr << "Boundary side not specified correcly.\n";}

            }

        }
        
    }

    // Corners
    Msh.vPhi[0][0] = 0.5 * (Msh.vPhi[1][0] + Msh.vPhi[0][1]);
    Msh.vPhi[0][Msh.N[1]-1] = 0.5 * (Msh.vPhi[1][Msh.N[1]-1] + Msh.vPhi[0][Msh.N[1]-2]);
    Msh.vPhi[Msh.N[0]-1][0] = 0.5 * (Msh.vPhi[Msh.N[0]-2][0] + Msh.vPhi[Msh.N[0]-1][1]);
    Msh.vPhi[Msh.N[0]-1][Msh.N[1]-1] = 0.5 * (Msh.vPhi[Msh.N[0]-2][Msh.N[1]-1] + Msh.vPhi[Msh.N[0]-1][Msh.N[1]-2]);

    // Corners
    int l{}, m{}, n{};

    l = 0; m = 0; n = l * Msh.N[1] + m; // SW
    Msh.vPhi[l][m] = 0.5 * (Msh.vPhi[l+1][m] + Msh.vPhi[l][m+1]);
    Msh.matA[n].ap = 2; Msh.matA[n].ae = -1; Msh.matA[n].an = -1;

    l = 0; m = Msh.N[1]-1; n = l * Msh.N[1] + m; // NW
    Msh.vPhi[l][m] = 0.5 * (Msh.vPhi[l+1][m] + Msh.vPhi[l][m-1]);
    Msh.matA[n].ap = 2; Msh.matA[n].ae = -1; Msh.matA[n].as = -1;

    l = Msh.N[0]-1; m = 0; n = l * Msh.N[1] + m; // SE
    Msh.vPhi[l][m] = 0.5 * (Msh.vPhi[l-1][m] + Msh.vPhi[l][m+1]);
    Msh.matA[n].ap = 2; Msh.matA[n].aw = -1; Msh.matA[n].an = -1;

    l = Msh.N[0]-1; m = Msh.N[1]-1; n = l * Msh.N[1] + m; // NE
    Msh.vPhi[l][m] = 0.5 * (Msh.vPhi[l-1][m] + Msh.vPhi[l][m-1]);
    Msh.matA[n].ap = 2; Msh.matA[n].aw = -1; Msh.matA[n].as = -1;

    /* // Control */
    /* if (bIgnore){Msh.nIgnore.push_back({0, 0}); Msh.nIgnore.push_back({0, Msh.N[1]-1}); Msh.nIgnore.push_back({Msh.N[0]-1, 0}); Msh.nIgnore.push_back({Msh.N[0]-1, Msh.N[1]-1});} */
    /* bIgnore = false; */

}


void Discretizer::newSetBoundaries(Material& Mat, Mesh& Msh, ExpressionParser& Prs, double t){

    // Control
    /* double Fw{}, Fe{}, Fs{}, Fn{}, Dw{}, De{}, Ds{}, Dn{}, Pw{}, Pe{}, Ps{}, Pn{}; */
    std::vector<int> Pos0{}, Pos1{}; Pos0.resize(Msh.N.size()); Pos1.resize(Msh.N.size()); double gamm{}, Dp{}, Fp{}, Pp{}, sumConv{}, sumDiff{}; int k{};

    // Boundary Conditions
    for (Boundary bC : Msh.boundaryConditions){

        // Positions (nD)
        for (size_t i = 0; i < Msh.N.size(); i++){
            Pos0[i] = std::lower_bound(Msh.Nodes[i].begin(), Msh.Nodes[i].end(), bC.x0[i] - epsFind) - Msh.Nodes[i].begin();
            Pos1[i] = std::lower_bound(Msh.Nodes[i].begin(), Msh.Nodes[i].end(), bC.x1[i] - epsFind) - Msh.Nodes[i].begin();
        }

        // Boundaries (Non-nD)
        if (bC.type == 0){

            // Update Value
            if (bC.bUpdate){
                bC.value = Prs.evaluateTime(bC.iExpr, t);
            }

            // Dirichlet
            if (Pos0[0] == Pos1[0]){

                // xBoundary
                for (int i = Pos0[1]; i <= Pos1[1]; i++){
                    // Ignore Corners
                    if (i == 0 || i == Msh.N[1]-1){continue;}

                    // Value
                    Msh.vPhi[Pos0[0]][i] = bC.value;

                    // Coefficients
                    k = Pos0[0] * Msh.N[1] + i;
                    Msh.matA[k].ap = 1;
                    Msh.bp[k] = bC.value;

                    /* // Control */
                    /* if (bIgnore){Msh.nIgnore.push_back({Pos0[0], i});} */
                }

            } else if (Pos0[1] == Pos1[1]){
                
                // yBoundary
                for (int i = Pos0[0]; i <= Pos1[0]; i++){
                    // Ignore Corners
                    if (i == 0 || i == Msh.N[0]-1){continue;}
                    
                    // Value
                    Msh.vPhi[i][Pos0[1]] = bC.value;

                    // Coefficients
                    k = i * Msh.N[1] + Pos0[1];
                    Msh.matA[k].ap = 1;
                    Msh.bp[k] = bC.value;

                    /* // Control */
                    /* if (bIgnore){Msh.nIgnore.push_back({i, Pos0[1]});} */
                }

            }
            
        } else if (bC.type == 1){

            // Neumann
            if (Pos0[0] == Pos1[0]){
                
                // xBoundary
                if (bC.side == 0){

                    // West Boundary // ae, ap, bp
                    for (int i = Pos0[1]; i <= Pos1[1]; i++){
                        // Ignore Corners
                        if (i == 0 || i == Msh.N[1]-1){continue;}

                        // Thermal Conductivity
                        gamm = Mat.vMat[Msh.nMat[Pos0[0]][i]].gamma;

                        // Value
                        /* Msh.vPhi[Pos0[0]][i] = bC.value * Msh.nd[0][Pos0[0]] / gamm + Msh.vPhi[Pos0[0]+1][i]; */

                        // Index
                        k = Pos0[0] * Msh.N[1] + i;

                        // Coefficients Conv-Diff
                        Dp = - gamm * Msh.Se[Pos0[0]][i] / Msh.nd[0][Pos0[0]]; Fp = Mat.vMat[Msh.nMat[Pos0[0]][i]].rho * Msh.Se[Pos0[0]][i] * Msh.vField.vConv[0][Pos0[0]+1][i]; Pp = Fp / Dp;

                        // Coefficients A
                        Msh.matA[k].ae = - beta * (Dp * funcScheme(std::abs(Pp)) + std::max(-Fp, 0.0));
                        Msh.matA[k].ap = - Msh.matA[k].ae;
                        
                        // Coefficients B
                        /* sumConv = Mat.vMat[Msh.nMat[Pos0[0]][i]].rho * (Msh.vField.vConv[0][Pos0[0]+1][i] * Msh.vPhi[Pos0[0]+1][i] * Msh.Se[Pos0[0]][i] + Msh.vField.vConv[0][Pos0[0]+1][i] * Msh.Se[Pos0[0]][i] * Msh.vPhi[Pos0[0]][i]); */
                        sumConv = Mat.vMat[Msh.nMat[Pos0[0]][i]].rho * Msh.vField.vConv[0][Pos0[0]+1][i] * Msh.Se[Pos0[0]][i] * Msh.vPhi[Pos0[0]+1][i];
                        sumDiff = gamm * (Msh.Se[Pos0[0]][i] * Msh.vPhi[Pos0[0]+1][i] / Msh.nd[0][Pos0[0]] - Msh.Se[Pos0[0]][i] * Msh.vPhi[Pos0[0]][i]) / Msh.nd[0][Pos0[0]];
                        Msh.bp[k] = bC.value * Msh.Se[Pos0[0]][i] + (1 - beta) * (sumConv + sumDiff); // Not entirely sure if this should have q_N * S or just q_N

  			            /* // Control (Testing) */
	                    /* if (bIgnore){Msh.nIgnore.push_back({Pos0[0], i});} */

                    }

                } else if (bC.side == 1){

                    // East Boundary // aw, ap, bp
                    for (int i = Pos0[1]; i <= Pos1[1]; i++){
                        // Ignore Corners
                        if (i == 0 || i == Msh.N[1]-1){continue;}
                        
                        // Thermal Conductivity
                        gamm = Mat.vMat[Msh.nMat[Pos0[0]][i]].gamma;

                        // Value
                        /* Msh.vPhi[Pos0[0]][i] = bC.value * Msh.nd[0][Pos0[0]-1] / gamm + Msh.vPhi[Pos0[0]-1][i]; */

                        // Index
                        k = Pos0[0] * Msh.N[1] + i;

                        // Coefficients Conv-Diff
                        Dp = - gamm * Msh.Sw[Pos0[0]][i] / Msh.nd[0][Pos0[0]-1]; Fp = Mat.vMat[Msh.nMat[Pos0[0]][i]].rho * Msh.Sw[Pos0[0]][i] * Msh.vField.vConv[0][Pos0[0]-1][i]; Pp = Fp / Dp;

                        // Coefficients A
                        Msh.matA[k].aw = - beta * (Dp * funcScheme(std::abs(Pp)) + std::max(Fp, 0.0));
                        Msh.matA[k].ap = - Msh.matA[k].aw; 

                        // Coefficients B
                        sumConv = Mat.vMat[Msh.nMat[Pos0[0]][i]].rho * Msh.vField.vConv[0][Pos0[0]-1][i] * Msh.Sw[Pos0[0]][i] * Msh.vPhi[Pos0[0]][i];
                        sumDiff = gamm * (Msh.Sw[Pos0[0]][i] * Msh.vPhi[Pos0[0]-1][i] / Msh.nd[0][Pos0[0]-1] - Msh.Sw[Pos0[0]][i] * Msh.vPhi[Pos0[0]][i] / Msh.nd[0][Pos0[0]-1]);
                        Msh.bp[k] = bC.value * Msh.Sw[Pos0[0]][i] + (1 - beta) * (sumConv + sumDiff);

          			/* // Control (Testing) */
	                /* if (bIgnore){Msh.nIgnore.push_back({Pos0[0], i});} */

                    }

                } else {std::cerr << "Boundary side not specified correcly.\n";}

            } else if (Pos0[1] == Pos1[1]){

                // yBoundary
                if (bC.side == 0){

                    // South Boundary // an, ap, bp
                    for (int i = Pos0[0]; i <= Pos1[0]; i++){
                        // Ignore Corners
                        if (i == 0 || i == Msh.N[0]-1){continue;}

                        // Thermal Conductivity
                        gamm = Mat.vMat[Msh.nMat[i][Pos0[1]]].gamma;

                        // Value
                        /* Msh.vPhi[i][Pos0[1]] = bC.value * Msh.nd[1][Pos0[1]] / gamm + Msh.vPhi[i][Pos0[1]+1]; */

                        // Index
                        k = i * Msh.N[1] + Pos0[1];
                        
                        // Coefficients Conv-Diff
                        Dp = - gamm * Msh.Sn[i][Pos0[0]] / Msh.nd[1][Pos0[0]]; Fp = Mat.vMat[Msh.nMat[i][Pos0[0]]].rho * Msh.Sn[i][Pos0[0]] * Msh.vField.vConv[1][i][Pos0[0]+1]; Pp = Fp / Dp;
                        
                        // Coefficients A
                        Msh.matA[k].an = - beta * (Dp * funcScheme(std::abs(Pp)) + std::max(-Fp, 0.0));
                        Msh.matA[k].ap = - Msh.matA[k].an;

                        // Coefficients B
                        sumConv = Mat.vMat[Msh.nMat[i][Pos0[0]]].rho * Msh.vField.vConv[1][i][Pos0[0]+1] * Msh.Sn[i][Pos0[0]] * Msh.vPhi[i][Pos0[0]+1];
                        sumDiff = gamm * (Msh.Sn[i][Pos0[0]] * Msh.vPhi[i][Pos0[0]+1] / Msh.nd[1][Pos0[0]] - Msh.Sn[i][Pos0[0]] * Msh.vPhi[i][Pos0[0]] / Msh.nd[1][Pos0[0]]);
                        Msh.bp[k] = bC.value * Msh.Sn[i][Pos0[0]] + (1 - beta) * (sumConv + sumDiff); // VALIDATE IF AREAS HERE HAVE A VALUE OR IF THEY'RE JUST 0 - IF 0 CHANGE

  			        /* // Control (Testing) */
	                /* if (bIgnore){Msh.nIgnore.push_back({i, Pos0[1]});} */

                    }

                } else if (bC.side == 1){

                    // North Boundary // as, ap, bp
                    for (int i = Pos0[0]; i <= Pos1[0]; i++){
                        // Ignore Corners
                        if (i == 0 || i == Msh.N[0]-1){continue;}
                        
                        // Thermal Conductivity
                        gamm = Mat.vMat[Msh.nMat[i][Pos0[1]]].gamma;

                        // Value
                        Msh.vPhi[i][Pos0[1]] = bC.value * Msh.nd[1][Pos0[1]-1] / gamm + Msh.vPhi[i][Pos0[1]-1];

                        // Index
                        k = i * Msh.N[1] + Pos0[1];

                        // Coefficients Conv-Diff
                        Dp = - gamm * Msh.Ss[i][Pos0[0]] / Msh.nd[1][Pos0[0]-1]; Fp = Mat.vMat[Msh.nMat[i][Pos0[0]]].rho * Msh.Ss[i][Pos0[0]] * Msh.vField.vConv[1][i][Pos0[0]-1]; Pp = Fp / Dp;

                        // Coefficients A
                        Msh.matA[k].as = - beta * (Dp * funcScheme(std::abs(Pp)) + std::max(Fp, 0.0));
                        Msh.matA[k].ap = - Msh.matA[k].as;

                        // Coefficients B
                        sumConv = Mat.vMat[Msh.nMat[i][Pos0[0]]].rho * Msh.vField.vConv[1][i][Pos0[0]-1] * Msh.Ss[i][Pos0[0]] * Msh.vPhi[i][Pos0[0]];
                        sumDiff = gamm * (Msh.Ss[i][Pos0[0]] * Msh.vPhi[i][Pos0[0]-1] / Msh.nd[1][Pos0[0]-1] - Msh.Ss[i][Pos0[0]] * Msh.vPhi[i][Pos0[0]] / Msh.nd[1][Pos0[0]-1]);
                        Msh.bp[k] = bC.value * Msh.Ss[i][Pos0[0]] + (1 - beta) * (sumConv + sumDiff);


            /* Ds = - gammas * Msh.Ss[i][j] / Msh.nd[1][j-1]; Fs = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Ss[i][j] * Msh.vField.vConv[1][i][j-1]; Ps = Fs / Ds; */
            /* Msh.matA[k].as = - beta * (Ds * funcScheme(std::abs(Ps)) + std::max(Fs, 0.0)); */






  	        		/* // Control (Testing) */
	                /* if (bIgnore){Msh.nIgnore.push_back({i, Pos0[1]});} */

                    }

                } else {std::cerr << "Boundary side not specified correcly.\n";}

            }



        } else if (bC.type == 2){

        } else {std::cerr << "Boundary side not specified correctly.\n";}

    }

}


void Discretizer::setCoefficients(Material& Mat, Mesh& Msh){ // DIFFUSION, NO CONVECTION
    
    // Control
    double lambw, lambe, lambs, lambn; int k;

    // Interior Nodes (Non-nD)
    for (size_t i = 1; i < Msh.N[0]-1; i++){
        for (size_t j = 1; j < Msh.N[1]-1; j++){

            // Harmonic Mean
            lambw = calcHarmonicMean(Msh.nd[0][i-1], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i-1][j]].gamma}, {Msh.ndelta[0][i], Msh.ndelta[0][i-1]});
            lambe = calcHarmonicMean(Msh.nd[0][i], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i+1][j]].gamma}, {Msh.ndelta[0][i], Msh.ndelta[0][i+1]});
            lambs = calcHarmonicMean(Msh.nd[1][j-1], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i][j-1]].gamma}, {Msh.ndelta[1][j], Msh.ndelta[1][j-1]});
            lambn = calcHarmonicMean(Msh.nd[1][j], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i][j+1]].gamma}, {Msh.ndelta[1][j], Msh.ndelta[1][j+1]});

            // Index
            k = i * Msh.N[1] + j;
            
            // Coefficients A
            Msh.matA[k].aw = - beta * lambw * Msh.Sw[i][j] / Msh.nd[0][i-1];
            Msh.matA[k].ae = - beta * lambe * Msh.Se[i][j] / Msh.nd[0][i];
            Msh.matA[k].as = - beta * lambs * Msh.Ss[i][j] / Msh.nd[1][j-1];
            Msh.matA[k].an = - beta * lambn * Msh.Sn[i][j] / Msh.nd[1][j];
            Msh.matA[k].ap = Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp * Msh.Vp[i][j] / dt - Msh.matA[k].aw - Msh.matA[k].ae - Msh.matA[k].as - Msh.matA[k].an;

           // Coefficients B
            Msh.bp[k] = Msh.sPhi[i][j] * Msh.Vp[i][j] + Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp * Msh.vPhi[i][j] * Msh.Vp[i][j] / dt + (1 - beta) * (lambw * Msh.Sw[i][j] * Msh.vPhi[i-1][j] / Msh.nd[0][i-1] + lambe * Msh.Se[i][j] * Msh.vPhi[i+1][j] / Msh.nd[0][i] + lambs * Msh.Ss[i][j] * Msh.vPhi[i][j-1] / Msh.nd[1][j-1] + lambn * Msh.Sn[i][j] * Msh.vPhi[i][j+1] / Msh.nd[1][j] - (lambw * Msh.Sw[i][j] / Msh.nd[0][i-1] + lambe * Msh.Se[i][j] / Msh.nd[0][i] + lambs * Msh.Ss[i][j] / Msh.nd[1][j-1] + lambn * Msh.Sn[i][j] / Msh.nd[1][j]) * Msh.vPhi[i][j]);

        }
    }

}


void Discretizer::newSetCoefficients(Material& Mat, Mesh& Msh){

    // Control
    double gammaw, gammae, gammas{}, gamman{}, sumConv{}, sumDiff{}; int k{};
    double Fw{}, Fe{}, Fs{}, Fn{}, Dw{}, De{}, Ds{}, Dn{}, Pw{}, Pe{}, Ps{}, Pn{};

    // Interior Nodes Loop
    for (size_t i = 1; i < Msh.N[0]-1; i++){
        for (size_t j = 1; j < Msh.N[1]-1; j++){
        
            // Harmonic Mean
            gammaw = calcHarmonicMean(Msh.nd[0][i-1], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i-1][j]].gamma}, {Msh.ndelta[0][i], Msh.ndelta[0][i-1]});
            gammae = calcHarmonicMean(Msh.nd[0][i], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i+1][j]].gamma}, {Msh.ndelta[0][i], Msh.ndelta[0][i+1]});
            gammas = calcHarmonicMean(Msh.nd[1][j-1], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i][j-1]].gamma}, {Msh.ndelta[1][j], Msh.ndelta[1][j-1]});
            gamman = calcHarmonicMean(Msh.nd[1][j], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i][j+1]].gamma}, {Msh.ndelta[1][j], Msh.ndelta[1][j+1]});

            // Index
            k = i * Msh.N[1] + j;

            // Coefficients Conv-Diff
            De = - gammae * Msh.Se[i][j] / Msh.nd[0][i]; Fe = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Se[i][j] * Msh.vField.vConv[0][i+1][j]; Pe = Fe / De;
            Dw = - gammaw * Msh.Sw[i][j] / Msh.nd[0][i-1]; Fw = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sw[i][j] * Msh.vField.vConv[0][i-1][j]; Pw = Fw / Dw;
            Dn = - gamman * Msh.Sn[i][j] / Msh.nd[1][j]; Fn = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sn[i][j] * Msh.vField.vConv[1][i][j+1]; Pn = Fn / Dn;
            Ds = - gammas * Msh.Ss[i][j] / Msh.nd[1][j-1]; Fs = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Ss[i][j] * Msh.vField.vConv[1][i][j-1]; Ps = Fs / Ds;

            // Coefficients A
            Msh.matA[k].ae = - beta * (De * funcScheme(std::abs(Pe)) + std::max(-Fe, 0.0));
            Msh.matA[k].aw = - beta * (Dw * funcScheme(std::abs(Pw)) + std::max(Fw, 0.0));
            Msh.matA[k].an = - beta * (Dn * funcScheme(std::abs(Pn)) + std::max(-Fn, 0.0));
            Msh.matA[k].as = - beta * (Ds * funcScheme(std::abs(Ps)) + std::max(Fs, 0.0));
            Msh.matA[k].ap = Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp * Msh.Vp[i][j] / dt - Msh.matA[k].ae - Msh.matA[k].aw - Msh.matA[k].an - Msh.matA[k].as - Msh.sPhi[i][j] * Msh.Vp[i][j];

            // Coefficients B
            sumConv = Mat.vMat[Msh.nMat[i][j]].rho * (Msh.vField.vConv[0][i+1][j] * Msh.vPhi[i+1][j] * Msh.Se[i][j] - Msh.vField.vConv[0][i-1][j] * Msh.vPhi[i-1][j] * Msh.Sw[i][j] + Msh.vField.vConv[1][i][j+1] * Msh.vPhi[i][j+1] * Msh.Sn[i][j] - Msh.vField.vConv[1][i][j-1] * Msh.vPhi[i][j-1] * Msh.Ss[i][j] + (Msh.vField.vConv[0][i+1][j] * Msh.Se[i][j] - Msh.vField.vConv[0][i-1][j] * Msh.Sw[i][j] + Msh.vField.vConv[1][i][j+1] * Msh.Sn[i][j] - Msh.vField.vConv[1][i][j-1] * Msh.Ss[i][j]) * Msh.vPhi[i][j]);
            sumDiff = (gammaw * Msh.Sw[i][j] * Msh.vPhi[i-1][j] / Msh.nd[0][i-1] + gammae * Msh.Se[i][j] * Msh.vPhi[i+1][j] / Msh.nd[0][i] + gammas * Msh.Ss[i][j] * Msh.vPhi[i][j-1] / Msh.nd[1][j-1] + gamman * Msh.Sn[i][j] * Msh.vPhi[i][j+1] / Msh.nd[1][j] - (gammaw * Msh.Sw[i][j] / Msh.nd[0][i-1] + gammae * Msh.Se[i][j] / Msh.nd[0][i] + gammas * Msh.Ss[i][j] / Msh.nd[1][j-1] + gamman * Msh.Sn[i][j] / Msh.nd[1][j]) * Msh.vPhi[i][j]);
            Msh.bp[k] = Msh.sPhi[i][j] * Msh.Vp[i][j] + Mat.vMat[Msh.nMat[i][j]].rho * Msh.Vp[i][j] * Msh.vPhi[i][j] / dt + (1 - beta) * (sumConv + sumDiff);

        }
    }

    
}


void Discretizer::setRHS(Material& Mat, Mesh& Msh){ // DIFFUSION, NO CONVECTION

    // Control
    double lambw, lambe, lambs, lambn; int k;

    // Interior Nodes (Non-nD)
    for (size_t i = 1; i < Msh.N[0]-1; i++){
        for (size_t j = 1; j < Msh.N[1]-1; j++){
            
            // Harmonic Mean
            lambw = calcHarmonicMean(Msh.nd[0][i-1], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i-1][j]].gamma}, {Msh.ndelta[0][i], Msh.ndelta[0][i-1]});
            lambe = calcHarmonicMean(Msh.nd[0][i], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i+1][j]].gamma}, {Msh.ndelta[0][i], Msh.ndelta[0][i+1]});
            lambs = calcHarmonicMean(Msh.nd[1][j-1], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i][j-1]].gamma}, {Msh.ndelta[1][j], Msh.ndelta[1][j-1]});
            lambn = calcHarmonicMean(Msh.nd[1][j], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i][j+1]].gamma}, {Msh.ndelta[1][j], Msh.ndelta[1][j+1]});

            // Index
            k = i * Msh.N[1] + j;

            // Coefficients B
            Msh.bp[k] = Msh.sPhi[i][j] * Msh.Vp[i][j] + Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp * Msh.Vp[i][j] * Msh.vPhi[i][j] / dt + (1 - beta) * (lambw * Msh.Sw[i][j] * Msh.vPhi[i-1][j] / Msh.nd[0][i-1] + lambe * Msh.Se[i][j] * Msh.vPhi[i+1][j] / Msh.nd[0][i] + lambs * Msh.Ss[i][j] * Msh.vPhi[i][j-1] / Msh.nd[1][j-1] + lambn * Msh.Sn[i][j] * Msh.vPhi[i][j+1] / Msh.nd[1][j] - (lambw * Msh.Sw[i][j] / Msh.nd[0][i-1] + lambe * Msh.Se[i][j] / Msh.nd[0][i] + lambs * Msh.Ss[i][j] / Msh.nd[1][j-1] + lambn * Msh.Sn[i][j] / Msh.nd[1][j]) * Msh.vPhi[i][j]);

        }
    }

}


void Discretizer::newSetRHS(Material& Mat, Mesh& Msh){

    // Control
    double gammaw{}, gammae{}, gammas{}, gamman{}, sumConv{}, sumDiff{}; int k;
    double Fw{}, Fe{}, Fs{}, Fn{}, Dw{}, De{}, Ds{}, Dn{}, Pw{}, Pe{}, Ps{}, Pn{};

    // Interior Nodes (Non-nD)
    for (size_t i = 1; i < Msh.N[0]-1; i++){
        for (size_t j = 1; j < Msh.N[1]-1; j++){
            
            // Harmonic Mean
            gammaw = calcHarmonicMean(Msh.nd[0][i-1], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i-1][j]].gamma}, {Msh.ndelta[0][i], Msh.ndelta[0][i-1]});
            gammae = calcHarmonicMean(Msh.nd[0][i], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i+1][j]].gamma}, {Msh.ndelta[0][i], Msh.ndelta[0][i+1]});
            gammas = calcHarmonicMean(Msh.nd[1][j-1], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i][j-1]].gamma}, {Msh.ndelta[1][j], Msh.ndelta[1][j-1]});
            gamman = calcHarmonicMean(Msh.nd[1][j], {Mat.vMat[Msh.nMat[i][j]].gamma, Mat.vMat[Msh.nMat[i][j+1]].gamma}, {Msh.ndelta[1][j], Msh.ndelta[1][j+1]});

            // Index
            k = i * Msh.N[1] + j;

            // Coefficients Conv-Diff
            De = - gammae * Msh.Se[i][j] / Msh.nd[0][i]; Fe = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Se[i][j] * Msh.vField.vConv[0][i+1][j]; Pe = Fe / De;
            Dw = - gammaw * Msh.Sw[i][j] / Msh.nd[0][i-1]; Fw = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sw[i][j] * Msh.vField.vConv[0][i-1][j]; Pw = Fw / Dw;
            Dn = - gamman * Msh.Sn[i][j] / Msh.nd[1][j]; Fn = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Sn[i][j] * Msh.vField.vConv[1][i][j+1]; Pn = Fn / Dn;
            Ds = - gammas * Msh.Ss[i][j] / Msh.nd[1][j-1]; Fs = Mat.vMat[Msh.nMat[i][j]].rho * Msh.Ss[i][j] * Msh.vField.vConv[1][i][j-1]; Ps = Fs / Ds;


            // Coefficients B 
            sumConv = Mat.vMat[Msh.nMat[i][j]].rho * (Msh.vField.vConv[0][i+1][j] * Msh.vPhi[i+1][j] * Msh.Se[i][j] - Msh.vField.vConv[0][i-1][j] * Msh.vPhi[i-1][j] * Msh.Sw[i][j] + Msh.vField.vConv[1][i][j+1] * Msh.vPhi[i][j+1] * Msh.Sn[i][j] - Msh.vField.vConv[1][i][j-1] * Msh.vPhi[i][j-1] * Msh.Ss[i][j] + (Msh.vField.vConv[0][i+1][j] * Msh.Se[i][j] - Msh.vField.vConv[0][i-1][j] * Msh.Sw[i][j] + Msh.vField.vConv[1][i][j+1] * Msh.Sn[i][j] - Msh.vField.vConv[1][i][j-1] * Msh.Ss[i][j]) * Msh.vPhi[i][j]);
            sumDiff = (gammaw * Msh.Sw[i][j] * Msh.vPhi[i-1][j] / Msh.nd[0][i-1] + gammae * Msh.Se[i][j] * Msh.vPhi[i+1][j] / Msh.nd[0][i] + gammas * Msh.Ss[i][j] * Msh.vPhi[i][j-1] / Msh.nd[1][j-1] + gamman * Msh.Sn[i][j] * Msh.vPhi[i][j+1] / Msh.nd[1][j] - (gammaw * Msh.Sw[i][j] / Msh.nd[0][i-1] + gammae * Msh.Se[i][j] / Msh.nd[0][i] + gammas * Msh.Ss[i][j] / Msh.nd[1][j-1] + gamman * Msh.Sn[i][j] / Msh.nd[1][j]) * Msh.vPhi[i][j]);
            Msh.bp[k] = Msh.sPhi[i][j] * Msh.Vp[i][j] + Mat.vMat[Msh.nMat[i][j]].rho * Msh.Vp[i][j] * Msh.vPhi[i][j] / dt + (1 - beta) * (sumConv + sumDiff);
            Msh.bp[k] = Msh.sPhi[i][j] * Msh.Vp[i][j] + Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp * Msh.Vp[i][j] * Msh.vPhi[i][j] / dt + (1 - beta) * (sumConv + sumDiff);

        }
    }

}
