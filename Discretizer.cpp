// Imports
#include <iostream>
#include <vector>
#include <json/json.h>
#include <cmath>
#include <algorithm>

// Self-Imports
#include "Material.h"
#include "Mesh.h"
#include "Discretizer.h"
#include "ExpressionParser.h"


double Discretizer::calcHarmonicMean(double dPF, std::vector<double> lambda, std::vector<double> deltaX) {

    // Denominator
    double A = 0;
    for (int i = 0; i < lambda.size(); i++){
        A += (deltaX[i] / 2) / lambda[i];
    }
    
    return dPF / A;

}


Discretizer::Discretizer(std::string scheme, double endTime, double dt, double epsFind) {
    
    // Time Parameters
    this->scheme = scheme;
    this->endTime = endTime; this->dt = dt;
    this->epsFind = epsFind;

}


void Discretizer::setSchemeParameters(Material& Mat, Mesh& Msh){

    // Scheme Selection
    if (scheme == "explicit") {
        
        // Beta
        beta = 0;
        
        // Calculate Time-step
        std::vector<double> dtNew(Msh.totNodes, 0); double dtMin;
        for (int i = 0; i < Msh.totNodes; i++){
            dtNew[i] = 0.5 * pow(Msh.deltaX[i], 2) / Mat.vMat[Msh.xMat[i]].alpha;
        }

        // Update Time-step
        dtMin = *std::min_element(dtNew.begin()+1, dtNew.end()-1);
        if (dtMin < dt) {dt = dtMin;}
        
    } else if (scheme == "crank-nicolson") {

        // Beta
        beta = 0.5;

    } else if (scheme == "implicit") {

        // Beta
        beta = 1;

    } else {
        std::cerr << "Error: Invalid discretization scheme " << scheme << std::endl;
    }

}


void Discretizer::newSetBoundaryConditions(Material& Mat, Mesh& Msh, ExpressionParser& Prs, double t){

    // Boundary Conditions
    std::vector<int> Pos0, Pos1; Pos0.resize(Msh.N.size()); Pos1.resize(Msh.N.size()); double lamb; int j;
    for (Boundary bC : Msh.newBoundaryConditions){

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
                for (int i = Pos0[1]; i < Pos1[1]; i++){
                    // Ignore Corners
                    if (i == 0 || i == Msh.N[1]-1){continue;}

                    // Value
                    Msh.nT[Pos0[0]][i] = bC.value;

                    // Coefficients
                    j = Pos0[0] * Msh.N[1] + i;
                    Msh.matA[j].ap = 1;
                    Msh.bp[j] = bC.value;

                    // Control
                    if (bIgnore){Msh.nIgnore.push_back({Pos0[0], i});}
                }

            } else if (Pos0[1] == Pos1[1]){
                
                // yBoundary
                for (int i = Pos0[0]; i < Pos1[0]; i++){
                    // Ignore Corners
                    if (i == 0 || i == Msh.N[0]-1){continue;}
                    
                    // Value
                    Msh.nT[i][Pos0[1]] = bC.value;

                    // Coefficients
                    j = i * Msh.N[1] + Pos0[1];
                    Msh.matA[j].ap = 1;
                    Msh.bp[j] = bC.value;

                    // Control
                    if (bIgnore){Msh.nIgnore.push_back({i, Pos0[1]});}
                }

            }

        } else if (bC.type == 1){
            
            // double stabCoeff = 1;

            // Neumann
            if (Pos0[0] == Pos1[0]){
                
                // xBoundary
                if (bC.side == 0){

                    // West Boundary // ae, ap, bp
                    for (int i = Pos0[1]; i < Pos1[1]; i++){
                        // Ignore Corners
                        if (i == 0 || i == Msh.N[1]-1){continue;}

                        // Thermal Conductivity
                        lamb = Mat.vMat[Msh.nMat[Pos0[0]][i]].lambda;

                        // Value
                        Msh.nT[Pos0[0]][i] = bC.value * Msh.nd[0][Pos0[0]] / lamb + Msh.nT[Pos0[0]+1][i];

                        // Coefficients
                        j = Pos0[0] * Msh.N[1] + i;
                        Msh.matA[j].ap = 1;
                        Msh.bp[j] = Msh.nT[Pos0[0]][i];

                        // std::cout << "\nCheck values for bp:\n";

                        // std::cout << "Node: " << Pos0[0] << " " << i << "\n";
                        // std::cout << "Q_Neumann: " << bC.value << "\n";
                        // std::cout << "beta: " << beta << "\n";
                        // std::cout << "lambda: " << lamb << "\n";
                        // std::cout << "T_E: " << Msh.nT[Pos0[0]+1][i] << "\n";
                        // std::cout << "nd: " << Msh.nd[0][Pos0[0]] << "\n";

                        // std::cout << "bp: " << Msh.bp[j] << "\n";

                        // std::system("pause");

                    }

                } else if (bC.side == 1){

                    // East Boundary // aw, ap, bp
                    for (int i = Pos0[1]; i < Pos1[1]; i++){
                        // Ignore Corners
                        if (i == 0 || i == Msh.N[1]-1){continue;}
                        
                        // Thermal Conductivity
                        lamb = Mat.vMat[Msh.nMat[Pos0[0]][i]].lambda;

                        // Value
                        Msh.nT[Pos0[0]][i] = bC.value * Msh.nd[0][Pos0[0]-1] / lamb + Msh.nT[Pos0[0]-1][i];

                        // Coefficients
                        j = Pos0[0] * Msh.N[1] + i;
                        Msh.matA[j].ap = 1;
                        Msh.bp[j] = Msh.nT[Pos0[0]][i];

                    }

                } else {std::cerr << "Boundary side not specified correcly.\n";}

            } else if (Pos0[1] == Pos1[1]){

                // yBoundary
                if (bC.side == 0){

                    // South Boundary // an, ap, bp
                    for (int i = Pos0[0]; i < Pos1[0]; i++){
                        // Ignore Corners
                        if (i == 0 || i == Msh.N[0]-1){continue;}

                        // Thermal Conductivity
                        lamb = Mat.vMat[Msh.nMat[i][Pos0[1]]].lambda;

                        // // Coefficients
                        // j = i * Msh.N[1] + Pos0[1];
                        // Msh.matA[j].an = - beta * lamb / Msh.nd[1][Pos0[1]];
                        // Msh.matA[j].ap = - stabCoeff * Msh.matA[j].an;
                        // Msh.bp[j] = bC.value + (1 - beta) * lamb * (Msh.nT[i][Pos0[1]] - Msh.nT[i][Pos0[1]+1]) / Msh.nd[1][Pos0[1]];
                        
                        // OLD COEFFICIENTS
                        // Msh.matA[iPos].aw = - beta * lamb / Msh.dx[iPos-1];
                        // Msh.matA[iPos].ap = - Msh.matA[iPos].aw;
                        // Msh.bp[iPos] = bC[2] + (1 - beta) * (lamb * Msh.TNodes[iPos-1] / Msh.dx[iPos-1] - lamb * Msh.TNodes[iPos] / Msh.dx[iPos-1]);


                        // DIRICHLET TYPE IMPLEMENTATION
                        // Value
                        Msh.nT[i][Pos0[1]] = bC.value * Msh.nd[1][Pos0[1]] / lamb + Msh.nT[i][Pos0[1]+1];

                        // Coefficients
                        j = i * Msh.N[1] + Pos0[1];
                        Msh.matA[j].ap = 1;
                        Msh.bp[j] = Msh.nT[i][Pos0[1]];

                    }

                } else if (bC.side == 1){

                    // North Boundary // as, ap, bp
                    for (int i = Pos0[0]; i < Pos1[0]; i++){
                        // Ignore Corners
                        if (i == 0 || i == Msh.N[0]-1){continue;}
                        
                        // Thermal Conductivity
                        lamb = Mat.vMat[Msh.nMat[i][Pos0[1]]].lambda;

                        // std::cout << "North BC: " << i << " " << Pos0[1];

                        // // Coefficients
                        // j = i * Msh.N[1] + Pos0[1];
                        // Msh.matA[j].as = - beta * lamb / Msh.nd[1][Pos0[1]-1];
                        // Msh.matA[j].ap = - stabCoeff * Msh.matA[j].an;
                        // Msh.bp[j] = bC.value + (1 - beta) * lamb * (Msh.nT[i][Pos0[1]-1] - Msh.nT[i][Pos0[1]]) / Msh.nd[1][Pos0[1]-1];


                        // Value
                        Msh.nT[i][Pos0[1]] = bC.value * Msh.nd[1][Pos0[1]-1] / lamb + Msh.nT[i][Pos0[1]-1];

                        // Coefficients
                        j = i * Msh.N[1] + Pos0[1];
                        Msh.matA[j].ap = 1;
                        Msh.bp[j] = Msh.nT[i][Pos0[1]];
                    }

                } else {std::cerr << "Boundary side not specified correcly.\n";}

            }

        } else if (bC.type == 2){
            
            // Convection
            if (Pos0[0] == Pos1[0]){
                
                // xBoundary
                if (bC.side == 0){

                    // West Boundary // ae, ap, bp
                    for (int i = Pos0[1]; i < Pos1[1]; i++){
                        // Ignore Corners
                        if (i == 0 || i == Msh.N[1]-1){continue;}

                        // Thermal Conductivity
                        lamb = Mat.vMat[Msh.nMat[Pos0[0]][i]].lambda;

                        // std::cout << "Test Convection: \n";
                        // std::cout << "Tinf: " << bC.value << "\n";
                        // std::cout << "alpha: " << bC.alpha << "\n";
                        // std::cout << "lamb: " << lamb << "\n";
                        // std::cout << "Te: " << Msh.nT[Pos0[0]+1][i] << "\n";
                        // std::cout << "dx: " << Msh.nd[0][Pos0[0]] << "\n";

                        // Value
                        Msh.nT[Pos0[0]][i] = (bC.alpha * bC.value + lamb * Msh.nT[Pos0[0]+1][i] / Msh.nd[0][Pos0[0]]) / (lamb/Msh.nd[0][Pos0[0]] + bC.alpha);

                        // std::cout << "Tp: " << Msh.nT[Pos0[0]][i] << "\n";

                        // std:system("pause");

                        // Coefficients
                        j = Pos0[0] * Msh.N[1] + i;
                        Msh.matA[j].ap = 1;
                        Msh.bp[j] = Msh.nT[Pos0[0]][i];

                    }

                } else if (bC.side == 1){

                    // East Boundary // aw, ap, bp
                    for (int i = Pos0[1]; i < Pos1[1]; i++){
                        // Ignore Corners
                        if (i == 0 || i == Msh.N[1]-1){continue;}
                        
                        // Thermal Conductivity
                        lamb = Mat.vMat[Msh.nMat[Pos0[0]][i]].lambda;

                        // Value
                        Msh.nT[Pos0[0]][i] = (bC.alpha * bC.value + lamb * Msh.nT[Pos0[0]-1][i] / Msh.nd[0][Pos0[0]-1]) / (lamb/Msh.nd[0][Pos0[0]-1] + bC.alpha);

                        // Coefficients
                        j = Pos0[0] * Msh.N[1] + i;
                        Msh.matA[j].ap = 1;
                        Msh.bp[j] = Msh.nT[Pos0[0]][i];

                    }

                } else {std::cerr << "Boundary side not specified correcly.\n";}

            } else if (Pos0[1] == Pos1[1]){
                
                // yBoundary
                if (bC.side == 0){

                    // South Boundary // an, ap, bp
                    for (int i = Pos0[0]; i < Pos1[0]; i++){
                        // Ignore Corners
                        if (i == 0 || i == Msh.N[0]-1){continue;}

                        // Thermal Conductivity
                        lamb = Mat.vMat[Msh.nMat[i][Pos0[1]]].lambda;

                        // Value
                        Msh.nT[i][Pos0[1]] = (bC.alpha * bC.value + lamb * Msh.nT[i][Pos0[1]+1] / Msh.nd[1][Pos0[1]]) / (lamb/Msh.nd[1][Pos0[1]] + bC.alpha);

                        // Coefficients
                        j = i * Msh.N[1] + Pos0[1];
                        Msh.matA[j].ap = 1;
                        Msh.bp[j] = Msh.nT[i][Pos0[1]];

                    }

                } else if (bC.side == 1){

                    // North Boundary // as, ap, bp
                    for (int i = Pos0[0]; i < Pos1[0]; i++){
                        // Ignore Corners
                        if (i == 0 || i == Msh.N[0]-1){continue;}
                        
                        // Thermal Conductivity
                        lamb = Mat.vMat[Msh.nMat[i][Pos0[1]]].lambda;

                        // Value
                        Msh.nT[i][Pos0[1]] = (bC.alpha * bC.value + lamb * Msh.nT[i][Pos0[1]-1] / Msh.nd[1][Pos0[1]-1]) / (lamb/Msh.nd[1][Pos0[1]-1] + bC.alpha);

                        // Coefficients
                        j = i * Msh.N[1] + Pos0[1];
                        Msh.matA[j].ap = 1;
                        Msh.bp[j] = Msh.nT[i][Pos0[1]];
                        
                    }

                } else {std::cerr << "Boundary side not specified correcly.\n";}

            }

        }
        
    }

    // Corners
    Msh.nT[0][0] = 0.5 * (Msh.nT[1][0] + Msh.nT[0][1]);
    Msh.nT[0][Msh.N[1]-1] = 0.5 * (Msh.nT[1][Msh.N[1]-1] + Msh.nT[0][Msh.N[1]-2]);
    Msh.nT[Msh.N[0]-1][0] = 0.5 * (Msh.nT[Msh.N[0]-2][0] + Msh.nT[Msh.N[0]-1][1]);
    Msh.nT[Msh.N[0]-1][Msh.N[1]-1] = 0.5 * (Msh.nT[Msh.N[0]-2][Msh.N[1]-1] + Msh.nT[Msh.N[0]-1][Msh.N[1]-2]);

    // Control
    if (bIgnore){Msh.nIgnore.push_back({0, 0}); Msh.nIgnore.push_back({0, Msh.N[1]-1}); Msh.nIgnore.push_back({Msh.N[0]-1, 0}); Msh.nIgnore.push_back({Msh.N[0]-1, Msh.N[1]-1});}
    bIgnore = false;

}


void Discretizer::newSetCoefficients(Material& Mat, Mesh& Msh){
    
    // Control
    double lambw, lambe, lambs, lambn; int k;

    // Interior Nodes (Non-nD)
    for (size_t i = 1; i < Msh.N[0]-1; i++){
        for (size_t j = 1; j < Msh.N[1]-1; j++){

            // Harmonic Mean
            lambw = calcHarmonicMean(Msh.nd[0][i-1], {Mat.vMat[Msh.nMat[i][j]].lambda, Mat.vMat[Msh.nMat[i-1][j]].lambda}, {Msh.ndelta[0][i], Msh.ndelta[0][i-1]});
            lambe = calcHarmonicMean(Msh.nd[0][i], {Mat.vMat[Msh.nMat[i][j]].lambda, Mat.vMat[Msh.nMat[i+1][j]].lambda}, {Msh.ndelta[0][i], Msh.ndelta[0][i+1]});
            lambs = calcHarmonicMean(Msh.nd[1][j-1], {Mat.vMat[Msh.nMat[i][j]].lambda, Mat.vMat[Msh.nMat[i][j-1]].lambda}, {Msh.ndelta[1][j], Msh.ndelta[1][j-1]});
            lambn = calcHarmonicMean(Msh.nd[1][j], {Mat.vMat[Msh.nMat[i][j]].lambda, Mat.vMat[Msh.nMat[i][j+1]].lambda}, {Msh.ndelta[1][j], Msh.ndelta[1][j+1]});

            // Index
            k = i * Msh.N[1] + j;
            
            // Coefficients A
            Msh.matA[k].aw = - beta * lambw * Msh.nSw[i][j] / Msh.nd[0][i-1];
            Msh.matA[k].ae = - beta * lambe * Msh.nSe[i][j] / Msh.nd[0][i];
            Msh.matA[k].as = - beta * lambs * Msh.nSs[i][j] / Msh.nd[1][j-1];
            Msh.matA[k].an = - beta * lambn * Msh.nSn[i][j] / Msh.nd[1][j];
            Msh.matA[k].ap = Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp * Msh.nVp[i][j] / dt - Msh.matA[k].aw - Msh.matA[k].ae - Msh.matA[k].as - Msh.matA[k].an;

            // Coefficients B
            Msh.bp[k] = Msh.nQv[i][j] * Msh.nVp[i][j] + Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp * Msh.nT[i][j] * Msh.nVp[i][j] / dt + (1 - beta) * (lambw * Msh.nSw[i][j] * Msh.nT[i-1][j] / Msh.nd[0][i-1] + lambe * Msh.nSe[i][j] * Msh.nT[i+1][j] / Msh.nd[0][i] + lambs * Msh.nSs[i][j] * Msh.nT[i][j-1] / Msh.nd[1][j-1] + lambn * Msh.nSn[i][j] * Msh.nT[i][j+1] / Msh.nd[1][j] - (lambw * Msh.nSw[i][j] / Msh.nd[0][i-1] + lambe * Msh.nSe[i][j] / Msh.nd[0][i] + lambs * Msh.nSs[i][j] / Msh.nd[1][j-1] + lambn * Msh.nSn[i][j] / Msh.nd[1][j]) * Msh.nT[i][j]);


            // if (i == 1 && j == 1) {
            //     std::cout << "Node (1,1) Coefficients:\n";
            //     std::cout << "ap=" << Msh.matA[k].ap << ", aw=" << Msh.matA[k].aw << ", ae=" << Msh.matA[k].ae << ", as=" << Msh.matA[k].as << ", an=" << Msh.matA[k].an << "\n";
            //     std::cout << "bp=" << Msh.bp[k] << "\n";
            //     std::cout << "nQv=" << Msh.nQv[i][j] << ", nVp=" << Msh.nVp[i][j] << "\n";
            //     std::cout << "lambw=" << lambw << ", nSw=" << Msh.nSw[i][j] << ", nd[0][i-1]=" << Msh.nd[0][i-1] << "\n";

            //     std::system("pause");
            // }


        }
    }

}


void Discretizer::newSetRHS(Material& Mat, Mesh& Msh){

    // Control
    double lambw, lambe, lambs, lambn; int k;

    // Interior Nodes (Non-nD)
    for (size_t i = 1; i < Msh.N[0]-1; i++){
        for (size_t j = 1; j < Msh.N[1]-1; j++){
            
            // Harmonic Mean
            lambw = calcHarmonicMean(Msh.nd[0][i-1], {Mat.vMat[Msh.nMat[i][j]].lambda, Mat.vMat[Msh.nMat[i-1][j]].lambda}, {Msh.ndelta[0][i], Msh.ndelta[0][i-1]});
            lambe = calcHarmonicMean(Msh.nd[0][i], {Mat.vMat[Msh.nMat[i][j]].lambda, Mat.vMat[Msh.nMat[i+1][j]].lambda}, {Msh.ndelta[0][i], Msh.ndelta[0][i+1]});
            lambs = calcHarmonicMean(Msh.nd[1][j-1], {Mat.vMat[Msh.nMat[i][j]].lambda, Mat.vMat[Msh.nMat[i][j-1]].lambda}, {Msh.ndelta[1][j], Msh.ndelta[1][j-1]});
            lambn = calcHarmonicMean(Msh.nd[1][j], {Mat.vMat[Msh.nMat[i][j]].lambda, Mat.vMat[Msh.nMat[i][j+1]].lambda}, {Msh.ndelta[1][j], Msh.ndelta[1][j+1]});

            // Index
            k = i * Msh.N[1] + j;

            // Coefficients B
            Msh.bp[k] = Msh.nQv[i][j] * Msh.nVp[i][j] + Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp * Msh.nVp[i][j] * Msh.nT[i][j] / dt + (1 - beta) * (lambw * Msh.nSw[i][j] * Msh.nT[i-1][j] / Msh.nd[0][i-1] + lambe * Msh.nSe[i][j] * Msh.nT[i+1][j] / Msh.nd[0][i] + lambs * Msh.nSs[i][j] * Msh.nT[i][j-1] / Msh.nd[1][j-1] + lambn * Msh.nSn[i][j] * Msh.nT[i][j+1] / Msh.nd[1][j] - (lambw * Msh.nSw[i][j] / Msh.nd[0][i-1] + lambe * Msh.nSe[i][j] / Msh.nd[0][i] + lambs * Msh.nSs[i][j] / Msh.nd[1][j-1] + lambn * Msh.nSn[i][j] / Msh.nd[1][j]) * Msh.nT[i][j]);


            double qVVV{}, qTTT{}, qDDD{};

            qVVV = Msh.nQv[i][j] * Msh.nVp[i][j];
            qTTT = Mat.vMat[Msh.nMat[i][j]].rho * Mat.vMat[Msh.nMat[i][j]].cp * Msh.nVp[i][j] * Msh.nT[i][j] / dt;
            qDDD = (1 - beta) * (lambw * Msh.nSw[i][j] * Msh.nT[i-1][j] / Msh.nd[0][i-1] + lambe * Msh.nSe[i][j] * Msh.nT[i+1][j] / Msh.nd[0][i] + lambs * Msh.nSs[i][j] * Msh.nT[i][j-1] / Msh.nd[1][j-1] + lambn * Msh.nSn[i][j] * Msh.nT[i][j+1] / Msh.nd[1][j] - (lambw * Msh.nSw[i][j] / Msh.nd[0][i-1] + lambe * Msh.nSe[i][j] / Msh.nd[0][i] + lambs * Msh.nSs[i][j] / Msh.nd[1][j-1] + lambn * Msh.nSn[i][j] / Msh.nd[1][j]) * Msh.nT[i][j]);

            // std::cout << "Node (1,1) RHS components:\n";
            // std::cout << "qV*V = " << qVVV << "\n";
            // std::cout << "transient = " << qTTT << "\n";
            // std::cout << "explicit_diffusion = " << qDDD << "\n";
            // std::cout << "Total bp = " << qVVV + qTTT + qDDD << "\n";
            // std::cout << "Original bp = " << Msh.bp[k] << "\n";


            // std::system("pause");

        }
    }

}
