// Imports
#include <cmath>
#include <complex>
#include <cstddef>
#include <ios>
#include <json/json.h>
#include <iostream>
#include <complex.h>
#include <iomanip>
#include <vector>

// Self-Imports
#include "o02_Mesh.h"
#include "o06_Burgers.h"

    // Burgers 1D
    // Function organization for funcBurgers() in the corresponding model files
    // Mesh -> Creates the data structures used for spectral analysis (Resizes)
    // Discretizer -> Discretizes the spatial model into the coefficients (Initial conditions)
    // Solver -> Solves the Equation (RK)
    // Probe -> Store vector

    // DNS
    // Initialize vectors
        // N -> Grid points -- This determines the loops for k
        // k -> Wave number index -- Goes from -N to N -- k.size() = 2 * N + 1 -- Can be simplified to k.size() = N + 1 due to symmetry (u_(-k) = conj(u_k))
    // Initial conditions (uk, Fk)
    // Calculate Ck
    // Calculate Rk
    // RK4: Calculate uk n+1
    // Energy Eq -- Check conservation

void meshBurgers(MeshModal<1>& Msh, size_t N) { // Done

    // Control
    Msh.totNodes = N+1; 

    // Resize
    Msh.uHat.resize(N+1); Msh.ouHat.resize(N+1);
    Msh.E.resize(N+1); Msh.R.resize(N+1);

}

void discretizerBurgers(MeshModal<1>& Msh, double Re) { // Done

    // Initial Values
    Msh.Re = Re;
    for (size_t k = 1; k < Msh.totNodes; k++) { Msh.uHat[k] = 1 / static_cast<double>(k); }
    for (size_t k = 0; k < Msh.totNodes; k++) { Msh.E[k] = std::real(Msh.uHat[k] * std::conj(Msh.uHat[k])); }

}

void discretizerCalculateConvection(std::vector<std::complex<double>>& C, std::vector<std::complex<double>> const& u) { // Done

    int p{}, N = static_cast<int>(u.size()) - 1; std::complex<double> i(0, 1), up{}, uq{};
    for (int k = 0; k <= N; k++) {
        for (int q = k - N; q <= N; q++) {
            p = k - q; if (std::abs(p) > u.size()) { continue; } // Control
            
            up = (p >= 0) ? u[p] : std::conj(u[-p]);
            uq = (q >= 0) ? u[q] : std::conj(u[-q]);

            C[k] += up * (i * static_cast<double>(q)) * uq;
        }
    }

}

std::vector<std::complex<double>> discretizerCalculateRHS(std::vector<std::complex<double>> const& u, double Re, double t) { // Done

    // Control
    std::vector<std::complex<double>> D(u.size()), C(u.size()), F(u.size()), R(u.size());

    // Diffusion
    for (size_t k = 0; k < u.size(); k++) {
        D[k] = 2 * static_cast<double>(k) * static_cast<double>(k) * u[k] / Re;
    }
    
    // Convection
    discretizerCalculateConvection(C, u);

    // Forcing Term
    F[1] = D[1] + C[1];

    // RHS
    for (size_t k = 0; k < u.size(); k++) {
        R[k] = - D[k] - C[k] + F[k];
    } R[0] = 0;

    return R;

}

void solverRK4(std::vector<std::complex<double>>& u, double Re, double t, double dt) {

    // Control
    std::vector<std::complex<double>> uTemp{}, STemp{}, SFinal{}; uTemp.resize(u.size()); STemp.resize(u.size()); SFinal.resize(u.size());
    
    // S1
    STemp = discretizerCalculateRHS(u, Re, t);
    for (size_t k = 0; k < u.size(); k++) { SFinal[k] += dt * STemp[k] / static_cast<double>(6); }

    // S2
    for (size_t k = 0; k < u.size(); k++) { uTemp[k] = u[k] + 0.5 * dt * STemp[k]; }
    STemp = discretizerCalculateRHS(uTemp, Re, t + 0.5 * dt);
    for (size_t k = 0; k < u.size(); k++) { SFinal[k] += dt * 2 * STemp[k] / static_cast<double>(6); }

    // S3
    for (size_t k = 0; k < u.size(); k++) { uTemp[k] = u[k] + 0.5 * dt * STemp[k]; }
    STemp = discretizerCalculateRHS(uTemp, Re, t + 0.5 * dt);
    for (size_t k = 0; k < u.size(); k++) { SFinal[k] += dt * 2 * STemp[k] / static_cast<double>(6); }

    // S4
    for (size_t k = 0; k < u.size(); k++) { uTemp[k] = u[k] + dt * STemp[k]; }
    STemp = discretizerCalculateRHS(uTemp, Re, t + dt);
    for (size_t k = 0; k < u.size(); k++) { SFinal[k] += dt * STemp[k] / static_cast<double>(6); }

    // Update
    for (size_t k = 0; k < u.size(); k++) { u[k] = u[k] + SFinal[k]; }

}

double calcErr(std::vector<std::complex<double>> matA, std::vector<std::complex<double>> matB) { // Done

    // Control
    double rsNew{};
    std::complex<double> errVal{};
    std::vector<double> errVec{}; errVec.resize(matA.size(), 0);

    // Error
    for (size_t k = 0; k < matA.size(); k++) {
        errVal = abs(matA[k] - matB[k]);
        rsNew += std::real(errVal * errVal);
    }

    return std::sqrt(rsNew);

}

void discretizerCalculateEnergy(std::vector<double>& E, std::vector<std::complex<double>> uHat, double Re){
    for (size_t k = 0; k < uHat.size(); k++) { E[k] = std::norm(uHat[k]); }
}

void probeStoreData(std::vector<std::complex<double>> uHat, std::vector<double> E) {
    std::cout << "\n";
    std::cout << "uHat: "; for (std::complex<double> val : uHat) { std::cout << val << " "; } std::cout << "\n";
    std::cout << "E: "; for (double val : E) { std::cout << val << " "; } std::cout << "\n";
}

void medicRunDiagnostics(std::vector<std::complex<double>> u) {

    double sumVal = 0;
    std::vector<std::complex<double>> C(u.size()); discretizerCalculateConvection(C, u);

    for (size_t k = 1; k < u.size(); k++) {
        sumVal += 2 * std::real(std::conj(u[k]) * C[k]);
    }

    std::cout << "Check: " << sumVal << "\n";
}

void runBurgers(Json::Value dBurg) { // Everything will be done from this function but once finished with DNS/LES organize it within the main structure of the code.

    // Time
    double dt = dBurg["dt"].asDouble(); double endTime = dBurg["endTime"].asDouble(); // Should I send this to discretizer?
    double tolTemporal = 1e-3, tolNumeric = 1e-6; std::cout << std::scientific << std::setprecision(3);

    // Mesh
    MeshModal<1> Msh; std::cout << "Initializing mesh.\n";
    meshBurgers(Msh, dBurg["N"].asInt());
    
    // Discretizer
    discretizerBurgers(Msh, dBurg["Re"].asDouble()); std::cout << "Initial conditions set.\n";

    // Solver Loop
    std::cout << "Processing ...\n";
    probeStoreData(Msh.uHat, Msh.E);

    double t{}; size_t iMax = endTime / dt;
    for (size_t i = 0; i < iMax; i++) {

        // Control
        t = i * dt; Msh.ouHat = Msh.uHat; Msh.oE = Msh.E;

        // RK Solver
        solverRK4(Msh.uHat, Msh.Re, t, dt); // Updates Msh.u

        // Energy Balance
        discretizerCalculateEnergy(Msh.E, Msh.uHat, Msh.Re);

        // Write Data
        /* probeStoreData(Msh.uHat, Msh.E); */
        std::cout << "\r" << double(100 * static_cast<double>(i) / iMax) << " %";

        // Diagnostic
        medicRunDiagnostics(Msh.uHat);

        // Convergence
        if (calcErr(Msh.uHat, Msh.ouHat) / dt < tolTemporal) { std::cout << "\nSteady-state achieved @ t = " << std::setprecision(2) << t << " seconds."; break; }

    } std::cout << "\n";

    // Steady-State print
    std::cout << std::scientific << std::setprecision(3);
    probeStoreData(Msh.uHat, Msh.E);

    // End
    std::cout << "Files saved to: \n";

}
