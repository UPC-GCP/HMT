// Imports
#include <cmath>
#include <cstddef>

// Self-Imports
#include "o03_DiscretizerSpectral.h"
#include "o02_MeshSpectral.h"
#include "o04_SolverSpectral.h"

void DiscretizerSpectral::initializeBurgers(MeshBurgers& Msh, double Re) {
    // Initialize Values
    for (size_t k = 1; k < Msh.totNodes; k++) { Msh.uHat[k] = 1 / static_cast<double>(k); }
    for (size_t k = 0; k < Msh.totNodes; k++) { Msh.E[k] = std::real(Msh.uHat[k] * std::conj(Msh.uHat[k])); Msh.visc[k] = 1 / Re; Msh.viscEff = Msh.visc; }
}

void DiscretizerSpectral::initializeBurgers(MeshLES& Msh, double Re) {
    // Base
    initializeBurgers(static_cast<MeshBurgers&>(Msh), Re);

    // LES
    for (size_t k = 1; k < Msh.totNodes; k++) { Msh.viscStar[k] = 1 + 34.5 * exp(-3.03 * Msh.N[0] / k); } Msh.viscStar.front() = 1;
}

void DiscretizerSpectral::calculateConvectiveTerm(std::vector<std::complex<double>>& C, const std::vector<std::complex<double>>& u) {
    // Control
    int p{}, N = static_cast<int>(u.size()) - 1; std::complex<double> i(0, 1), up{}, uq{};
    
    // Triadic Interactions
    for (int k = 0; k <= N; k++) {
        for (int q = k - N; q <= N; q++) {
            // Control
            p = k - q; if (std::abs(p) > u.size()) { continue; }
            
            // Calculate
            up = (p >= 0) ? u[p] : std::conj(u[-p]); uq = (q >= 0) ? u[q] : std::conj(u[-q]);
            C[k] += up * (i * static_cast<double>(q)) * uq;
        }
    }
}

void DiscretizerSpectral::calculateEnergy(std::vector<double>& E, const std::vector<std::complex<double>>& uHat) {
    for (size_t k = 0; k < uHat.size(); k++) { E[k] = std::norm(uHat[k]); }
}

template <typename MeshClass> void DiscretizerSpectral::calculateVEff(MeshClass& Msh) {
    if constexpr (std::is_same_v<MeshClass, MeshBurgers>) { Msh.viscEff = Msh.visc; }
    else if constexpr (std::is_same_v<MeshClass, MeshLES>) { for (size_t k = 0; k < Msh.totNodes; k++) { Msh.viscEff[k] = Msh.visc[k] + Msh.viscInf * Msh.viscStar[k] * std::sqrt(Msh.E.back() / Msh.N[0]); }}
}

template <typename MeshClass> std::vector<std::complex<double>> DiscretizerSpectral::calculateRHS(const std::vector<std::complex<double>>& u, double Re, double t, solSpectral algorithm, MeshClass& Msh) {
    // Control
    std::vector<std::complex<double>> D(u.size()), C(u.size()), F(u.size()), R(u.size());

    // Dissipation
    calculateVEff(Msh);
    for (size_t k = 0; k < u.size(); k++) { D[k] = 2 * static_cast<double>(k) * static_cast<double>(k) * u[k] * Msh.viscEff[k]; }
    
    // Convection
    calculateConvectiveTerm(C, u);

    // Forcing Term
    F[1] = D[1] + C[1];

    // RHS
    for (size_t k = 0; k < u.size(); k++) { R[k] = - D[k] - C[k] + F[k]; } R[0] = 0;

    return R;
}

// Instances
template std::vector<std::complex<double>> DiscretizerSpectral::calculateRHS<MeshBurgers>(const std::vector<std::complex<double>>&, double, double, solSpectral, MeshBurgers&);
template std::vector<std::complex<double>> DiscretizerSpectral::calculateRHS<MeshLES>(const std::vector<std::complex<double>>&, double, double, solSpectral, MeshLES&);
template void DiscretizerSpectral::calculateVEff<MeshBurgers>(MeshBurgers&);
template void DiscretizerSpectral::calculateVEff<MeshLES>(MeshLES&);
