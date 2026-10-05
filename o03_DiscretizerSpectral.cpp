// Self-Imports
#include "o03_DiscretizerSpectral.h"

void DiscretizerSpectral::initializeBurgers(MeshBurgers Msh) {
    // Initialize Values
    for (size_t k = 1; k < Msh.totNodes; k++) { Msh.uHat[k] = 1 / static_cast<double>(k); }
    for (size_t k = 0; k < Msh.totNodes; k++) { Msh.E[k] = std::real(Msh.uHat[k] * std::conj(Msh.uHat[k])); }
}

void calculateConvectiveTerm(std::vector<std::complex<double>>& C, std::vector<std::complex<double>> const& u) {
    // Control
    int p{}, N = static_cast<int>(u.size()) - 1; std::complex<double> i(0, 1), up{}, uq{};
    
    // Triadic Interactions
    for (int k = 0; k <= N; k++) {
        for (int q = k - N; q <= N; q++) {
            // Control
            p = k - q; if (std::abs(p) > u.size()) { continue; }
            
            // Calculate
            up = (p >= 0) ? u[p] : std::conj(u[-p]);
            uq = (q >= 0) ? u[q] : std::conj(u[-q]);
            C[k] += up * (i * static_cast<double>(q)) * uq;
        }
    }
}

std::vector<std::complex<double>> DiscretizerSpectral::calculateRHS(std::vector<std::complex<double>> const& u, double Re, double t) {
    // Control
    std::vector<std::complex<double>> D(u.size()), C(u.size()), F(u.size()), R(u.size());

    // Diffusion
    for (size_t k = 0; k < u.size(); k++) {
        D[k] = 2 * static_cast<double>(k) * static_cast<double>(k) * u[k] / Re;
    }
    
    // Convection
    calculateConvectiveTerm(C, u);

    // Forcing Term
    F[1] = D[1] + C[1];

    // RHS
    for (size_t k = 0; k < u.size(); k++) {
        R[k] = - D[k] - C[k] + F[k];
    } R[0] = 0;

    return R;
}
    
void DiscretizerSpectral::calculateEnergy(std::vector<double>& E, const std::vector<std::complex<double>>& uHat, double Re) {
    for (size_t k = 0; k < uHat.size(); k++) { E[k] = std::norm(uHat[k]); }
}
