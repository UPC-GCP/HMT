// Imports
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

// Self-Imports
#include "o06_MedicSpectral.h"
#include "o01_Material.h"
#include "o02_MeshSpectral.h"
#include "o03_Discretizer.h"
#include "o03_DiscretizerSpectral.h"
#include "o04_Solver.h"
#include "o05_Probe.h"

MedicSpectral::MedicSpectral(const MeshBurgers& Burg, const Probe<1>& Prb, double timeStep) : Medic(Prb) {
    // Control
    for (std::vector<double>* vec : {&dE, &vErr}) { vec->resize(Burg.totNodes); } dt = timeStep;
}

void MedicSpectral::checkDiagnostics(const Material& Mat, const MeshBurgers& Burg) {
    // Transfer Spectrum
    std::vector<std::complex<double>> C(Burg.totNodes); DiscretizerSpectral::calculateConvectiveTerm(C, Burg.uHat);
    double sumVal = 0; for (size_t k = 1; k < Burg.totNodes; k++) { sumVal += 2 * std::real(std::conj(Burg.uHat[k]) * C[k]); }
    if (std::abs(sumVal) > epsTol) { throw std::runtime_error("Transfer spectrum: " + std::to_string(sumVal)); }

    // Energy Transport Equation
    std::vector<std::complex<double>> D(Burg.totNodes), F(Burg.totNodes);
    for (size_t k = 1; k < Burg.totNodes; k++) { D[k] = 2 * static_cast<double>(k) * static_cast<double>(k) * Burg.uHat[k] / Mat.Re; } DiscretizerSpectral::calculateConvectiveTerm(C, Burg.uHat); F[1] = D[1] + C[1];
    for (size_t k = 1; k < Burg.totNodes; k++) { dE[k] = 2 * k * k * Burg.oE[k] / Mat.Re - 2 * std::real(std::conj(Burg.uHat[k]) * C[k]) + 2 * std::real(std::conj(Burg.uHat[k]) * F[k]); vErr[k] = Burg.E[k] - (Burg.oE[k] + dt * dE[k]); }
    double eNorm = calcErr(vErr, vErr); if (eNorm > epsTol) { throw std::runtime_error("Energy transport equation: " + std::to_string(eNorm)); }
}
