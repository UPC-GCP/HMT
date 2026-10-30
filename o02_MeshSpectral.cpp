// Imports
#include <cmath>
#include <math.h>
#include <cstddef>

// Self-Imports
#include "o02_MeshSpectral.h"

void MeshSpectral::generateMeshBurgers(MeshBurgers& Msh, size_t N) {
    // Control
    Msh.totNodes = N+1; Msh.N[0] = N;

    // Resize
    Msh.uHat.resize(Msh.totNodes); Msh.ouHat.resize(Msh.totNodes); Msh.E.resize(Msh.totNodes); Msh.visc.resize(Msh.totNodes); Msh.viscEff.resize(Msh.totNodes);
}

void MeshSpectral::generateMeshBurgers(MeshLES& Msh, size_t N, double smagorinskyConstant, double slope) {
    // Base
    generateMeshBurgers(Msh, N);

    // Control
    Msh.Cs = smagorinskyConstant; Msh.m = slope; Msh.filterWidth = M_PI / N;

    // Resize
    Msh.viscStar.resize(Msh.totNodes);

    // Constants
    Msh.ls = Msh.Cs * Msh.filterWidth; Msh.viscInf = 0.31 * (5 - Msh.m) * std::sqrt(3 - Msh.m) * std::pow(Msh.Cs, -1.5) / (Msh.m + 1);
}
