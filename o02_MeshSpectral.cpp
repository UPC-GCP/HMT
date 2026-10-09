// Self-Imports
#include "o02_MeshSpectral.h"

void MeshSpectral::generateMeshBurgers(MeshBurgers& Msh, size_t N) {
    // Control
    Msh.totNodes = N+1; Msh.N[0] = N;

    // Resize
    Msh.uHat.resize(Msh.totNodes); Msh.ouHat.resize(Msh.totNodes);
    Msh.E.resize(Msh.totNodes); Msh.R.resize(Msh.totNodes);
}

