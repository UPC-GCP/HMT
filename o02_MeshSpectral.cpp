// Self-Imports
#include "o02_MeshSpectral.h"

void MeshSpectral::generateMeshBurgers(MeshBurgers Msh, size_t N) {
    // Control
    Msh.totNodes = N+1; 

    // Resize
    Msh.uHat.resize(N+1); Msh.ouHat.resize(N+1);
    Msh.E.resize(N+1); Msh.R.resize(N+1);
}

