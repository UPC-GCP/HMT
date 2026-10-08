// Self-Imports
#include "o06_MedicSpectral.h"
#include "o02_MeshSpectral.h"
#include "o05_Probe.h"

MedicSpectral::MedicSpectral(const MeshBurgers& Burg, const Probe<1>& Prb) : Medic(Prb) {

    // How to configure this
    // Needs to check if sum(T) = 0
    // Also calculate Energy transport equation and validate against energy values
    // Probably should also take the mesh and resize according to that

}

void MedicSpectral::checkDiagnostics(const MeshBurgers& Burg) {

    // This will perform diagnostics 
    // Needs to read the mesh since it compares against it

    // Will not save anything in this case but will keep it as warnings in case it exceeds thresholds

}
