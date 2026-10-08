#ifndef MEDICSPECTRAL_H_
#define MEDICSPECTRAL_H_

// Self-Imports
#include "o06_Medic.h"
#include "o02_MeshSpectral.h"

// Class
class MedicSpectral : Medic<1> {
private:

public:
    // Constructor
    MedicSpectral(const MeshBurgers& Burg, const Probe<1>& Prb);

    // Headers
    void checkDiagnostics(const MeshBurgers& Burg);
};

#endif
