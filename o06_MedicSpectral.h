#ifndef MEDICSPECTRAL_H_
#define MEDICSPECTRAL_H_

// Self-Imports
#include "o06_Medic.h"
#include "o01_Material.h"
#include "o02_MeshSpectral.h"

// Class
class MedicSpectral : Medic<1> {
private:

public:
    // Variables
    double dt{};
    std::vector<double> dE{}, vErr{};

    // Constructor
    MedicSpectral(const MeshBurgers& Burg, const Probe<1>& Prb, double timeStep);

    // Headers
    void checkDiagnostics(const Material& Mat, const MeshBurgers& Burg);
};

#endif
