#ifndef MESHSpectral_H_
#define MESHSpectral_H_

// Imports
#include <cstddef>

// Self-Imports
#include "o02_Mesh.h"

// Types
struct MeshBurgers : MeshSimplified<1> {
    std::vector<double> E{}, oE{};
    std::vector<std::complex<double>> uHat{}, ouHat{}, R{};
};

// Class
class MeshSpectral : Mesh<1> {
private:
    
public:
    // Headers
    void generateMeshBurgers(MeshBurgers Msh, size_t N);
};

#endif
