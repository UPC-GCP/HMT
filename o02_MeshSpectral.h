#ifndef MESHSpectral_H_
#define MESHSpectral_H_

// Self-Imports
#include "o02_MeshDEV.h"
#include <cstddef>

// Types
template <size_t Dim> struct MeshBurgers : MeshSimplified<Dim> {
    double Re{}; std::vector<double> E{}, oE{};
    std::vector<std::complex<double>> uHat{}, ouHat{}, R{};
};

// Class
class MeshSpectral : Mesh<1> {

};

#endif
