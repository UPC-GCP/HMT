#ifndef MESHDEV_H_
#define MESHDEV_H_

// Imports
#include <cstddef>
#include <optional>
#include <complex.h>
#include <json/json.h>

namespace COMP { // Compass
    constexpr size_t W = 0, E = 1, S = 2, N = 3, B = 4, T = 5; // Direction indexes
    constexpr size_t X = 0, Y = 1, Z = 2; // Dimension indexes
}

// Types
template <size_t Dim> struct MeshSimplified {
    size_t totNodes=1; std::array<size_t, Dim> N{}; // Nodes -> [nAxis]
};

// Class
template <size_t Dim> class Mesh {
private:

public:

};

// Functions
inline size_t calcIndex(size_t iX, size_t iY=0, size_t Ny=1, size_t iZ=0, size_t Nx=1) { return static_cast<size_t>(iY + Ny * (iX + Nx * iZ)); };

template <size_t Dim, typename Func> void runLoopMesh(std::array<size_t, Dim> N, Func lamb, std::optional<std::array<size_t, Dim>> i0 = std::nullopt, std::optional<std::array<size_t, Dim>> i1 = std::nullopt) {
    // Control
    if (!i0) { for (size_t i = 0; i < Dim; i++) { (*i0)[i] = 0; (*i1)[i] = N[i]; } }
    size_t nLoop=1; for (size_t i = 0; i < Dim; i++) { nLoop *= ((*i1)[i] - (*i0)[i]); }

    // Pragma
    if constexpr (Dim == 1) { // 1D
        #pragma omp parallel for if (nLoop > 10000)
        for (size_t i = (*i0)[0]; i < (*i1)[0]; i++) {
            lamb({i, 0, 0}, 1, 1);
        }
    } else if constexpr (Dim == 2) { // 2D
        #pragma omp parallel for collapse(2) if (nLoop > 10000)
        for (size_t i = (*i0)[0]; i < (*i1)[0]; i++) {
            for (size_t j = (*i0)[1]; j < (*i1)[1]; j++) {
                lamb({i, j, 0}, N[1], 1);
            }
        }
    } else if constexpr (Dim == 3) { // 3D
        #pragma omp parallel for collapse(3) if (nLoop > 10000)
        for (size_t k = (*i0)[2]; k < (*i1)[2]; k++) {
            for (size_t i = (*i0)[0]; i < (*i1)[0]; i++) {
                for (size_t j = (*i0)[1]; j < (*i1)[1]; j++) {
                    lamb({i, j, k}, N[1], N[0]);
                }
            }
        }
    }
}

#endif
