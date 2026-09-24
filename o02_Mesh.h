#ifndef MESH3D_H_
#define MESH3D_H_

// Imports
#include <cstddef>
#include <json/json.h>
#include <optional>
#include <complex.h>

// Self-Imports
#include "o01_Material.h"
#include "o09_Parser.h"

// Types
namespace COMP { // Compass
    constexpr size_t W = 0, E = 1, S = 2, N = 3, B = 4, T = 5; // Direction indexes
    constexpr size_t X = 0, Y = 1, Z = 2; // Dimension indexes
}

template <size_t Dim> struct Matrix {
    double ap{}; // Node coefficient
    std::array<double, 2 * Dim> ak{}; // Neighbour coefficients
};

template <size_t Dim> struct Boundary {
    size_t type{}, side{}, iExpr{}, iEq{}, iExprA{}, iEqA{}; double value{}, alpha{}; // General
    bool bUpdate = false, bA = false; std::string expression{}, expressionA{}; // Expression
    std::array<size_t, Dim> i0{}, i1{}; // Indexes -> [nAxis]
    std::vector<double> Phi{}, oPhi{}; // Phi, oPhi -> [m] -> m = 1D flattened for 2D/3D
    std::vector<double> A{}, oA{}; // Alpha, oAlpha -> [m] -> Robin BC
};

template <size_t Dim> struct Obstacle {
    std::array<size_t, Dim> i0{}, i1{}; // Indexes -> [nAxis]
};

template <size_t Dim> struct MeshSimplified {
    size_t totNodes=1;
};

template <size_t Dim> struct MeshModal : MeshSimplified<Dim> { // Not now but change to other definition afterwards
/* struct MeshModal : MeshSimplified<1> { */
    double Re{}; std::vector<double> E{}, oE{};
    std::vector<std::complex<double>> uHat{}, ouHat{}, R{};
};

template <size_t Dim> struct MeshBase : MeshSimplified<Dim> {
    std::array<size_t, Dim> N{}; // Nodes -> [nAxis]
    std::array<std::vector<double>, Dim> Faces{}, Nodes{}, deltaX{}, dX{}; // Coordinates, distances -> [nAxis][index]
    std::vector<Matrix<Dim>> matA{}; std::vector<double> matB{}, oR{}; // Ax = b -> [l] -> l = i + Nx * (j + Ny * k)
    std::array<std::vector<double>, Dim> S{}; // Surfaces -> [nAxis][l]
    std::vector<double> Vp{}, Phi{}, oPhi{}; // Volume -> [l]
    std::vector<Boundary<Dim>> BC{}; // Boundary Conditions -> [index]
    std::vector<Obstacle<Dim>> Obs{}; // Obstacle -> [index]
};

template <size_t Dim> struct MeshSolver : MeshBase<Dim> {
    std::vector<Matrix<Dim>> tempA{}; std::vector<double> tempB{}; // Ax = b -> [l]
    std::vector<size_t> nMat{}; // Material -> [l]
    std::vector<double> sPhi{}; // Source term -> [l]
    std::vector<bool> bObs{}; // Obstacle -> [l]
};

// Class
template <size_t Dim> class Mesh {
private:
    // Headers
    void calculateFaces(std::array<size_t, Dim> cNode, Json::Value refData, std::array<std::vector<double>, Dim>& nFaces); // Mesh refinement

public:
    // Variables
    double epsFind=1e-8; // Config 

    // Headers
    void generateMeshSolver(MeshSolver<Dim>& Msh, Json::Value qNode, Json::Value sections, Json::Value refinement, Json::Value obstacles); // Generate MeshSolver
    void addBoundariesSolver(MeshSolver<Dim>& Msh, Material Mat, Parser& Prs, Json::Value boundaries, double dInit, std::string sInit); // Boundaries MeshSolver
    void deriveMeshBase(MeshSolver<Dim> p, std::array<MeshBase<Dim>, Dim>& V); // Generate MeshBase
    // Was it worth it to pass the std::array<MeshBase<Dim>, Dim>& V ? Or should I just pass MeshBase<Dim>
    // Test with addBoundariesBase() -- The code doesn't know what index to access if I send them independently, would need to pass that as well -- in main() -> for (size_t i = 0; i < Dim; i++) { Msh.addBoundariesBase(i, V, Mat, Prs, boundaries, bInit, sInit) }
    void addBoundariesBase(size_t i, std::array<MeshBase<Dim>, Dim>& V, Material Mat, Parser& Prs, Json::Value boundaries, std::vector<double> dInit, std::vector<std::string> sInit); // Boundaries MeshBase
};

// Functions
inline size_t calcIndex(size_t iX, size_t iY=0, size_t Ny=1, size_t iZ=0, size_t Nx=1) {return static_cast<size_t>(iY + Ny * (iX + Nx * iZ));};

template <size_t Dim, typename Func> void runLoopMesh(std::array<size_t, Dim> N, Func lamb, std::optional<std::array<size_t, Dim>> i0 = std::nullopt, std::optional<std::array<size_t, Dim>> i1 = std::nullopt) {
     
    // Control
    if (!i0) { for (size_t i = 0; i < Dim; i++) { (*i0)[i] = 0; (*i1)[i] = N[i]; } }
    size_t nLoop=1; for (size_t i = 0; i < Dim; i++) { nLoop *= ((*i1)[i] - (*i0)[i]); }

    // Pragma
    if constexpr (Dim == 1) { // 1D
        #pragma omp parallel for if (nLoop > 10000)
        for (size_t i = (*i0)[0]; i < (*i1)[1]; i++) {
            lamb({i, 0, 0}, 1, 1);
        }
    } else if constexpr (Dim == 2) { // 2D
        #pragma omp parallel for collapse(2) if (nLoop > 10000)
        for (size_t i = (*i0)[0]; i < (*i1)[1]; i++) {
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
