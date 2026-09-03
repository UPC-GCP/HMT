#ifndef MESH3D_H_
#define MESH3D_H_

#include <cstddef>
#include <json/value.h>
#include <json/json.h>
/* #include <optional> */
#include <iostream>
#include <optional>
#include <omp.h>

#include "o01_Material.h"
#include "o09_Parser.h"

namespace Compass {
    constexpr size_t W = 0, E = 1, S = 2, N = 3, B = 4, T = 5; // Direction indexes
    constexpr size_t X = 0, Y = 1, Z = 2; // Dimension indexes
}

template <size_t Dim> struct Matrix {
    double ap{}; // Node coefficient
    std::array<double, 2 * Dim> ak{}; // Neighbour coefficients
};

template <size_t Dim> struct Boundary {
    size_t type{}, side{}, iExpr{}, iEq{}, iExprA{}, iEqA{}; double value{}, alpha{};
    bool bUpdate = false, bA = false; std::string expression{}, expressionA{};
    std::array<size_t, Dim> i0{}, i1{}; // Indexes -> [nAxis]
    std::vector<double> Phi{}, oPhi{}; // Phi, oPhi -> [m] -> m = 1D flattened for 2D/3D
    std::vector<double> A{}, oA{}; // Alpha, oAlpha -> [m] -> Robin BC
};

template <size_t Dim> struct Obstacle {
    std::array<size_t, Dim> i0{}, i1{}; // Indexes -> [nAxis]
};

template <size_t Dim> struct MeshBase {
    size_t totNodes=1; std::array<size_t, Dim> N{}; // Nodes -> [nAxis]
    std::array<std::vector<double>, Dim> Faces{}, Nodes{}, deltaX{}, dX{}; // Coordinates, distances -> [nAxis][index]
    std::vector<Matrix<Dim>> matA{}; std::vector<double> matB{}, oR{}; // Ax = b -> [l] -> l = i + Nx * (j + Ny * k)
    std::array<std::vector<double>, Dim> S{}; // Surfaces -> [nAxis][l]
    std::vector<double> Vp{}, Phi{}, oPhi{}; // Volume, Phi, oPhi -> [l]
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
    // Functions
    void calculateFaces(std::array<size_t, Dim> cNode, Json::Value refData, std::array<std::vector<double>, Dim>& nFaces); // Mesh refinement

public:
    // Variables
    double epsFind=1e-8; // Config 

    /* // Constructor */
    /* Mesh(); */

    // Functions
    void generateMeshSolver(MeshSolver<Dim>& Msh, Material Mat, Json::Value qNode, Json::Value sections, Json::Value refinement, Json::Value obstacles); // Generate MeshSolver
    void addBoundariesSolver(MeshSolver<Dim>& Msh, Material Mat, Parser& Prs, Json::Value boundaries, double bInit, std::string sInit); // Boundaries MeshSolver
    void generateMeshBase(); // Generate MeshBase
    void addBoundariesMase(); // Boundaries MeshBase

    /* void addBoundariesTemperature(MeshSolver<Dim>& Msh, Material Mat, Parser& Prs, Json::Value boundaries); // Boundaries (T) */ 
    /* void addBoundariesPhi(MeshSolver<Dim>& Msh, Material Mat, Parser& Prs, Json::Value boundaries); // Boundaries (Phi) */

    /* template <size_t Dim> void generateMeshVelocity(); // Generate MeshBase */
    /* template <size_t Dim> void generateMeshVelocity(Material Mat, MeshSolver<Dim> p, MeshBase<Dim>& u, MeshBase<Dim>& v); // generate u, v, w */
    /* void addBoundariesVelocity(Json::Value boundaries, Material Mat); // Boundaries (u, v, w) */

    
};


// General Utilities
inline size_t calcIndex(size_t iX, size_t iY=0, size_t Ny=1, size_t iZ=0, size_t Nx=1) {return static_cast<size_t>(iY + Ny * (iX + Nx * iZ));};

template <size_t Dim, typename Func> void runLoopMesh(Func lamb, std::optional<std::array<size_t, Dim>> i0 = std::nullopt, std::optional<std::array<size_t, Dim>> i1 = std::nullopt) {
     
    std::cout << "Loop\n";

    // PENDING -- DEBUG HERE
    // PONER DEFAULT Msh.N para casos sin i0, i1

    // Control
    if (i0) { std::cout << "i0: "; for (size_t val : *i0) {std::cout << val << " ";} std::cout << "\n"; }
    else {std::cout << "i0: null\n";} 

    size_t nLoop=1; for (size_t i = 0; i < Dim; i++) { nLoop *= ((*i1)[i] - (*i0)[i]); std::cout << "Axis " << i << ": " << (*i0)[i] << " " << (*i1)[i] << "\n";}

    std::cout << "Total: " << nLoop;

    throw std::invalid_argument("Check .json");
    return;

    if constexpr (Dim == 1) { // 1D
        #pragma omp parallel for if (nLoop > 10000)
        for (size_t i = (*i0)[0]; i < (*i1)[1]; i++) {
            lamb(i, 0, 0);
        }
    } else if constexpr (Dim == 2) { // 2D
        #pragma omp parallel for collapse(2) if (nLoop > 10000)
        for (size_t i = (*i0)[0]; i < (*i1)[1]; i++) {
            for (size_t j = (*i0)[1]; j < (*i1)[1]; j++) {
                lamb(i, j, 0);
            }
        }
    } else if constexpr (Dim == 3) { // 3D
        #pragma omp parallel for collapse(3) if (nLoop > 10000)
        for (size_t k = (*i0)[2]; k < (*i1)[2]; k++) {
            for (size_t i = (*i0)[0]; i < (*i1)[0]; i++) {
                for (size_t j = (*i0)[1]; j < (*i1)[1]; j++) {
                    lamb(i, j, k);
                }
            }
        }
    }

}

/* template <size_t Dim, typename Func> void loopMesh(Func lamb, std::optional<std::array<size_t, Dim>> i0 = std::nullopt, std::optional<std::array<size_t, Dim>> i1 = std::nullopt); */
/* template <size_t Dim, typename Func> void loopMesh(Func lamb, std::array<size_t, Dim> i0, std::array<size_t, Dim> i1) { */

///// Deprecated (Delete later, keep for now in case it serves as reference)

/* template <size_t Dim> struct boundVelocity{ */
/*     int type{}, side{}; double uVal{}, vVal{}, wVal{}; */
/*     bool bUpdate = false; std::string expression{}; */
/*     std::array<int, Dim> i0{}, i1{}; // Indexes -> [nAxis] */
/*     std::vector<double> Phi{}, oPhi{}; // Phi, oPhi -> [m] -> m = 2D flattened for the respective dimension */
/* }; */

#endif
