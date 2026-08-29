#ifndef MESH3D_H_
#define MESH3D_H_

#include <cstddef>
#include <json/value.h>
#include <json/json.h>

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

inline size_t calcIndex(size_t iX, size_t iY=0, size_t Ny=1, size_t iZ=0, size_t Nx=1) {return static_cast<size_t>(iY + Ny * (iX + Nx * iZ));}; // General utility

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
    void addBoundariesPressure(MeshSolver<Dim>& Msh, Material Mat, Parser& Prs, Json::Value boundaries); // Boundaries (p)
    void addBoundariesTemperature(MeshSolver<Dim>& Msh, Material Mat, Parser& Prs, Json::Value boundaries); // Boundaries (T) 
    /* void addBoundariesPhi(MeshSolver<Dim>& Msh, Material Mat, Parser& Prs, Json::Value boundaries); // Boundaries (Phi) */

    /* template <size_t Dim> void generateMeshVelocity(); // Generate MeshBase */
    /* template <size_t Dim> void generateMeshVelocity(Material Mat, MeshSolver<Dim> p, MeshBase<Dim>& u, MeshBase<Dim>& v); // generate u, v, w */
    /* void addBoundariesVelocity(Json::Value boundaries, Material Mat); // Boundaries (u, v, w) */

    
};

///// Deprecated (Delete later, keep for now in case it serves as reference)
/* template <size_t Dim> struct boundVelocity{ */
/*     int type{}, side{}; double uVal{}, vVal{}, wVal{}; */
/*     bool bUpdate = false; std::string expression{}; */
/*     std::array<int, Dim> i0{}, i1{}; // Indexes -> [nAxis] */
/*     std::vector<double> Phi{}, oPhi{}; // Phi, oPhi -> [m] -> m = 2D flattened for the respective dimension */
/* }; */

#endif
