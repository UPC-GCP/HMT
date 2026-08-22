#ifndef MESH3D_H_
#define MESH3D_H_

#include <cstddef>
#include <json/value.h>
#include <json/json.h>

#include "o01_Material.h"
/* #include "o09_ExpressionParser.h" */

namespace Compass {
    constexpr size_t W = 0, E = 1, S = 2, N = 3, B = 4, T = 5; // Direction indexes
    constexpr size_t X = 0, Y = 1, Z = 2; // Dimension indexes
}

template <size_t Dim> struct Matrix {
    double ap{}; // Node coefficient
    std::array<double, 2 * Dim> ak{}; // Neighbour coefficients
};

template <size_t Dim> struct Boundary {
    size_t type{}, side{}, iExpr{}, iEq{}; double value{};
    bool bUpdate = false; std::string expression{};
    std::array<size_t, Dim> i0{}, i1{}; // Indexes -> [nAxis]
    std::vector<double> Phi{}, oPhi{}; // Phi, oPhi -> [m] -> m = 2D flattened for the respective dimension
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
    std::vector<Boundary<Dim>> BC{};
};

template <size_t Dim> struct MeshSolver : MeshBase<Dim> {
    std::vector<Matrix<Dim>> tempA{}; std::vector<double> tempB{}; // Ax = b -> [l]
    std::vector<size_t> nMat{}; // Material -> [l]
    std::vector<double> sPhi{}; // Source term -> [l]
    std::vector<bool> bObs{}; // Obstacle -> [l]
    std::vector<Obstacle<Dim>> Obs{}; // Obstacle -> [index] // Not entirely convinced about having this here -- May change in the future
};

class Mesh {
private:
    // Functions
    void calculateFaces(std::vector<size_t> cNode, Json::Value refData, std::vector<std::vector<double>>& nFaces); // Mesh refinement
    template <size_t Dim> void calculateMeshGeometry(MeshBase<Dim>& Msh, double valInit); // What was this for?

public:
    // Variables
    double epsFind=1e-5; // Config 

    // Constructor
    Mesh();

    // Functions
    template <size_t Dim> void generateMeshSolver(MeshSolver<Dim>& Msh, Material Mat, Json::Value qNode, Json::Value sections, Json::Value refinement, Json::Value obs); // Generate MeshSolver
    template <size_t Dim> void generateMeshVelocity(); // Generate MeshBase

    // PENDING
    /* template <size_t Dim> void generateMeshVelocity(Material Mat, MeshSolver<Dim> p, MeshBase<Dim>& u, MeshBase<Dim>& v); // generate u, v, w */
    /* void addBoundariesPressure(Json::Value boundaries, Material Mat, ExpressionParser& Prs); // Boundaries (p) */
    /* void addBoundariesEnergy(Json::Value boundaries, Material Mat, ExpressionParser& Prs); // Boundaries (T) */ 
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
