#ifndef MESH3D_H_
#define MESH3D_H_

#include <array>
#include <cstddef>
#include <json/value.h>
#include <json/json.h>

/* #include "exprtk.hpp" */
#include "o01_Material.h"
#include "o09_ExpressionParser.h"

namespace Compass{
    constexpr size_t W = 0, E = 1, S = 2, N = 3, B = 4, T = 5; // Direction indexes
    constexpr size_t X = 0, Y = 0, Z = 0; // Dimension indexes
}

template <size_t Dim> struct Matrix{
    double ap{}; // Node coefficient
    std::array<double, 2 * Dim> ak{}; // Neighbour coefficients
};

template <size_t Dim> struct MeshBase{
    int totNodes=1; std::vector<size_t> N{}; // Nodes -> [nAxis]
    std::array<std::vector<double>, Dim> Faces{}, Nodes{}, deltaX{}, dX{}; // Coordinates, distances -> [nAxis][index]
    std::vector<Matrix<Dim>> matA{}; std::vector<double> matB{}, oR{}; // Ax = b -> [l] -> l = i + Nx * (j + Ny * k)
    std::array<std::vector<double>, Dim> S{}; // Surfaces -> [nAxis][l]
    std::vector<double> Vp{}, Phi{}, oPhi{}; // Volume, Phi, oPhi -> [l]
};

template <size_t Dim> struct MeshSolver : MeshBase<Dim>{
    std::vector<Matrix<Dim>> tempA{}; std::vector<double> tempB{}; // Ax = b -> [l]
    std::vector<size_t> nMat{}; // Material -> [l]
    std::vector<double> sPhi{}; // Source term -> [l]
    std::vector<bool> bObs{}; // Obstacle -> [l]
};

template <size_t Dim> struct boundMain{
    int type{}, side{}, iExpr{}, iEq{}; double value{};
    bool bUpdate = false; std::string expression{};
    std::array<int, Dim> i0{}, i1{}; // Indexes -> [nAxis]
    std::vector<double> Phi{}, oPhi{}; // Phi, oPhi -> [m] -> m = 2D flattened for the respective dimension
};

template <size_t Dim> struct boundVelocity{
    int type{}, side{}; double uVal{}, vVal{}, wVal{};
    bool bUpdate = false; std::string expression{};
    std::array<int, Dim> i0{}, i1{}; // Indexes -> [nAxis]
    std::vector<double> Phi{}, oPhi{}; // Phi, oPhi -> [m] -> m = 2D flattened for the respective dimension
};

template <size_t Dim> struct Obstacle{
    std::array<int, Dim> i0{}, i1{}; // Indexes -> [nAxis]
};

class Mesh
{
private:
    // Functions
    void calculateFaces(std::vector<size_t> cNode, Json::Value refData, std::vector<std::vector<double>>& nFaces); // Refines each region with their own algorithm
    template <size_t Dim> void calculateMeshGeometry(MeshBase<Dim>& Msh, double valInit);

public:
    // Variables
    double epsFind=1e-5; // Config 
    /* MeshBase u{}, v{}, w{}; MeshSolver p{}, T{}; // Meshes */
    // Meshes should be created directly in o00 not stored here
    // u, v, w need to be defined in a single vector that stores each MeshBase and loops over them

    // Vectors
    /* std::vector<int> vExpr{}; // Velocity Expression -> [nAxis] */ 
    /* std::vector<boundMain> boundaryPressure{}, boundaryEnergy{}; // Boundaries (Pressure, Temperature) */
    /* std::vector<boundVelocity> boundaryVelocity{}; // Boundaries (Velocity Components) */
    /* std::vector<Obstacle> obstacles{}; // Obstacles */
    // This should also be created in o00 and populated by o02

    // Constructor
    Mesh();
    /* Mesh(int algorithm, double A = 0, double xC = 0.5, double kappa = 1, double delta = 0.001); // Initialize Mesh */
    // Don't want to define this as a general parameter anymore
    // Define in "refinement" and separate it for each region instead

    // Functions
    template <size_t Dim> void generateMesh(MeshSolver<Dim>& Msh, double Phi0, Json::Value qNode, Json::Value sections, Json::Value refinement, Json::Value obs); // generate p, T
    template <size_t Dim> void generateMeshVelocity(Material Mat, MeshSolver<Dim> p, MeshBase<Dim>& u, MeshBase<Dim>& v); // generate u, v, w
    void addBoundariesPressure(Json::Value boundaries, Material Mat, ExpressionParser& Prs); // Boundaries (p)
    void addBoundariesEnergy(Json::Value boundaries, Material Mat, ExpressionParser& Prs); // Boundaries (T) 
    void addBoundariesVelocity(Json::Value boundaries, Material Mat); // Boundaries (u, v, w)
};

#endif
