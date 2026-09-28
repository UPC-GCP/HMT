#ifndef MESHFVM_H_
#define MESHFVM_H_

// Self-Imports
#include "o01_Material.h"
#include "o02_MeshDEV.h"
#include "o09_Parser.h"

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
template <size_t Dim> class MeshFVM : Mesh<Dim> {
private: 
    // Headers
    void calculateFaces(std::array<size_t, Dim> cNode, Json::Value refData, std::array<std::vector<double>, Dim>& nFaces); // Mesh refinement

public:
    // Headers
    void generateMeshSolver(MeshSolver<Dim>& Msh, Json::Value qNode, Json::Value sections, Json::Value refinement, Json::Value obstacles); // Generate MeshSolver
    void addBoundariesSolver(MeshSolver<Dim>& Msh, Material Mat, Parser& Prs, Json::Value boundaries, double dInit, std::string sInit); // Boundaries MeshSolver
    void deriveMeshBase(MeshSolver<Dim> p, std::array<MeshBase<Dim>, Dim>& V); // Generate MeshBase
    // Was it worth it to pass the std::array<MeshBase<Dim>, Dim>& V ? Or should I just pass MeshBase<Dim>
    // Test with addBoundariesBase() -- The code doesn't know what index to access if I send them independently, would need to pass that as well -- in main() -> for (size_t i = 0; i < Dim; i++) { Msh.addBoundariesBase(i, V, Mat, Prs, boundaries, bInit, sInit) }
    void addBoundariesBase(size_t i, std::array<MeshBase<Dim>, Dim>& V, Material Mat, Parser& Prs, Json::Value boundaries, std::vector<double> dInit, std::vector<std::string> sInit); // Boundaries MeshBase
};

#endif
