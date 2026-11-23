#ifndef MESHFVM_H_
#define MESHFVM_H_

// Self-Imports
#include "o01_Material.h"
#include "o02_Mesh.h"
#include "o09_Parser.h"

template <size_t Dim> struct Matrix {
    double ap{}; // Node coefficient
    std::array<double, 2 * Dim> ak{}; // Neighbour coefficients
};

template <size_t Dim> struct Obstacle {
    std::array<size_t, Dim> i0{}, i1{}; // Indexes -> [nAxis]
};

template <size_t Dim> struct Boundary {
    size_t type{}, side{}, iExpr{}, iEq{}, iExprA{}, iEqA{}; double value{}, alpha{}; // General
    bool bUpdate = false, bA = false; std::string expression{}, expressionA{}; // Expression
    std::vector<size_t> iExprV{}, iEqV{}; std::vector<bool> bV{}; std::vector<double> valueV{}; std::vector<std::string> expressionV{}; // Velocity Field
    // Velocity field variables are there just for Velocity BCs
    // They shouldn't use Phi since it has dimensions
    // Phi should be a nD array as well if you wanted it to follow the same capacity as MeshSolver
    // I could just impose it as value for now and add it later
    std::array<size_t, Dim> i0{}, i1{}; // Indexes -> [nAxis]
    std::vector<double> Phi{}, oPhi{}; // Phi, oPhi -> [m] -> m = 1D flattened for 2D/3D
    std::vector<double> A{}, oA{}; // Alpha, oAlpha -> [m] -> Robin BC
};

template <size_t Dim> struct MeshBase : MeshSimplified<Dim> {
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
    void calculateFaces(const std::array<size_t, Dim>& cNode, const Json::Value& refData, std::array<std::vector<double>, Dim>& nFaces); // Mesh refinement

public:
    // Headers
    void generateMeshSolver(MeshSolver<Dim>& Msh, const Json::Value& qNode, const Json::Value& sections, const Json::Value& refinement, const Json::Value& obstacles);
    void deriveMeshBase(const MeshSolver<Dim>& p, std::array<MeshBase<Dim>, Dim>& V);
    void addBoundariesSolver(MeshSolver<Dim>& Msh, Material Mat, Parser& Prs, Json::Value boundaries, double dInit, std::string sInit);
    void addBoundariesBase(size_t i, std::array<MeshBase<Dim>, Dim>& V, const Material& Mat, Parser& Prs, const Json::Value& boundaries, std::vector<double> dInit, std::vector<std::string> sInit);
};

#endif
