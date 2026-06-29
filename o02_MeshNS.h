#ifndef MESHNS_H_
#define MESHNS_H_

#include <vector>
#include <json/json.h>

#include "o01_Material.h"
#include "o09_ExpressionParser.h"

struct Matrix{
    double ap=0, aw=0, ae=0, as=0, an=0;
};

struct MeshBase{
    int totNodes=1;
    std::vector<int> N{};
    std::vector<std::vector<double>> Faces{}, Nodes{}, deltaX{}, dX{}; // Coordinates, distances
    std::vector<std::vector<double>> Sw{}, Se{}, Ss{}, Sn{}, Vp{}; // Geometry
    std::vector<std::vector<double>> Phi{}, oPhi{}; // Phi
    std::vector<Matrix> matA{}; // A
    std::vector<double> matB{}, oR{}; // b
};

struct MeshSolver : MeshBase{
    std::vector<std::vector<int>> nMat{}; // Material
    std::vector<std::vector<double>> sPhi{}; // Source-term
    std::vector<Matrix> tempA{}; // A
    std::vector<double> tempB{}; // b
};

struct boundMain{
    int type{}, side{}, iExpr{}, iEq{};
    std::vector<double> x0{}, x1{}, Phi{}, oPhi{}; // Phi, oPhi = 1D for their respective dimensionA
    double value{}, alpha{};
    bool bUpdate = false;
    std::string expression;
};

struct boundVelocity{
    int type{}, side{};
    double uVal{}, vVal{};
    std::vector<int> i0{}, i1{}; // each one stores the x-y coordinates of their position
};

class Mesh
{
private:
    void calculateFaces(int cNode, int NSec, double x0, double x1, std::vector<double>& fVec);
    void calculateMeshGeometry(MeshBase& Msh, double valInit);

public:
    // Variables
    int algorithm{};
    double W{}, strength{}, centering{}, kStrength{}, delta{}; // Config
    MeshBase u{}, v{}; MeshSolver p{}; // Meshes

    // Vectors
    std::vector<int> vExpr{}; // dimension
    std::vector<boundMain> boundaryConditions{};
    std::vector<boundVelocity> boundaryVelocity{};

    // Constructor
    Mesh(int algo, double W = 1, double A = 0, double xC = 0.5, double kStr = 1, double delta = 0.001);

    // Functions
    void generateMesh(MeshSolver& Msh, Material Mat, Json::Value qNode, Json::Value sections, Json::Value refinement); // generate p
    void generateMeshVelocity(Material Mat, MeshSolver p, MeshBase& u, MeshBase& v); // generate u, v
    void addBoundaryConditions(Json::Value boundaries, Material Mat, ExpressionParser& Prs);
    void addBoundariesVelocity(Json::Value boundaries, Material Mat);

    // How to organize Functions
    // generateMesh already creates MeshSolver. Should just use that one as is and let it return the solver. p as input;
    // generateVelocityField (NEW) - create Mesh for velocity fields (u, v) use p as an input since it already has all values needed.
    // addBoundaryConditions should now take all three meshes and update them all
};

#endif
