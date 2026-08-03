#ifndef MESH3D_H_
#define MESH3D_H_

#include <cstddef>
#include <json/value.h>
#include <json/json.h>

#include "o01_Material.h"
#include "o09_ExpressionParser.h"

struct Matrix{
    double ap=0, aw=0, ae=0, as=0, an=0, ab=0, at=0; // Neighbour coefficients
};

struct MeshBase{
    int totNodes=1; std::vector<size_t> N{}; // Nodes -> [nAxis]
    std::vector<std::vector<double>> Faces{}, Nodes{}, deltaX{}, dX{}; // Coordinates, distances -> [nAxis][index]
    std::vector<Matrix> matA{}; std::vector<double> matB{}, oR{}; // Ax = b -> [l] -> l = i + N[1] * (j + N[2] * k)

    std::vector<double> Sw{}, Se{}, Ss{}, Sn{}, Sb{}, St{}, Vp{}; // Geometry -> [l] -> l = i + N[1] * (j + N[2] * k)
    std::vector<double> Phi{}, oPhi{}; // Phi -> [l] -> l = i + N[1] * (j + N[2] * k)

    /* std::vector<std::vector<std::vector<double>>> Sw{}, Se{}, Ss{}, Sn{}, Vp{}; // Geometry -> [x][y][z] */ 
    /* std::vector<std::vector<std::vector<double>>> Phi{}, oPhi{}; // Phi -> [x][y][z] */
};

struct MeshSolver : MeshBase{
    std::vector<std::vector<std::vector<int>>> nMat{}; // Material -> [x][y][z]
    std::vector<std::vector<std::vector<double>>> sPhi{}; // Source-term -> [x][y][z]
    std::vector<Matrix> tempA{}; // A -> [l] -> l = i + N[1] * j + N[2] * k
    std::vector<double> tempB{}; // b -> [l] -> l = i * N[1] * j + N[2] * k
    std::vector<std::vector<std::vector<bool>>> bObs{}; // Obstacles -> [x][y][z] 
};

struct boundMain{
    int type{}, side{}, iExpr{}, iEq{}; double value{};
    bool bUpdate = false; std::string expression{};
    std::vector<int> i0{}, i1{}; // Indexes -> [nAxis]
    std::vector<std::vector<double>> Phi{}, oPhi{}; // Phi, oPhi -> [a][b] -> 2D for their respective dimension
};

struct boundVelocity{
    int type{}, side{}; double uVal{}, vVal{}, wVal{};
    bool bUpdate = false; std::string expression{};
    std::vector<int> i0{}, i1{}; // Indexes -> [nAxis]
    std::vector<std::vector<double>> Phi{}, oPhi{}; // Phi, oPhi -> [a][b] -> 2D for their respective dimension
};

struct Obstacle{
    std::vector<int> i0{}, i1{}; // Indexes -> [nAxis]
};

class Mesh
{
private:
    // Functions
    /* void calculateFaces(int cNode, int NSec, double x0, double x1, std::vector<double>& fVec); */
    void calculateFaces(std::vector<size_t> cNode, Json::Value refData, std::vector<std::vector<double>>& nFaces); // Refines each region with their own algorithm
    // Still need to pass Msh.Faces so it knows which mesh to work on.


    // Keep for now, haven't checked
    void calculateMeshGeometry(MeshBase& Msh, double valInit);

public:
    // Variables
    /* int algorithm{}; double A{}, xC{}, kappa{}, delta{}, epsFind=1e-5; // Config */
    double epsFind=1e-5; // Config 
    MeshBase u{}, v{}, w{}; MeshSolver p{}, T{}; // Meshes

    // Vectors
    std::vector<int> vExpr{}; // Velocity Expression -> [nAxis] 
    std::vector<boundMain> boundaryPressure{}, boundaryEnergy{}; // Boundaries (Pressure, Temperature)
    std::vector<boundVelocity> boundaryVelocity{}; // Boundaries (Velocity Components)
    std::vector<Obstacle> obstacles{}; // Obstacles

    // Constructor
    Mesh();
    /* Mesh(int algorithm, double A = 0, double xC = 0.5, double kappa = 1, double delta = 0.001); // Initialize Mesh */
    // Don't want to define this as a general parameter anymore
    // Define in "refinement" and separate it for each region instead

    // Functions
    void generateMesh(MeshSolver& Msh, double Phi0, Json::Value qNode, Json::Value sections, Json::Value refinement, Json::Value obs); // generate p, T
    void generateMeshVelocity(Material Mat, MeshSolver p, MeshBase& u, MeshBase& v); // generate u, v, w
    void addBoundariesPressure(Json::Value boundaries, Material Mat, ExpressionParser& Prs); // Boundaries (p)
    void addBoundariesEnergy(Json::Value boundaries, Material Mat, ExpressionParser& Prs); // Boundaries (T) 
    void addBoundariesVelocity(Json::Value boundaries, Material Mat); // Boundaries (u, v, w)
};

#endif
