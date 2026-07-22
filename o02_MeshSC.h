#ifndef MESHSC_H_
#define MESHSC_H_

#include <json/value.h>
/* #include <vector> */
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
    std::vector<std::vector<bool>> bObs{}; // initialized
};

struct boundMain{
    int type{}, side{}, iExpr{}, iEq{}; double value{};
    std::vector<int> i0{}, i1{};
    std::vector<double> Phi{}, oPhi{}; // Phi, oPhi = 1D for their respective dimension
    bool bUpdate = false;
    std::string expression{};
};

struct boundVelocity{
    int type{}, side{};
    double uVal{}, vVal{};
    std::vector<int> i0{}, i1{}; // indexes
    std::vector<double> Phi{}, oPhi{}; // bound velocity modified to check the other 
    bool bUpdate = false;
    std::string expression{};
};

struct Obstacle{
    std::vector<int> i0{}, i1{}; // Indexes for blocked region
};

class Mesh
{
private:
    void calculateFaces(int cNode, int NSec, double x0, double x1, std::vector<double>& fVec);
    void calculateMeshGeometry(MeshBase& Msh, double valInit);

public:
    // Variables
    int algorithm{};
    double W{}, strength{}, centering{}, kStrength{}, delta{}, epsFind=1e-5; // Config
    MeshBase u{}, v{}; MeshSolver p{}, T{}; // Meshes

    // Vectors
    std::vector<int> vExpr{}; // dimension
    std::vector<boundMain> boundaryPressure{}, boundaryEnergy{};
    std::vector<boundVelocity> boundaryVelocity{};
    std::vector<Obstacle> obstacles{}; 

    // Constructor
    Mesh(int algo, double W = 1, double A = 0, double xC = 0.5, double kStr = 1, double delta = 0.001);

    // Functions
    void generateMesh(MeshSolver& Msh, double Phi0, Json::Value qNode, Json::Value sections, Json::Value refinement, Json::Value obs); // generate p, T
    void generateMeshVelocity(Material Mat, MeshSolver p, MeshBase& u, MeshBase& v); // generate u, v
    void addBoundariesPressure(Json::Value boundaries, Material Mat, ExpressionParser& Prs);
    void addBoundariesEnergy(Json::Value boundaries, Material Mat, ExpressionParser& Prs);
    void addBoundariesVelocity(Json::Value boundaries, Material Mat);
};

#endif

