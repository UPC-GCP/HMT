#ifndef MESH_H_
#define MESH_H_

#include <vector>
#include <json/json.h>

#include "o01_Material.h"
#include "o09_ExpressionParser.h"

struct Matrix{
    double ap=0, aw=0, ae=0, as=0, an=0;
};

struct Boundary{
    int type{}, side{}, iExpr{}, iEq{};
    std::vector<double> x0{}, x1{};
    double value{}, alpha{};
    bool bUpdate = false;
    std::string expression;
};

struct VelocityField{
    double Vw=0, Ve=0, Vs=0, Vn=0;
};

class Mesh
{
private:
    void calculateFaces(int cNode, int NSec, double x0, double x1, std::vector<double>& fVec);
public:
    // Variables
    double W{}, strength{}, centering{}, kStrength{}, delta{};
    int totNodes=1, algorithm{};
    
    // Vectors
    std::vector<int> N{}, vExpr{}; // dimension
    std::vector<double> bp{}, tempB{}; // k = x-axis * N[1] + y-axis 
    std::vector<Matrix> matA{}, tempA{}; // k = x-axis * N[1] + y-axis 
    std::vector<Boundary> boundaryConditions{};
    std::vector<std::vector<int>> nMat{}, nIgnore{}; // x-axis, y-axis
    std::vector<std::vector<double>> Faces{}, Nodes{}, DeltaX{}, dX{}; // dimension, position
    std::vector<std::vector<double>> sPhi{}, vPhi{}, Sw{}, Se{}, Ss{}, Sn{}, Vp{}; // x-axis, y-axis // Eventually get rid of this and store it in a single variable struct Areas{double Aw=0, Ae=0, As=0, An=0;}
    std::vector<std::vector<std::vector<VelocityField>>> vConv{}; // dimension, x-axis, y-axis

    // Constructor
    Mesh(int algo, double W = 1, double A = 0, double xC = 0.5, double kStr = 1, double delta = 0.001);

    // Functions
    void generateMesh(Material& Mat, Json::Value qNode, Json::Value sections, Json::Value refinement);
    void addBoundaryConditions(Json::Value boundaries, ExpressionParser& Prs);
    void addVelocityField(Json::Value nField, ExpressionParser& Prs);
};

#endif
