#ifndef MESH_H_
#define MESH_H_

#include <vector>
#include <json/json.h>

#include "Material.h"
#include "ExpressionParser.h"

struct Matrix{
    double ap=0, aw=0, ae=0, as=0, an=0;
};

struct Boundary{
    int type{}, side{}, iExpr{};
    std::vector<double> x0{}, x1{};
    double value{}, alpha{};
    bool bUpdate = false;
    std::string expression;
};

class Mesh
{
private:

public:
    // Variables
    double W{}, strength{}, centering{}, kStrength{}, delta{};
    int totNodes=1, algorithm{};

    // Vectors
    std::vector<double> bp{};
    std::vector<Matrix> matA{};
    std::vector<Boundary> newBoundaryConditions{};
    std::vector<int> N{};
    std::vector<std::vector<double>> Faces{}, Nodes{}, ndelta{}, nd{}; // dimensions, values
    std::vector<std::vector<int>> nMat{}, nIgnore{}; // x-axis, y-axis
    std::vector<std::vector<double>> nQv{}, nT{}, nSw{}, nSe{}, nSs{}, nSn{}, nVp{}; // x-axis, y-axis

    // Constructor
    Mesh(int algo, double W = 1, double A = 0, double xC = 0.5, double kStr = 1, double delta = 0.001);

    // Functions
    bool isFormula(std::string value);
    std::vector<double> splitString(std::string strSplit, char delimiter);
    void newCalculateFaces(int cNode, int NSec, double x0, double x1, std::vector<double>& fVec);
    void newGenerateMesh(Material& Mat, Json::Value qNode, Json::Value sections, Json::Value refinement);
    void newAddBoundaryConditions(Json::Value boundaries, ExpressionParser& Prs);
};

#endif