#ifndef MESH_H_
#define MESH_H_

#include <vector>
#include <json/json.h>

#include "Material.h"
#include "ExpressionParser.h"

// PENDING CLEAN VARIABLES AND FUNCTIONS

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
    double W{}, H{}, strength{}, centering{}, kStrength{}, delta{};

    // Vectors
    std::vector<int> ignoreBC{}, xMat{};
    std::vector<std::vector<double>> boundaryConditions{};
    std::vector<std::string> boundaryExpr{};
    std::vector<double> xFaces{}, xNodes{}, TNodes{}, Sw{}, Se{}, dx{}, deltaX{}, Vp{}, bp{}, qV{};
    std::vector<Matrix> matA{};


    // New Variables
    int totNodes=1, algorithm{};

    // New Vectors
    std::vector<Boundary> newBoundaryConditions{};
    std::vector<int> N{};
    std::vector<std::vector<double>> Faces{}, Nodes{}, ndelta{}, nd{};
    std::vector<std::vector<int>> nMat{}, nIgnore{};
    std::vector<std::vector<double>> nQv{}, nT{}, nSw{}, nSe{}, nSs{}, nSn{}, nVp{};


    // Constructor
    Mesh(int algo, double W = 1, double H = 1, double A = 0, double xC = 0.5, double kStr = 1, double delta = 0.001);

    // Functions
    bool isFormula(std::string value);
    std::vector<double> splitString(std::string strSplit, char delimiter);
    void newCalculateFaces(int cNode, int NSec, double x0, double x1, std::vector<double>& fVec); // Pueden ser locales

    void newAddBoundaryConditions(Json::Value boundaries, ExpressionParser& Prs);
    void newGenerateMesh(Material& Mat, Json::Value qNode, Json::Value sections, Json::Value refinement);
    
};

#endif