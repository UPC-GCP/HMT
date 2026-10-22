#ifndef PARSER_H_ 
#define PARSER_H_

#include <array>
#include <cmath>
#include <cstddef>
#define _USE_MATH_DEFINES
#include <math.h>

#define exprtk_enable_all_features
#define exprtk_disable_string_capabilities
#include "exprtk.hpp"

class Parser
{
private:

public:
    // Variables
    const double varPi = M_PI;
    double varTime{};
    exprtk::parser<double> parser;
    exprtk::symbol_table<double> symbol_table;
    /* double varTime{}, varX{}, varY{}, varZ{}; */

    // Vectors
    std::array<double, 3> varK{};
    std::vector<std::string> sExpr;
    std::vector<exprtk::expression<double>> vExpr;

    // Constructor
    Parser();

    // Functions
    size_t registerExpression(std::string exprStr);
    double evaluateTime(size_t i, double nVal) { varTime = nVal; return vExpr[i].value(); };
    double evaluateCoordinates(size_t i, std::array<double, 3> iX, size_t iRmv=3) { for (size_t k = 0, kk = 0; k < 3; k++) { if (k == iRmv) { continue; } varK[k] = iX[kk++]; } return vExpr[i].value(); };
};

#endif
