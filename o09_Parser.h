#ifndef PARSER_H_ 
#define PARSER_H_

#include <cmath>
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
    exprtk::parser<double> parser;
    exprtk::symbol_table<double> symbol_table;
    double varTime{}, varX{}, varY{}, varZ{};
    // Alternative: template <size_t Dim> std::array<double, Dim> varK{};

    // Vectors
    std::vector<std::string> sExpr;
    std::vector<exprtk::expression<double>> vExpr;

    // Constructor
    Parser();

    // Functions
    size_t registerExpression(std::string exprStr);
    double evaluateTime(int i, double nVal){
        varTime = nVal; return vExpr[i].value();
    };
    double evaluateCoordinates(int i, double xCoord, double yCoord=NAN, double zCoord=NAN){
        varX = xCoord; if (!std::isnan(yCoord)) {varY = yCoord;} if (!std::isnan(zCoord)) {varZ = zCoord;} return vExpr[i].value();
    };
};

#endif
