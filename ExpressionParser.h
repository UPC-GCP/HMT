#ifndef EXPRESSIONPARSER_H_
#define EXPRESSIONPARSER_H_

#include <vector>

#define _USE_MATH_DEFINES
#include <math.h>

#define exprtk_enable_all_features
#define exprtk_disable_string_capabilities
#include "exprtk.hpp"

class ExpressionParser
{
private:

public:
    // Variables
    exprtk::symbol_table<double> symbol_table;
    exprtk::parser<double> parser;
    double varTime;
    const double varPi = M_PI;

    // Vectors
    std::vector<exprtk::expression<double>> vExpr;

    // Constructor
    ExpressionParser();

    // Functions
    int registerExpression(std::string exprStr);
    double evaluateTime(int i, double nVal);
};

#endif