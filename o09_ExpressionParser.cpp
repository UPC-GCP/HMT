// Imports
#include <vector>

#define _USE_MATH_DEFINES
#include <math.h>

#define exprtk_enable_all_features
#define exprtk_disable_string_capabilities
#include "exprtk.hpp"

// Self Imports
#include "o09_ExpressionParser.h"

ExpressionParser::ExpressionParser(){

    // Constants
    symbol_table.add_constant("pi", varPi);

    // Variables
    symbol_table.add_variable("t", varTime);
    symbol_table.add_variable("x", varX);
    symbol_table.add_variable("y", varY);

}

int ExpressionParser::registerExpression(std::string exprStr){

    // Parse Expression
    exprtk::expression<double> exprTemp;
    exprTemp.register_symbol_table(symbol_table);
    parser.compile(exprStr, exprTemp);

    // Store
    vExpr.push_back(exprTemp);

    return vExpr.size() - 1;

}

double ExpressionParser::evaluateTime(int i, double nVal){

    // Update
    varTime = nVal; return vExpr[i].value();

}

double ExpressionParser::evaluateCoordinates(int i, double xCoord, double yCoord){

	// Update
	varX = xCoord; varY = yCoord; return vExpr[i].value();

}
