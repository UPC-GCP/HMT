// Imports
#include <cstddef>
#include <vector>

#define _USE_MATH_DEFINES
#include <math.h>

#define exprtk_enable_all_features
#define exprtk_disable_string_capabilities
#include "exprtk.hpp"

// Self Imports
#include "o09_Parser.h"

Parser::Parser(){

    // Constants
    symbol_table.add_constant("pi", varPi);

    // Variables
    symbol_table.add_variable("t", varTime);
    symbol_table.add_variable("x", varX);
    symbol_table.add_variable("y", varY);
    symbol_table.add_variable("z", varZ);

}


size_t Parser::registerExpression(std::string exprStr){

    // Control
    size_t iPos = std::find(sExpr.begin(), sExpr.end(), exprStr) - sExpr.begin();
    if (iPos < sExpr.size()){return iPos;}

    // Parse Expression
    exprtk::expression<double> exprTemp;
    exprTemp.register_symbol_table(symbol_table);
    parser.compile(exprStr, exprTemp);

    // Store
    vExpr.push_back(exprTemp); sExpr.push_back(exprStr);

    return vExpr.size() - 1;

}


/* double Parser::evaluateTime(int i, double nVal){ */

/*     // Update */
/*     varTime = nVal; return vExpr[i].value(); */

/* } */


/* double Parser::evaluateCoordinates(int i, double xCoord, double yCoord, double zCoord){ */
	
/*     // Update */
/* 	varX = xCoord; varY = yCoord; varZ = zCoord; return vExpr[i].value(); */

/* } */
