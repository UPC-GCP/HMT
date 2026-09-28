// Imports
#include <cmath>
#include <math.h>
#include <iostream>
#include <stdexcept>
#include <json/json.h>
#include <string>

// Self-Imports
#include "o01_Material.h"
#include "o02_Mesh.h"

namespace {
    template <size_t nDim> void runSolverLoop(const Json::Value& data){
        
        ///// Parser /////
        std::cout << "Initializing Parser ...\n";
        Parser Prs; std::cout << "Parser configured.\n";

        ///// Material /////
        std::cout << "Initializing Materials ...\n";
        Material Mat(data["materials"], data["g"].isNull() ? 9.81 : data["g"].asDouble()); std::cout << "Material properties set.\n";
        
        // Initial Conditions
        if (data["PHI0"].isDouble()) {
            // Control
            for (Json::Value::ArrayIndex i = 0; i < data["V0"].size(); i++) {if (!data["V0"][i].isDouble() || !data["V0"][i].isString()) {std::cerr << "Initial conditions (V) not defined properly.\n"; throw std::invalid_argument("Check .json");}}

            // Store data
            Mat.setInitialConditions(data["PHI0"].asDouble(), data["V0"]);

        } else if (data["PHI0"].isString()) {
            // Control
            for (Json::Value::ArrayIndex i = 0; i < data["V0"].size(); i++) {if (!data["V0"][i].isDouble() || !data["V0"][i].isString()) {std::cerr << "Initial conditions (V) not defined properly.\n"; throw std::invalid_argument("Check .json");}}

            // Store data
            Mat.setInitialConditions(data["PHI0"].asString(), data["V0"]);
        } else {std::cerr << "Initial conditions not defined correctly.\n"; throw std::invalid_argument("Check .json");} std::cout << "Initial conditions logged.\n";

        ///// Mesh /////
        std::cout << "Initializing mesh ...\n"; 
        MeshSolver<nDim> PHI{};

        // PENDING: FINISH SCALAR SOLVER STRUCTURE, SHOULD BE QUICK ONCE ALL OBJECTS ARE DONE
        
    }
}

void runFVMScalar(const Json::Value& data) {

    ///// Simulation /////
    size_t nDim = data["N"].size();
    switch (nDim) {
        case 1: runSolverLoop<1>(data); break;
        case 2: runSolverLoop<2>(data); break;
        case 3: runSolverLoop<3>(data); break;
        default: throw std::runtime_error("Unrecognized number of dimensions: " + std::to_string(nDim));
    }

}
