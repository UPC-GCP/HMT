// Imports
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <json/json.h>
#include <math.h>
#include <algorithm>
#include <chrono>
#include <thread>
#include <iomanip>

// Self-Imports
#include "Material.h"
#include "Mesh.h"
#include "Discretizer.h"
#include "Solver.h"
#include "GS.h"
#include "CG.h"
#include "ExpressionParser.h"
#include "Probe.h"
#include "Medic.h"

Json::Value getParsedData(std::string fileName){
    
    // Open File
    std::ifstream file(fileName, std::ifstream::binary);

    // Filter
    if (!file.is_open()) {
    std::cerr << "Error: Could not open the file " << fileName << std::endl;
    return 1;
    }

    // Parsing
    Json::Value data;
    Json::CharReaderBuilder readerBuilder;
    std::string errs;
    Json::parseFromStream(readerBuilder, file, &data, &errs);

    // Close File
    file.close();

    return data;

}

int main(int argc, char* argv[]){

    // Time
    auto t1 = std::chrono::high_resolution_clock::now();
    
    std::cout << "Initializing model ... \n";
    
    ////////// Config File //////////
    std::cout << "Reading data ... \n";

    // Read Config
    Json::Value data = getParsedData(argv[1]);
    std::cout << "Data parsed successfully. \n";



    ////////// Model Implementation //////////

    ///// Material /////
    std::cout << "Initializing Materials ...\n";
    Material Mat(data["materials"]); std::cout << "Material properties set.\n";
    Mat.setInitialConditions(data["T0"].asDouble()); std::cout << "Initial conditions set.\n";


    ///// Parser /////
    std::cout << "Initializing Parser ...\n";
    ExpressionParser Prs; std::cout << "Parser configured.\n";


    ///// Mesh /////
    std::cout << "Initializing mesh ...\n";
    Mesh Msh(data["meshAlgorithm"].asInt(), data["width"].asDouble(), data["height"].asDouble(), data["strength"].asDouble(), data["centering"].asDouble(), data["kappa"].asDouble(), data["delta"].asDouble()); std::cout << "Mesh parameters set.\n";
    Msh.newGenerateMesh(Mat, data["N"], data["sections"], data["refinement"]); std::cout << "Mesh created with " << Msh.totNodes << " nodes.\n";
    Msh.newAddBoundaryConditions(data["boundaries"], Prs); std::cout << Msh.newBoundaryConditions.size() << " boundary conditions added.\n";

    
    ///// Discretizer /////
    std::cout << "Initializing discretizer ...\n";
    Discretizer Dsc(data["scheme"].asString(), data["endTime"].asDouble(), data["timeStep"].asDouble());
    Dsc.setSchemeParameters(Mat, Msh); std::cout << "Temporal parameters set.\n";
    Dsc.newSetBoundaryConditions(Mat, Msh, Prs); std::cout << "Boundary conditions set.\n";
    Dsc.newSetCoefficients(Mat, Msh); std::cout << "Discretized coefficients set.\n";

    
    ///// Probes /////
    std::cout << "Initializing probes ...\n";
    Probe Prb(Msh, data["probes"], Dsc.scheme, argv[1]); std::cout << "Files configured.\n";
    Prb.checkProbes(Msh);


    ///// Medic /////
    std::cout << "Initializing medic ...\n";
    Medic Mdc(Msh, Prb); std::cout << "Diagnostic tools configured.";


    ///// Solver /////
    std::cout << "Initializing solver ... \n";
    Solver* Sol = nullptr;
    if (data["solver"] == "CG"){
        Sol = new CG(Dsc.scheme, data["maxIterations"].asDouble(), data["tolNumeric"].asDouble(), data["tolTemporal"].asDouble(), argv[1], data["solver"].asString());
    } else if (data["solver"] == "GS"){
        // Sol = new GS(Dsc.scheme, data["maxIterations"].asDouble(), data["tolNumeric"].asDouble(), data["tolTemporal"].asDouble(), argv[1], data["solver"].asString());
        std::cerr << "Currently unavailable.\n";
    } else {
        std::cerr << "Error: Invalid linear solver selected " << data["solver"].asString() << "\n";
    } std::cout << "Solver configured.\n";


    ////////// Temporal Loop //////////
    std::cout << "Processing ...\n";
    std::cout << std::fixed << std::setprecision(2);

    std::vector<std::vector<double>> cTemp{};
    for (double t = Dsc.dt; t <= Dsc.endTime; t += Dsc.dt){

        // Control
        cTemp = Msh.nT;

        // Update Coefficients
        Dsc.newSetBoundaryConditions(Mat, Msh, Prs, t);
        Dsc.newSetRHS(Mat, Msh);

        // Solver
        Sol->newSolve(Msh.matA, Msh.nT, Msh.bp, Msh.nIgnore);

        // Update Coefficients
        Dsc.newSetBoundaryConditions(Mat, Msh, Prs, t);
        Dsc.newSetRHS(Mat, Msh);

        // Diagnostics
        Mdc.getDiagnostic(Mat, Msh, Dsc, cTemp, t);
        Mdc.getSystemResidual(Mat, Msh, Dsc);

        // Write Data
        Prb.checkProbes(Msh, t);
        std::cout << "\r" << double(100 * t / Dsc.endTime) << " %";

        // Convergence
        if (std::sqrt(Sol->calcErr(cTemp, Msh.nT)) < data["tolTemporal"].asDouble()){std::cout << "\nSteady-state achieved @ t = " << t << " seconds."; break;}

    } std::cout << "\n";

    // Global Energy Balance
    Mdc.getGlobalBalance(Mat, Msh, Dsc);

    // Time
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> msDoub = t2 - t1;
    double tTime = msDoub.count()/1000/60;

    // End
    std::cout << "Time elapsed: " << int(tTime) << " minutes and " << (tTime - int(tTime))*60 << " seconds.\n";
    std::cout << "Files saved to: " << Prb.dirName;

}
