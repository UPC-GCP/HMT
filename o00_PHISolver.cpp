// Imports
#include <ios>
#include <iostream>
#include <string>
#include <fstream>
#include <json/json.h>
#include <math.h>
#include <chrono>
#include <iomanip> // (setprecision)

// Self-Imports
#include "o01_Material.h"
#include "o02_Mesh.h"
#include "o03_Discretizer.h"
#include "o04_Solver.h"
#include "o04_GS.h"
#include "o04_CG.h"
#include "o04_BCG.h"
#include "o05_Probe.h"
#include "o09_ExpressionParser.h"
#include "o09_Medic.h" // Diagnostic Tool: Activate if needed

Json::Value getParsedData(std::string fileName){
    
    // Open File
    std::ifstream file(fileName, std::ifstream::binary);

    // Filter
    if (!file.is_open()) {
    std::cerr << "Error: Could not open the file " << fileName << std::endl;
    return 0;
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
    std::cout << std::defaultfloat << std::setprecision(3);
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Initializing model ... \n";
    

    ////////// Config File //////////
    std::cout << "Reading data ... \n";
    Json::Value data = getParsedData(argv[1]);
    std::cout << "Data parsed successfully. \n";



    ////////// Model Implementation //////////

    ///// Material /////
    std::cout << "Initializing Materials ...\n";
    Material Mat(data["materials"]); std::cout << "Material properties set.\n";
    Mat.setInitialConditions(data["PHI0"].asDouble()); std::cout << "Initial conditions set.\n";


    ///// Parser /////
    std::cout << "Initializing Parser ...\n";
    ExpressionParser Prs; std::cout << "Parser configured.\n";


    ///// Mesh /////
    std::cout << "Initializing mesh ...\n";
    Mesh Msh(data["meshAlgorithm"].asInt(), data["width"].asDouble(), data["strength"].asDouble(), data["centering"].asDouble(), data["kappa"].asDouble(), data["delta"].asDouble()); std::cout << "Mesh parameters set.\n";
    Msh.generateMesh(Mat, data["N"], data["sections"], data["refinement"]); std::cout << "Mesh created with " << Msh.totNodes << " nodes.\n";
    Msh.addBoundaryConditions(data["boundaries"], Mat, Prs); std::cout << Msh.boundaryConditions.size() << " boundary conditions added.\n";
    Msh.addVelocityField(data["velocityField"], Prs); std::cout << "Velocity field generated.\n";

    
    ///// Discretizer /////
    std::cout << "Initializing discretizer ...\n";
    Discretizer Dsc(data["tempScheme"].asString(), data["spatScheme"].asString(), data["endTime"].asDouble(), data["timeStep"].asDouble());
    Dsc.setSchemeParameters(Mat, Msh); std::cout << "Temporal parameters set.\n";
    Dsc.newSetBoundaries(Mat, Msh, Prs); std::cout << "Boundary node coefficiens set.\n";
    Dsc.newSetCoefficients(Mat, Msh); std::cout << "Interior node coefficients set.\n";

    
    ///// Solver /////
    std::cout << "Initializing solver ... \n";
    Solver* Sol = nullptr;
    if (data["solver"] == "CG"){
        Sol = new CG(Dsc.tempScheme, data["maxIterations"].asDouble(), data["tolNumeric"].asDouble(), data["tolTemporal"].asDouble(), argv[1], data["solver"].asString());
    } else if (data["solver"] == "GS"){
        Sol = new GS(Dsc.tempScheme, data["maxIterations"].asDouble(), data["tolNumeric"].asDouble(), data["tolTemporal"].asDouble(), argv[1], data["solver"].asString());
    } else if (data["solver"] == "BCG"){
        Sol = new BCG(Dsc.tempScheme, data["maxIterations"].asDouble(), data["tolNumeric"].asDouble(), data["tolTemporal"].asDouble(), argv[1], data["solver"].asString());

    } else {
        std::cerr << "Error: Invalid linear solver selected " << data["solver"].asString() << "\n";
    } std::cout << "Solver configured.\n";


    ///// Probes /////
    std::cout << "Initializing probes ...\n";
    Probe Prb(Msh, data["probes"], Dsc.tempScheme, Dsc.spatScheme, argv[1]); std::cout << "Files configured.\n";
    Prb.checkProbes(Msh, Sol);


    ///// Medic /////
    std::cout << "Initializing medic ...\n";
    bool bMdc = data["medicOn"].asBool();
    Medic Mdc(Msh, Prb, bMdc); std::cout << "Diagnostic tools configured.\n";



    ////////// Temporal Loop //////////
    std::cout << "Processing ...\n";
    
    std::vector<std::vector<double>> cPhi{}, cbPhi{};
    for (double t = Dsc.dt; t <= Dsc.endTime; t += Dsc.dt){

        // Control
        cPhi = Msh.vPhi; cbPhi = Msh.bcPhi;
        Msh.oPhi = Msh.vPhi; Msh.obPhi = Msh.bcPhi;
        
        // Solver
        if (!Sol->newSolve(Msh.matA, Msh.vPhi, Msh.bp, Msh.nIgnore)){std::cerr << "Simulation diverges @ t = " << t; break;}

        // Diagnostics
        if (bMdc){
            Mdc.getDiagnostic(Mat, Msh, Dsc, Prs, t);
            Mdc.getSystemResidual(Mat, Msh, Dsc, t);
        }

        // Write Data
        Prb.checkProbes(Msh, Sol, t);
        std::cout << "\r" << double(100 * t / Dsc.endTime) << " %";
	
	    // Update Coefficients
        Dsc.newSetBoundaries(Mat, Msh, Prs, t);
        Dsc.newSetCoefficients(Mat, Msh);

        // Convergence
        if (std::sqrt(Sol->calcErr(cPhi, Msh.vPhi)) < data["tolTemporal"].asDouble()){std::cout << "\nSteady-state achieved @ t = " << std::setprecision(2) << t << " seconds."; break;}

    } std::cout << "\n";

    // Global Energy Balance
    /* Mdc.getGlobalBalance(Mat, Msh, Dsc); */

    // Time
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> msDoub = t2 - t1;
    double tTime = msDoub.count()/1000/60;

    // End
    std::cout << "Time elapsed: " << int(tTime) << " minutes and " << (tTime - int(tTime))*60 << " seconds.\n";
    std::cout << "Files saved to: " << Prb.dirName << "\n";

}

