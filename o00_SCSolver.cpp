// Imports
#include <iostream>
#include <json/forwards.h>
#include <string>
#include <fstream>
#include <json/json.h>
#include <math.h>
#include <chrono>
#include <iomanip> // (setprecision)

// Self-Imports
#include "o01_Material.h"
#include "o02_MeshSC.h"
/* #include "o03_DiscretizerNS.h" */
/* #include "o04_SolverNS.h" */
/* #include "o04_CGNS.h" */
/* #include "o04_BCG.h" */
/* #include "o05_ProbeNS.h" */
/* #include "o09_ExpressionParser.h" */
/* #include "o09_MedicNS.h" */


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
    std::cout << std::fixed << std::setprecision(3); // std::fixed, std::defaultfloat, std::scientific? (not sure about the last one)
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
    Mat.setInitialConditions(data["PHI0"].asDouble(), data["VF0"]); std::cout << "Initial conditions set.\n";


    ///// Parser /////
    std::cout << "Initializing Parser ...\n";
    ExpressionParser Prs; std::cout << "Parser configured.\n";


    ///// Mesh /////
    std::cout << "Initializing mesh ...\n";
    Mesh Msh(data["meshAlgorithm"].asInt(), data["width"].asDouble(), data["strength"].asDouble(), data["centering"].asDouble(), data["kappa"].asDouble(), data["delta"].asDouble()); std::cout << "Mesh parameters set.\n";
    Msh.generateMesh(Msh.p, data["PHI0"].asDouble(), data["N"], data["sections"], data["refinement"], data["obstacles"]); std::cout << "Primary mesh created with " << Msh.p.totNodes << " nodes.\n";
    /* Msh.generateMeshVelocity(Mat, Msh.p, Msh.u, Msh.v); std::cout << "Secondary meshes created with " << Msh.u.totNodes << " and " << Msh.v.totNodes << " nodes.\n"; */
    /* Msh.addBoundariesPressure(data["boundaries"], Mat, Prs); std::cout << Msh.boundaryConditions.size() << " boundary conditions added.\n"; */
    /* Msh.addBoundariesVelocity(data["boundariesVelocity"], Mat); std::cout << Msh.boundaryVelocity.size() << " velocity boundary conditions added.\n"; */

    std::cout << "nMat:\n";
    for (std::vector<int> vec : Msh.p.nMat){
        for (int val : vec){std::cout << val << " ";} std::cout << "\n";
    } std::cout << "\n";

    std::cout << "bObs:\n";
    for (std::vector<bool> vec : Msh.p.bObs){
        for (bool val : vec){std::cout << val << " ";} std::cout << "\n";
    } std::cout << "\n";

    return 0;
    
    /* /1* ///// Discretizer ///// *1/ */
    /* std::cout << "Initializing discretizer ...\n"; */
    /* Discretizer Dsc(data["tempScheme"].asString(), data["spatScheme"].asString(), data["endTime"].asDouble(), data["timeStep"].asDouble()); std::cout << "Discretizer parameters set.\n"; */
    /* Dsc.setSchemeParameters(Mat, Msh); std::cout << "Scheme parameters set.\n"; */
    /* Dsc.setMomentumBoundaries(Mat, Msh); std::cout << "Velocity boundaries set.\n"; */
    /* Dsc.setMomentumCoefficients(Mat, Msh); Dsc.setMomentumBoundaries(Mat, Msh); std::cout << "Velocity predictor set.\n"; */
    /* Dsc.setPressureBoundaries(Mat, Msh); std::cout << "Pressure boundaries set.\n"; */
    /* Dsc.setPressureCoefficients(Mat, Msh); Dsc.setPressureBoundaries(Mat, Msh); std::cout << "Pressure coefficients set.\n"; */
    

    /* ///// Solver ///// */
    /* std::cout << "Initializing solver ... \n"; */
    /* Solver* Sol = nullptr; */
    /* if (data["solver"] == "CG"){ */
    /*     Sol = new CG(Dsc.tempScheme, data["maxIterations"].asDouble(), data["tolNumeric"].asDouble(), data["tolTemporal"].asDouble(), argv[1], data["solver"].asString()); */
    /* } else if (data["solver"] == "GS"){ */
    /*     // Sol = new GS(Dsc.scheme, data["maxIterations"].asDouble(), data["tolNumeric"].asDouble(), data["tolTemporal"].asDouble(), argv[1], data["solver"].asString()); */
    /*     std::cerr << "Currently unavailable.\n"; */
    /* } else if (data["solver"] == "BCG"){ */
    /*     /1* Sol = new BCG(Dsc.tempScheme, data["maxIterations"].asDouble(), data["tolNumeric"].asDouble(), data["tolTemporal"].asDouble(), argv[1], data["solver"].asString()); *1/ */
    /* } else { */
    /*     std::cerr << "Error: Invalid linear solver selected " << data["solver"].asString() << "\n"; */
    /* } std::cout << "Solver configured.\n"; */


    /* /1* ///// Probes ///// *1/ */
    /* std::cout << "Initializing probes ...\n"; */
    /* Probe Prb(Msh, data["probes"], Dsc.tempScheme, Dsc.spatScheme, argv[1]); std::cout << "Files configured.\n"; */
    /* Prb.checkProbes(Msh, Sol); std::cout << "Initial conditions stored.\n"; */


    /* /1* ///// Medic ///// *1/ */
    /* std::cout << "Initializing medic ...\n"; */
    /* bool bMdc = data["medicOn"].asBool(); */
    /* Medic Mdc(Msh, Prb, bMdc); std::cout << "Diagnostic tools configured.\n"; */



    /* ////////// Temporal Loop ////////// */
    /* std::cout << "Processing ...\n"; */
    
    /* for (double t = Dsc.dt; t <= Dsc.endTime; t += Dsc.dt){ */

    /*     // Control */
    /*     Msh.p.oPhi = Msh.p.Phi; */
        
    /*     // Solver */
    /*     if (!Sol->newSolve(Msh.p.matA, Msh.p.Phi, Msh.p.matB)){std::cerr << "Simulation diverges @ t = " << t; break;} */
    /*     if (std::sqrt(Sol->lastRes) >= data["tolNumeric"].asDouble()){std::cerr << "\nWARN: CG unconverged @ t=" << t << " lastIter=" << Sol->lastIter << " lastRes=" << std::sqrt(Sol->lastRes);} */

    /*     // Correct Velocity */
    /*     Dsc.correctVelocity(Mat, Msh); */
    /*     Msh.u.oPhi = Msh.u.Phi; Msh.v.oPhi = Msh.v.Phi; */

    /*     // Diagnostics */
    /*     if (bMdc){ */
    /*         Mdc.getDiagnostic(Mat, Msh, Dsc, t); */
    /*         Mdc.getSystemResidual(Mat, Msh, Dsc, t); */
    /*     } */

    /*     // Write Data */
    /*     Prb.checkProbes(Msh, Sol, t); */
    /*     std::cout << "\r" << double(100 * t / Dsc.endTime) << " %"; */
	
	    /* // Update Coefficients */
    /*     Dsc.checkStability(Mat, Msh); */
    /*     Dsc.setMomentumBoundaries(Mat, Msh); */
    /*     Dsc.setMomentumCoefficients(Mat, Msh); */
    /*     Dsc.setPressureBoundaries(Mat, Msh); */
    /*     Dsc.setPressureCoefficients(Mat, Msh); */

    /*     // Convergence */
    /*     if (std::sqrt(Sol->calcErr(Msh.p.oPhi, Msh.p.Phi)) < data["tolTemporal"].asDouble()){std::cout << "\nSteady-state achieved @ t = " << std::setprecision(2) << t << " seconds."; break;} */

    /* } std::cout << "\n"; */

    /* // Global Energy Balance */
    /* Mdc.getGlobalBalance(Mat, Msh, Dsc); */

    /* // Time */
    /* auto t2 = std::chrono::high_resolution_clock::now(); */
    /* std::chrono::duration<double, std::milli> msDoub = t2 - t1; */
    /* double tTime = msDoub.count()/1000/60; */

    /* // End */
    /* std::cout << "Time elapsed: " << int(tTime) << " minutes and " << (tTime - int(tTime))*60 << " seconds.\n"; */
    /* std::cout << "Files saved to: " << Prb.dirName << "\n"; */

}
