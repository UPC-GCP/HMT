// Imports
/* #include <cstddef> */
#include <cmath>
#include <cstddef>
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
#include "o02_Mesh.h"
/* #include "o03_Discretizer.h" */
/* #include "o04_Solver.h" */
/* #include "o04_CG.h" */
/* #include "o04_BCG.h" */
/* #include "o05_Probe.h" */
#include "o09_Parser.h"
/* #include "o09_Medic.h" */


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

    ///// SOME IMPORTANT DETAILS ABOUT THIS MODEL IN PARTICULAR (AND HOW TO ADDRESS THEM) (HOPEFULLY) (WE'LL SEE AT THE END HOW MUCH OF THIS INITIAL PLANNING WAS SUCCESSFUL)
    // THIS ONE WILL BE THE FULL 3D MODEL AND SHOULD BE CAPABLE OF WORKING AT ALL DIMENSION LEVELS AND FOR ALL EXEPTIONS
    // CURRENTLY THE MODEL IS SEPARATED IN 4 SOLVERS:
        // MainSolver: Solver for a single scalar variable. Primitive architecture compared to the remaining models, but works fine.
        // NSSolver: Basic 2D Navier-Stokes solver. Solves a simplified version of the model, solves MASS EQUATION and MOMENTUM EQUATION.
        // DHCSolver: NSSolver with the addition of the ENERGY EQUATION.
        // SCSolver: NSSolver with the addition of OBSTACLES.
    // 3DSolver: 3D Navier-Stokes solver. Includes the full Navier-Stokes set of equations with 1 MASS EQUATION, 3 MOMENTUM EQUATIONS, and 1 ENERGY EQUATION. Will also be able to include obstacles.

    // Mesh will be rewritten to move everything to flattened arrays and just use l = i + N[1] * j + N[2] * k
    // Geometry arrays will be stored in dimension dependent vectors that grow with N.size()
    // Discretizer will have functions separated for each dimension and will activate the ones that are required
    // Numerical solver will also depend on the amount of dimensions

    // Need to generate Templates for all of these functions that depend on N.size()
    // Variables for each dimension need to be defined outside the array instead of inside

    // TRYING TO ORGANIZE THE IDEA
    // Index everything in dimensional arrays and have it set so iD = {0:W, 1:E, 2:S, 3:N, 4:B, 4:T}
    // Use that as reference and switch between each case
    // Single functions on each side that index with iD for each dimension
    // For each operation it should be evaluated across all active dimensions
    // Activate dimensions at the beginning by resizing everything accordingly (?)
    // Boundary conditions: for (i < iD.size()) {Evaluate boundary for i side}
    // Need coefficients and neighbours to be indexable with ease (inline int getIndex(i, j=1, k=1, Ny=1, Nz=1))
    // Generate mesh objects in Main -> Templates create efficient data structures for each without wasting space
    // Discretizer, Solver, Probe, Medic need to use the new structure and identify each neighbour within the arrays

    // Mesh
    // Objects in Mesh need to use Templates

    // Discretizer
    // Functions in Mesh need to use Templates

    // Use arrays for fixed dimensional objects
    // CANNOT DEFINE A TEMPLATE OBJECT IN THE HEADER FILE, NEED TO MOVE ALL OF THAT TO MAIN SOLVER
    // Can I have an array that calls different functions depending on what I want to write

template <size_t nDim> void runSimulation(Json::Value data){

    ///// Parser /////
    std::cout << "Initializing Parser ...\n";
    Parser Prs; std::cout << "Parser configured.\n";

    ///// Material /////
    std::cout << "Initializing Materials ...\n";
    Material Mat(data["materials"], data["g"].isNull() ? 9.81 : data["g"].asDouble()); std::cout << "Material properties set.\n";

    /* Mat.setInitialConditions(data["PHI0"].asDouble(), data["VF0"]); std::cout << "Initial conditions set.\n"; */

    /* if (!data["PHI0"].isNull()) {Mat.setInitialConditions(data["PHI0"].asDouble());} */

    // Initial conditions has 3 cases 
    // 1. Phi - Scalar map, takes a double or a .csv
    // 2. V - Fixed value, function, variable
    // 3. p, T - active/inactive



    return; 

    ///// Mesh /////
    std::cout << "Initializing mesh ...\n"; 
    Mesh Msh; std::cout << "Mesh initialized.\n";

    /* if (data["boundariesPressure"].size() != 0) { */
    /*     // Create Pressure (MeshSolver) */
    /* } */
    /* if (data["boundariesTemperature"].size() != 0) { */
    /*     // Create Temperature (MeshSolver) */
    /* } */
    /* if (data["boundariesVelocity"].size() != 0) { */
    /*     // Create Velocity (MeshBase) */
    /* } */


    // NEED TO CREATE MY VARIABLES HERE FOR THE DIFFERENT CASES DEPENDING ON THE CONFIGURATION FILE
    // Solver Variables: p, T
    // Explicit Variables: u, v, w
    // Boundaries could be stored within each template for the mesh
    // Obstacles defined independently
    
    // Velocity components needs to be in a single vector
    /* std::array<MeshBase<Dim>, Dim> V{}; */




    /* Mesh Msh(data["meshAlgorithm"].asInt(), data["width"].asDouble(), data["strength"].asDouble(), data["centering"].asDouble(), data["kappa"].asDouble(), data["delta"].asDouble()); std::cout << "Mesh parameters set.\n"; */



    /* Msh.generateMesh(Msh.p, data["P0"].asDouble(), data["N"], data["sections"], data["refinement"], data["obstacles"]); std::cout << "Primary mesh created with " << Msh.p.totNodes << " nodes and " << Msh.obstacles.size() << " obstacles.\n"; */
    /* Msh.generateMeshVelocity(Mat, Msh.p, Msh.u, Msh.v); std::cout << "Secondary meshes created with " << Msh.u.totNodes << " and " << Msh.v.totNodes << " nodes.\n"; */
    /* Msh.addBoundariesPressure(data["boundaries"], Mat, Prs); std::cout << Msh.boundaryPressure.size() << " boundary conditions added.\n"; */
    /* Msh.addBoundariesVelocity(data["boundariesVelocity"], Mat); std::cout << Msh.boundaryVelocity.size() << " velocity boundary conditions added.\n"; */



    /* ///// Discretizer ///// */
    /* std::cout << "Initializing discretizer ...\n"; */
    /* Discretizer Dsc(data["tempScheme"].asString(), data["spatScheme"].asString(), data["endTime"].asDouble(), data["timeStep"].asDouble()); std::cout << "Discretizer parameters set.\n"; */
    /* Dsc.setSchemeParameters(Mat, Msh); std::cout << "Scheme parameters set.\n"; */
    /* Dsc.setMomentumBoundaries(Mat, Msh); std::cout << "Velocity boundaries set.\n"; */
    /* Dsc.setMomentumCoefficients(Mat, Msh); Dsc.setMomentumBoundaries(Mat, Msh); std::cout << "Velocity predictor set.\n"; Dsc.setPressureBoundaries(Mat, Msh); std::cout << "Pressure boundaries set.\n"; Dsc.setPressureCoefficients(Mat, Msh); Dsc.setPressureBoundaries(Mat, Msh); std::cout << "Pressure coefficients set.\n"; */
    
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

    /* ///// Temporal Loop ///// */
    /* std::cout << "Processing ...\n"; */
    
    /* for (double t = Dsc.dt; t <= Dsc.endTime; t += Dsc.dt){ */

    /*     // Control */
    /*     Msh.p.oPhi = Msh.p.Phi; */
        
    /*     // Solver */
    /*     if (!Sol->newSolve(Msh.p.matA, Msh.p.Phi, Msh.p.matB, Msh.p.bObs)){std::cerr << "Simulation diverges @ t = " << t; break;} */
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
    /*     if (std::sqrt(Sol->calcErr(Msh.p.oPhi, Msh.p.Phi, Msh.p.bObs)) < data["tolTemporal"].asDouble()){std::cout << "\nSteady-state achieved @ t = " << std::setprecision(2) << t << " seconds."; break;} */

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


int main(int argc, char* argv[]){
    
    ////////// Configuration //////////
    
    ///// Setup /////
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Initializing model ... \n" << std::fixed << std::setprecision(3); // std::fixed, std::defaultfloat, std::scientific? (not sure about the last one)

    ///// Data /////
    std::cout << "Reading data ... \n";
    Json::Value data = getParsedData(argv[1]); std::cout << "Data parsed successfully. \n";

    ///// Simulation /////
    size_t nDim = data["N"].size();
    switch (nDim) {
        case 1:
            runSimulation<1>(data); break;
        case 2:
            runSimulation<2>(data); break;
        case 3:
            runSimulation<3>(data); break;
        default:
            std::cerr << "Configuration file not defined properly.\n"; return 1;
    }

}

    /* std::cout << "xFaces: "; */
    /* for (double val : Msh.p.Faces[0]){std::cout << val << " ";} std::cout << "\n"; */

    /* std::cout << "xNodes: "; */
    /* for (double val : Msh.p.Nodes[0]){std::cout << val << " ";} std::cout << "\n"; */

    /* std::cout << "xDelta: "; */
    /* for (double val : Msh.p.deltaX[0]){std::cout << val << " ";} std::cout << "\n"; */

    /* std::cout << "xD :"; */
    /* for (double val : Msh.p.dX[0]){std::cout << val << " ";} std::cout << "\n"; */

    /* std::cout << "yFaces: "; */
    /* for (double val : Msh.p.Faces[1]){std::cout << val << " ";} std::cout << "\n"; */

    /* std::cout << "yNodes: "; */
    /* for (double val : Msh.p.Nodes[1]){std::cout << val << " ";} std::cout << "\n"; */

    /* std::cout << "yDelta: "; */
    /* for (double val : Msh.p.deltaX[1]){std::cout << val << " ";} std::cout << "\n"; */

    /* std::cout << "yD :"; */
    /* for (double val : Msh.p.dX[1]){std::cout << val << " ";} std::cout << "\n"; */

    /* std::cout << "zFaces: "; */
    /* for (double val : Msh.p.Faces[2]){std::cout << val << " ";} std::cout << "\n"; */

    /* std::cout << "zNodes: "; */
    /* for (double val : Msh.p.Nodes[2]){std::cout << val << " ";} std::cout << "\n"; */

    /* std::cout << "zDelta: "; */
    /* for (double val : Msh.p.deltaX[2]){std::cout << val << " ";} std::cout << "\n"; */

    /* std::cout << "zD :"; */
    /* for (double val : Msh.p.dX[2]){std::cout << val << " ";} std::cout << "\n"; */
