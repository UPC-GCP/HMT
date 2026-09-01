// Imports
#include <cstddef>
#include <iostream>
#include <json/json.h>
#include <chrono>
#include <iomanip> // (setprecision)
#include <cstddef>
#include <cmath>
#include <cstdlib>
/* #include <json/forwards.h> */
#include <stdexcept>
/* #include <string> */
/* #include <fstream> */
#include <math.h>

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
#include "o09_Debugger.h"

////////// JSON Parser ///////////
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

////////// Debugging Tools //////////

template <size_t nDim> void printDebug(MeshSolver<nDim> Msh, debugOptions dOps) {

    // General Mesh
    for (size_t i = 0; i < nDim; i++) {
        std::cout << "\nAxis: " << i << "\n";
        std::cout << "Faces: "; for (double val : Msh.Faces[i]) {std::cout << val << " ";} std::cout << "\n";
        std::cout << "Nodes: "; for (double val : Msh.Nodes[i]) {std::cout << val << " ";} std::cout << "\n";
        std::cout << "DeltaX: "; for (double val : Msh.deltaX[i]) {std::cout << val << " ";} std::cout << "\n";
        std::cout << "dX: "; for (double val : Msh.dX[i]) {std::cout << val << " ";} std::cout << "\n";
    } std::cout << "\n";

    // Dimensional
    if constexpr (nDim == 1) {print1D(Msh, dOps);}
    else if constexpr (nDim == 2) {print2D(Msh, dOps);}
    else if constexpr (nDim == 3) {print3D(Msh, dOps);}

    // General Boundaries
    for (size_t i = 0; i < Msh.BC.size(); i++) {
        std::cout << "\nBC: " << i << "\n";
        std::cout << "Type - Side: " << Msh.BC[i].type << " " << Msh.BC[i].side << "\n";
        std::cout << "Expression: " << Msh.BC[i].bUpdate << " - " << Msh.BC[i].iExpr << " " << Msh.BC[i].iEq << " " << Msh.BC[i].expression << "\n";
        std::cout << "ExpressionA: " << Msh.BC[i].bA << " - " << Msh.BC[i].iExprA << " " << Msh.BC[i].iEqA << " " << Msh.BC[i].expressionA << "\n";
        std::cout << "Pos0: "; for (size_t val : Msh.BC[i].i0) { std::cout << val << " "; } std::cout << "\n";
        std::cout << "Pos1: "; for (size_t val : Msh.BC[i].i1) { std::cout << val << " "; } std::cout << "\n";

        // Dimensional
        // TERMINAR ESTO MAÑANA

        
    }

}

////////// Solvers //////////
bool bRun = false;

template <size_t nDim> void runNSSolver(Json::Value data){
    // Control
    bRun = true;

    ///// Parser /////
    std::cout << "Initializing Parser ...\n";
    Parser Prs; std::cout << "Parser configured.\n";

    ///// Material /////
    std::cout << "Initializing Materials ...\n";
    Material Mat(data["materials"], data["g"].isNull() ? 9.81 : data["g"].asDouble()); std::cout << "Material properties set.\n";
    
    // Initial Conditions
    if (data["P0"].isDouble()) {
        // Control
        if (!data["T0"].isDouble() && !data["T0"].isNull()) {std::cerr << "Initial conditions (p/T) not defined properly.\n"; throw std::invalid_argument("Check .json");}
        for (Json::Value::ArrayIndex i = 0; i < data["V0"].size(); i++) {if (!data["V0"][i].isDouble()) {std::cerr << "Initial conditions (V) not defined properly.\n"; throw std::invalid_argument("Check .json");}}

        // Store data
        Mat.setInitialConditions(data["T0"].isNull() ? 0 : data["T0"].asDouble(), data["P0"].asDouble(), data["VF0"]);
    } else if (data["P0"].isString()) {
        // Control
        if (!data["T0"].isString() && !data["T0"].isNull()) {std::cerr << "Initial conditions (p/T) not defined properly.\n"; throw std::invalid_argument("Check .json");}
        for (Json::Value::ArrayIndex i = 0; i < data["V0"].size(); i++) {if (!data["V0"][i].isString()) {std::cerr << "Initial conditions (V) not defined properly.\n"; throw std::invalid_argument("Check .json");}}

        // Store data
        Mat.setInitialConditions(data["T0"].asString(), data["P0"].asString(), data["VF0"]);
    } else {std::cerr << "Initial conditions not defined correctly.\n"; throw std::invalid_argument("Check .json");} std::cout << "Initial conditions logged.\n";

    ///// Mesh /////
    std::cout << "Initializing mesh ...\n"; 
    Mesh<nDim> Msh;

    // Pressure
    std::cout << data["obstacles"].size() << " obstacles identified.\n";
    MeshSolver<nDim> p{}; Msh.generateMeshSolver(p, Mat, data["N"], data["sections"], data["refinement"], data["obstacles"]); std::cout << "Pressure object created with " << p.totNodes << " nodes and " << p.Obs.size() << " obstacles.\n";
    Msh.addBoundariesSolver(p, Mat, Prs, data["boundariesPressure"], Mat.P0, Mat.sP0); std::cout << p.BC.size() << " Pressure boundary conditions added.\n";

    /// Debug Current -- DEBUGGING BOUNDARY CONDITIONS
    // Options
    debugOptions dOps{}; dOps.bPhi = true;
    printDebug(p, dOps);

    std::cout << "Test end\n";
    return;
    
    // Temperature
    if (!data["T0"].isNull()) {
        MeshSolver<nDim> T{}; Msh.generateMeshSolver(T, Mat, data["N"], data["sections"], data["refinement"], data["obstacles"]); std::cout << "Temperature object created with " << T.totNodes << " nodes and " << T.Obs.size() << " obstacles.\n";
        Msh.addBoundariesSolver(T, Mat, Prs, data["boundariesTemperature"], Mat.T0, Mat.sT0); std::cout << T.BC.size() << " Temperature boundary conditions added\n";
    }
    
    /* // Velocity */
    /* std::array<MeshBase<nDim>, nDim> V{}; // V.generateMeshVelocity(); std::cout << "Velocity objects created with "; for(MeshBase<nDim> Vk : V) {std::cout << Vk.totNodes << " ";} std::cout << " nodes.\n"; */



    /* Msh.addBoundariesVelocity(data["boundariesVelocity"], Mat); std::cout << Msh.boundaryVelocity.size() << " velocity boundary conditions added.\n"; */




    return;


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

    /* // End */
    /* std::cout << "Files saved to: " << Prb.dirName << "\n"; */

}

template <size_t nDim> void runPHISolver(Json::Value data){
    // Cases
    // 1D Rod, 4Materials: Phi0 = Fixed, V = 0
    // Smith-Hutton: Phi0 = Fixed, V = function

    ///// Control /////
    bRun = true;
    
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
    
}

////////// MAIN //////////
int main(int argc, char* argv[]){
    
    ///// Setup /////
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Initializing model ... \n" << std::fixed << std::setprecision(3); // std::fixed, std::defaultfloat, std::scientific? (not sure about the last one)

    ///// Data /////
    std::cout << "Reading data ... \n";
    Json::Value data = getParsedData(argv[1]); std::cout << "Data parsed successfully. \n";

    ///// Simulation /////
    if (data["PHI0"].isNull() && data["P0"].isNull()) {std::cerr << "Configuration file not defined properly.\n"; return 1;}

    try {
        size_t nDim = data["N"].size();
        switch (nDim) {
            case 1:
                data["PHI0"].isNull() ? runNSSolver<1>(data) : runPHISolver<1>(data); break;
            case 2:
                data["PHI0"].isNull() ? runNSSolver<2>(data) : runPHISolver<2>(data); break;
            case 3:
                data["PHI0"].isNull() ? runNSSolver<3>(data) : runPHISolver<3>(data); break;
            default:
                std::cerr << "Number of dimensions not defined properly.\n"; return 1;
        }
    }
    catch (const std::exception& e) {std::cerr << "Program shutdown...\n"; return EXIT_FAILURE;}

    ///// Control /////
    if (!bRun) {std::cerr << "Configuration file not defined properly.\n"; return 1;}
    auto t2 = std::chrono::high_resolution_clock::now(); std::chrono::duration<double, std::milli> msDoub = t2 - t1; double tTime = msDoub.count()/1000/60; std::cout << "Time elapsed: " << int(tTime) << " minutes and " << (tTime - int(tTime))*60 << " seconds.\n";

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
    // CANNOT DEFINE A TEMPLATE OBJECT IN THE HEADER FILE, NEED TO MOVE ALL OF THAT TO MAIN SOLVER, not entirely true
    // Can I have an array that calls different functions depending on what I want to write

    // Model follows Patankar's discretization so it ALWAYS NEEDS VELOCITY
