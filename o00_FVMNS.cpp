// Imports
#include <cmath>
#include <math.h>
#include <iostream>
#include <stdexcept>
#include <json/json.h>

// Self-Imports
#include "o01_Material.h"
#include "o02_Mesh.h"
/* #include "o03_Discretizer.h" */
/* #include "o04_Solver.h" */
/* #include "o04_CG.h" */
/* #include "o04_BCG.h" */
/* #include "o05_Probe.h" */
#include "o09_Parser.h"
#include "o09_Debugger.h"
/* #include "o09_Medic.h" */

namespace {
    template <size_t nDim> void runSolverLoop(const Json::Value& data){

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
        MeshSolver<nDim> p{}; Msh.generateMeshSolver(p, data["N"], data["sections"], data["refinement"], data["obstacles"]); std::cout << "Pressure object created with " << p.totNodes << " nodes and " << p.Obs.size() << " obstacles.\n";
        Msh.addBoundariesSolver(p, Mat, Prs, data["boundariesPressure"], Mat.P0, Mat.sP0); std::cout << p.BC.size() << " Pressure boundary conditions added.\n";
        
        // Temperature
        if (!data["T0"].isNull()) {
            MeshSolver<nDim> T{}; Msh.generateMeshSolver(T, data["N"], data["sections"], data["refinement"], data["obstacles"]); std::cout << "Temperature object created with " << T.totNodes << " nodes and " << T.Obs.size() << " obstacles.\n";
            Msh.addBoundariesSolver(T, Mat, Prs, data["boundariesTemperature"], Mat.T0, Mat.sT0); std::cout << T.BC.size() << " Temperature boundary conditions added\n"; }
        
        // Velocity
        std::array<MeshBase<nDim>, nDim> V{}; Msh.deriveMeshBase(p, V); std::cout << "Velocity objects created with "; for (MeshBase<nDim> Vk : V) {std::cout << Vk.totNodes << ", ";} std::cout << "\b\b nodes.\n";
        /* Msh.addBoundariesBase(V, Mat, Prs, data["boundariesVelocity"], Mat.VF0, Mat.sVF0); for (MeshBase<nDim> Vk : V) {std::cout << Vk.BC.size() << ", ";} std::cout << "\b boundary conditions added.\n"; */

        for (size_t i = 0; i < nDim; i++) { Msh.addBoundariesBase(i, V, Mat, Prs, data["boundariesVelocity"], Mat.VF0, Mat.sVF0); }

        /// Debug Current -- BOUNDARIES VELOCITY
        // Options
        debugOptions dOps{}; dOps.bGeneral = true; dOps.bBoundaries = true;

        // Print
        for (MeshBase<nDim> Vk : V) { printDebug(Vk, dOps); }

        return;

        // PENDING -- FINISH NAVIER-STOKES SOLVER BY COMPLETING ALL OTHER OBJECTS

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
        
        // temporal loop with doubles creates floating point errors after ~ 1e3 iterations
        // change to for (size_t i = 1; i < Dsc.endTime / Dsc.dt; i++) {t = dt * i;}
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
}

void runFVMNS(const Json::Value& data) {

    ///// Simulation /////
    size_t nDim = data["N"].size();
    switch (nDim) {
        case 1: runSolverLoop<1>(data); break;
        case 2: runSolverLoop<2>(data); break;
        case 3: runSolverLoop<3>(data); break;
        default: throw std::runtime_error("Unrecognized number of dimensions: " + std::to_string(nDim));
    }

}
