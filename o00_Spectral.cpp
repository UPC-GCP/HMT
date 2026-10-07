// Imports
#include <cmath>
#include <cstddef>
#include <math.h>
#include <iostream>
/* #include <stdexcept> */
#include <json/json.h>

// Self-Imports
#include "o01_Material.h"
#include "o02_MeshSpectral.h"
#include "o03_DiscretizerSpectral.h"
#include "o04_SolverSpectral.h"
#include "o05_ProbeSpectral.h"
/* #include "o06_MedicSpectral.h" */
/* #include "o09_Debugger.h" */

namespace {
    void runSolverLoop(const Json::Value& data, std::string configName) {

        ///// Material /////
        Material Mat(data["Re"].asDouble());

        ///// Mesh /////
        std::cout << "Initializing mesh ...\n"; 
        MeshSpectral Msh;

        // Burgers
        MeshBurgers Burg; Msh.generateMeshBurgers(Burg, data["N"].asInt()); std::cout << "Burgers mesh truncated at " << Burg.totNodes << " modes.\n";

        ///// Discretizer /////
        std::cout << "Initializing discretizer ...\n";
        DiscretizerSpectral Dsc; Dsc.initializeBurgers(Burg); std::cout << "Model initialized.\n";

        ///// Solver /////
        std::cout << "Initializing solver ...\n";
        SolverSpectral Sol(data["solver"].asString(), data["tolNumeric"].asDouble(), data["maxIterations"].asDouble()); std::cout << "Numerical solver configured.\n";

        ///// Probe /////
        std::cout << "Initializing probe ...\n";
        ProbeSpectral Prb(Burg, data["probes"], configName); std::cout << "Files stored at: " << Prb.fPath << "\n";

        // DEBUGGING HERE -- Finish adding probes

        return;

        // Add checkProbe() and store t = 0
        /* probeStoreData(Msh.uHat, Msh.E); */

        ///// Temporal Loop /////
        double dt = data["dt"].asDouble(), endTime = data["endTime"].asDouble(), tolTemporal = data["tolTemporal"].asDouble(), t{}; size_t iMax = endTime / dt;

        std::cout << "Processing ...\n";
        for (size_t i = 0; i < iMax; i++) {

            // Control
            t = i * dt; Burg.ouHat = Burg.uHat; Burg.oE = Burg.E;
            
            // Solver
            Sol.solveRK4(Burg.uHat, Mat.Re, t, dt);

            // Energy Balance
            Dsc.calculateEnergy(Burg.E, Burg.uHat, Mat.Re);

            // PENDING -- SETUP PROBE, SETUP MEDIC, END BURGERS DNS

            /* // Write Data */
            /* /1* probeStoreData(Msh.uHat, Msh.E); *1/ */
            /* std::cout << "\r" << double(100 * static_cast<double>(i) / iMax) << " %"; */

            /* // Diagnostic */
            // This one will also check Energy Transport Equation 
            /* medicRunDiagnostics(Msh.uHat); */

            // Convergence
            if (calcErr(Burg.uHat, Burg.ouHat) / dt < tolTemporal) { std::cout << "\nSteady-state achieved @ t = " << t << " seconds."; break; }

        } std::cout << "\n";

        // End
        /* std::cout << "Files saved to: \n"; */
    }

}

void runSpectral(const Json::Value& data, std::string configName) {
    ///// Simulation /////
    runSolverLoop(data, configName);
}
