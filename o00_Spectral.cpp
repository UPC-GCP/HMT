// Imports
#include <cmath>
#include <cstddef>
#include <iostream>
#include <json/json.h>
#include <optional>
#include <stdexcept>
#include <string>
#include <variant>
/* #include <type_traits> */

// Self-Imports
#include "o01_Material.h"
#include "o02_MeshSpectral.h"
#include "o03_DiscretizerSpectral.h"
#include "o04_Solver.h"
#include "o04_SolverSpectral.h"
#include "o05_ProbeSpectral.h"
#include "o06_MedicSpectral.h"

using MeshVariant = std::variant<MeshBurgers, MeshLES>;

namespace {
    void runSolverLoop(const Json::Value& data, std::string configName) {
        ///// Material /////
        Material Mat(data["Re"].asDouble());

        ///// Solver /////
        std::cout << "Initializing solver ...\n";
        SolverSpectral Sol(data["solver"].asString(), data["tolNumeric"].asDouble(), data["maxIterations"].asDouble()); std::cout << "Numerical solver configured.\n";

        ///// Mesh /////
        std::cout << "Initializing mesh ...\n"; 
        MeshSpectral Msh;

        // Burgers
        MeshVariant Burg;
        if (Sol.algorithm == solSpectral::DNS) { auto& vMsh = Burg.emplace<MeshBurgers>(); Msh.generateMeshBurgers(vMsh, data["N"].asInt()); }
        else if (Sol.algorithm == solSpectral::LES) { auto& vMsh = Burg.emplace<MeshLES>(); Msh.generateMeshBurgers(vMsh, data["N"].asInt(), data["CK"].asDouble(), data["slope"].asDouble()); }
        std::visit([&](auto& vMsh) { std::cout << "Burgers mesh truncated at " << vMsh.totNodes << " modes.\n"; }, Burg);

        ///// Discretizer /////
        std::cout << "Initializing discretizer ...\n"; DiscretizerSpectral Dsc; std::visit([&](auto& vMsh) { Dsc.initializeBurgers(vMsh, Mat.Re); }, Burg); std::cout << "Model initialized.\n";

        ///// Probe /////
        std::cout << "Initializing probe ...\n"; std::optional<ProbeSpectral> Prb; std::visit([&](auto& vMsh) { Prb.emplace(vMsh, data["probes"], configName); (*Prb).checkProbes(vMsh, 0); }, Burg); std::cout << "Directory and files created.\n";

        ///// Temporal Parameters /////
        double dt = data["timeStep"].asDouble(), endTime = data["endTime"].asDouble(), tolTemporal = data["tolTemporal"].asDouble(), t{}; 
        std::visit([&](auto& vMsh) { Dsc.evaluateTimeStability(Mat, vMsh, dt, data["C1"].asDouble()); }, Burg); size_t iMax = endTime / dt;

        ///// Medic /////
        std::cout << "Initializing medic ...\n";
        bool bMdc = data["medicOn"].asBool(); std::optional<MedicSpectral> Mdc; if (bMdc) { std::visit([&](auto& vMsh) { Mdc.emplace(vMsh, *Prb, dt); }, Burg); }

        ///// Temporal Loop /////
        std::cout << "Processing ...\n";
        std::visit([&](auto& vMsh) {
                for (size_t i = 1; i < iMax; i++) { 
                // Control
                t = i * dt; vMsh.ouHat = vMsh.uHat; vMsh.oE = vMsh.E;

                // Solver
                if (!Sol.solveModel(vMsh.uHat, Mat.Re, t, dt, vMsh)) { throw std::logic_error("Simulation diverges @ t = " + std::to_string(t)); };

                // Energy Balance
                Dsc.calculateEnergy(vMsh.E, vMsh.uHat);

                // Probe & Medic
                (*Prb).checkProbes(vMsh, t);
                if (Mdc.has_value()) { (*Mdc).checkDiagnostics(Mat, vMsh); }

                // Convergence
                if (calcErr(vMsh.uHat, vMsh.ouHat) / dt < tolTemporal) { std::cout << "\nSteady-state achieved @ t = " << t << " seconds."; break; }

                std::cout << "\r" << double(100 * static_cast<double>(i) / iMax) << " %";
                } std::cout << "\n"; }, Burg);

        // End
        std::cout << "Files saved to: " << (*Prb).fPath << "\n";
    }
}

void runSpectral(const Json::Value& data, std::string configName) {
    ///// Simulation /////
    runSolverLoop(data, configName);
}
