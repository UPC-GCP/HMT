// Imports
#include <cmath>
#include <cstddef>
#include <iostream>
#include <json/json.h>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>

// Self-Imports
#include "o01_Material.h"
#include "o02_MeshSpectral.h"
#include "o03_DiscretizerSpectral.h"
#include "o04_SolverSpectral.h"
#include "o05_ProbeSpectral.h"
#include "o06_MedicSpectral.h"

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
        std::unique_ptr<MeshBurgers> Burg = nullptr;
        if (Sol.algorithm == solSpectral::DNS) { std::cout << "MeshBurg - creation\n"; Burg = std::make_unique<MeshBurgers>(); Msh.generateMeshBurgers(*Burg, data["N"].asInt()); }
        else if (Sol.algorithm == solSpectral::LES) { std::cout << "MeshLES - creation\n"; Burg = std::make_unique<MeshLES>(); Msh.generateMeshBurgers(static_cast<MeshLES&>(*Burg), data["N"].asInt(), data["Cs"].asDouble(), data["slope"].asDouble()); }
        std::cout << "Burgers mesh truncated at " << (*Burg).totNodes << " modes.\n";

        std::cout << "Mesh type: " << typeid(*Burg).name() << "\n";

        ///// Discretizer /////
        std::cout << "Initializing discretizer ...\n";
        DiscretizerSpectral Dsc; Dsc.initializeBurgers(*Burg, Mat.Re); std::cout << "Model initialized.\n";

        ///// Probe /////
        std::cout << "Initializing probe ...\n";
        ProbeSpectral Prb(*Burg, data["probes"], configName); Prb.checkProbes(*Burg, 0); std::cout << "Directory and files created.\n"; // Probe only stores E, add uHat later

        ///// Temporal Parameters /////
        double dt = data["timeStep"].asDouble(), endTime = data["endTime"].asDouble(), tolTemporal = data["tolTemporal"].asDouble(), t{}; size_t iMax = endTime / dt;

        ///// Medic /////
        std::cout << "Initializing medic ...\n";
        bool bMdc = data["medicOn"].asBool(); std::optional<MedicSpectral> Mdc; if (bMdc) { Mdc.emplace(*Burg, Prb, dt); }

        ///// Temporal Loop /////
        std::cout << "Processing ...\n";
        for (size_t i = 1; i < iMax; i++) {
            // Control
            t = i * dt; (*Burg).ouHat = (*Burg).uHat; (*Burg).oE = (*Burg).E;

            // Solver
            if (!Sol.solveModel((*Burg).uHat, Mat.Re, t, dt, *Burg)) { throw std::logic_error("Simulation diverges @ t = " + std::to_string(t)); };

            // Energy Balance
            Dsc.calculateEnergy((*Burg).E, (*Burg).uHat);

            // Probe & Medic
            Prb.checkProbes(*Burg, t);
            if (Mdc.has_value()) { Mdc->checkDiagnostics(Mat, *Burg); }

            // Convergence
            if (calcErr((*Burg).uHat, (*Burg).ouHat) / dt < tolTemporal) { std::cout << "\nSteady-state achieved @ t = " << t << " seconds."; break; }

            std::cout << "\r" << double(100 * static_cast<double>(i) / iMax) << " %";
        } std::cout << "\n";

        // End
        std::cout << "Files saved to: " << Prb.fPath << "\n";
    }
}

void runSpectral(const Json::Value& data, std::string configName) {
    ///// Simulation /////
    runSolverLoop(data, configName);
}
