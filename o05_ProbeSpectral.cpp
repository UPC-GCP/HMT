// Imports
#include <string>
#include <cstddef>

// Self-Imports
#include "o02_Mesh.h"
#include "o05_ProbeSpectral.h"

ProbeSpectral::~ProbeSpectral() {
    // pModal
    if (!pModal.empty()) { for (size_t k = 0; k < pModal.size(); k++) { pModal[k].file.close(); } }
}

ProbeSpectral::ProbeSpectral(const MeshBurgers& Burg, const Json::Value& probes, std::string fName) : Probe(probes, fName) {
    // Control
    std::string tempString{}; probeRange<1> pTemp{};
    
    // Add Probes
    for (Json::Value::ArrayIndex k = 0; k < probes.size(); k++) {
        if (probes[k]["type"].asString() == "Modal") { // Modal
            // Create File
            tempString = "Probe_" + std::to_string(iProbes + 1) + "_Modal.csv";
            /* pTemp.file = createFile(newPath / tempString); */

            // newPath not recognized, need to keep it as a class variable and check how to initialize it

            // Time
            pTemp.t = { probes[k]["t"][0].asDouble(), probes[k]["t"][1].asDouble() };
	        if (!probes[k]["nWrite"].isNull()) { pTemp.nWrite = probes[k]["nWrite"].asInt(); }

            // Position
            pTemp.i0 = { static_cast<size_t>(probes[k]["x0"][0].asInt()) };
            pTemp.i1 = { static_cast<size_t>(probes[k]["x1"][0].asInt()) };

            // Header
            runLoopMesh<1>(Burg.N, [&](std::array<size_t, 3> iX, size_t Ny, size_t Nx) {
                    pTemp.file << "," << iX[0];
                    }, pTemp.i0, pTemp.i1); pTemp.file << "\n";

            // Control
            pModal.push_back(std::move(pTemp));
            pTemp = {};
        }
    }
        // Also test what happens when nWrite is empty/null
}
