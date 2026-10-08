// Imports
#include <string>
#include <cstddef>
#include <iostream>

// Self-Imports
#include "o02_Mesh.h"
#include "o02_MeshSpectral.h"
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
            pTemp.file = createFile(fPath / tempString);

            // Time
            pTemp.t = { probes[k]["t"][0].asDouble(), probes[k]["t"][1].asDouble() };
	        if (!probes[k]["nWrite"].isNull()) { pTemp.nWrite = probes[k]["nWrite"].asInt(); }

            // Position
            pTemp.i0 = { static_cast<size_t>(probes[k]["x0"][0].asInt()) };
            pTemp.i1 = { static_cast<size_t>(probes[k]["x1"][0].asInt()) }; for (size_t& val : pTemp.i1) { val += 1; }

            // Header
            runLoopMesh<1>(Burg.N, [&](std::array<size_t, 3> iX, size_t Ny, size_t Nx) { pTemp.file << "," << iX[0]; }, pTemp.i0, pTemp.i1); pTemp.file << "\n";

            // Control
            pModal.push_back(std::move(pTemp));
            pTemp = {};
        }
    }
}

void ProbeSpectral::checkProbes(const MeshBurgers& Burg, double t) {
    // Modal
    for (size_t k = 0; k < pModal.size(); k++) {
        // Control
        if (t < pModal[k].t[0] || t > pModal[k].t[1]) { continue; }
        if (pModal[k].nCount++ % pModal[k].nWrite != 0) { continue; }

        // Write Data
        pModal[k].file << t; 
        runLoopMesh<1>(Burg.N, [&](std::array<size_t, 3> iX, size_t Ny, size_t Nx) { pModal[k].file << "," << Burg.E[iX[0]]; }, pModal[k].i0, pModal[k].i1); pModal[k].file << "\n";
    }

}
