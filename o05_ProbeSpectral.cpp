// Self-Imports
#include "o05_ProbeSpectral.h"
#include "o05_Probe.h"
#include <string>

ProbeSpectral::~ProbeSpectral() {
    // pModal

}

ProbeSpectral::ProbeSpectral(const Json::Value& probes, std::string fName) : Probe(probes, fName) {
    // Control
    std::string tempString{}; probeRange<1> pTemp{};
    
    // Add Probes
    for (Json::Value::ArrayIndex k = 0; k < probes.size(); k++) {
        
        // Also test what happens when nWrite is empty/null

        if (probes[k]["type"].asString() == "Modal") {

            // Create File
            tempString = "Probe_" + std::to_string(iProbes + 1) + "_Modal.csv";
            /* pTemp.file = createFile(newPath / tempString); */

            // newPath not recognized, need to keep it as a class variable and check how to initialize it

            // Time
            pTemp.t = {probes[k]["t"][0].asDouble(), probes[k]["t"][1].asDouble()};


            /* ///// Energy  ///// */
            /* if (!Msh.T.Phi.empty()){ */
            /*     // Create File */
            /*     tempString = "Probe_" + std::to_string(probeMap.size() + 1) + "_Map.csv"; */
            /*     tempMap.file = createFile(newPath / tempString); */

            /*     // Time */
            /*     tempMap.t = {probes[i]["t"][0].asDouble(), probes[i]["t"][1].asDouble()}; */
            /*     tempMap.nWrite = probes[i]["nWrite"].isNull() ? 1 : probes[i]["nWrite"].asInt(); */

            /*     // Position */
            /*     tempMap.xPos = {static_cast<size_t>(std::lower_bound(Msh.T.Nodes[0].begin(), Msh.T.Nodes[0].end(), probes[i]["x0"][0].asDouble()) - Msh.T.Nodes[0].begin()), static_cast<size_t>(std::lower_bound(Msh.T.Nodes[0].begin(), Msh.T.Nodes[0].end(), probes[i]["x1"][0].asDouble()) - Msh.T.Nodes[0].begin())}; */
            /*     tempMap.yPos = {static_cast<size_t>(std::lower_bound(Msh.T.Nodes[1].begin(), Msh.T.Nodes[1].end(), probes[i]["x0"][1].asDouble()) - Msh.T.Nodes[1].begin()), static_cast<size_t>(std::lower_bound(Msh.T.Nodes[1].begin(), Msh.T.Nodes[1].end(), probes[i]["x1"][1].asDouble()) - Msh.T.Nodes[1].begin())}; */

            /*     // Header */
            /*     for (int j = tempMap.xPos[0]; j < tempMap.xPos[1]; j++){ */
            /*         for (int k = tempMap.yPos[0]; k < tempMap.yPos[1]; k++){ */
            /*             tempMap.file << "," << Msh.T.Nodes[0][j] << " " << Msh.T.Nodes[1][k]; */
            /*         } */
            /*     } tempMap.file << "\n"; */

            /*     // Control */
            /*     probeMap.push_back(std::move(tempMap)); */
            /*     tempMap = {}; */
            /* } */

            /* ///// Pressure ///// */
            /* // Create File */
            /* tempString = "Probe_" + std::to_string(probeMap.size() + 1) + "_Map.csv"; */
            /* tempMap.file = createFile(newPath / tempString); */

            /* // Time */
            /* tempMap.t = {probes[i]["t"][0].asDouble(), probes[i]["t"][1].asDouble()}; */
            /* tempMap.nWrite = probes[i]["nWrite"].isNull() ? 1 : probes[i]["nWrite"].asInt(); */

            /* // Position */
            /* tempMap.xPos = {static_cast<size_t>(std::lower_bound(Msh.p.Nodes[0].begin(), Msh.p.Nodes[0].end(), probes[i]["x0"][0].asDouble()) - Msh.p.Nodes[0].begin()), static_cast<size_t>(std::lower_bound(Msh.p.Nodes[0].begin(), Msh.p.Nodes[0].end(), probes[i]["x1"][0].asDouble()) - Msh.p.Nodes[0].begin())}; */
            /* tempMap.yPos = {static_cast<size_t>(std::lower_bound(Msh.p.Nodes[1].begin(), Msh.p.Nodes[1].end(), probes[i]["x0"][1].asDouble()) - Msh.p.Nodes[1].begin()), static_cast<size_t>(std::lower_bound(Msh.p.Nodes[1].begin(), Msh.p.Nodes[1].end(), probes[i]["x1"][1].asDouble()) - Msh.p.Nodes[1].begin())}; */
            
            /* // Header */
            /* for (int j = tempMap.xPos[0]; j < tempMap.xPos[1]; j++){ */
            /*     for (int k = tempMap.yPos[0]; k < tempMap.yPos[1]; k++){ */
            /*         tempMap.file << "," << Msh.p.Nodes[0][j] << " " << Msh.p.Nodes[1][k]; */
            /*     } */
            /* } tempMap.file << "\n"; */

            /* // Control */
            /* probeMap.push_back(std::move(tempMap)); */
            /* tempMap = {}; */

        }

    }


}
