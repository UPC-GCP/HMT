// Imports
#include <iostream>
#include <vector>
#include <json/json.h>
#include <fstream>
#include <ctime>
#include <filesystem>

#include "Probe.h"
#include "Mesh.h"

std::string createFolder(std::string scheme, std::string fName){

    // Directory
    std::filesystem::path dName = std::filesystem::current_path(); 
    std::string pBase = dName.string() + "\\TestData\\";

    // Timestamp
    time_t timeStamp = std::time(nullptr);
    struct tm datetime = *localtime(&timeStamp);
    char oName[35]; strftime(oName, sizeof(oName), "%Y%m%d%H%M%S_", &datetime);
    
    // Folder Name
    int iPos = fName.find(".json"); std::string dirName = oName + fName.substr(0, iPos) + "_" + scheme;
    pBase += "\\" + dirName + "\\";

    // Create Folder
    std::filesystem::create_directories(pBase);

    return pBase;

}

std::ofstream createFile(std::string fName){

    // Open File 
    std::ofstream file(fName);
    if (!file.is_open()){
        std::cerr << "Failed to open file. \n";
    }

    // File Header
    file << "Time";

    return file;

}

Probe::Probe(Mesh Msh, Json::Value probes, std::string scheme, std::string fName){

    // Create Folder
    pathBase = createFolder(scheme, fName.substr(12, fName.size()-1));
    
    // Add Probes and Create Files
    pMap tempMap{};
    for (Json::Value::ArrayIndex i = 0; i < probes.size(); i++){

        if (probes[i]["type"].asString() == "Point"){

            // Create File
            if (probePoint.bFile){
                probePoint.file = createFile(pathBase + "Probe_0_Point.csv"); probePoint.bFile = false;
            }

            // Header
            probePoint.file << ",[" << probes[i]["x0"][0].asDouble() << " " << probes[i]["x0"][1].asDouble() << "]";
            
            // Time
            probePoint.t0.push_back(probes[i]["t"][0].asDouble()); probePoint.t1.push_back(probes[i]["t"][1].asDouble());

            // Position
            probePoint.xPos.push_back(std::lower_bound(Msh.Nodes[0].begin(), Msh.Nodes[0].end(), probes[i]["x0"][0].asDouble()) - Msh.Nodes[0].begin());
            probePoint.yPos.push_back(std::lower_bound(Msh.Nodes[1].begin(), Msh.Nodes[1].end(), probes[i]["x0"][1].asDouble()) - Msh.Nodes[1].begin());

        } else if (probes[i]["type"].asString() == "Map"){
            
            // Create File
            tempMap.file = createFile(pathBase + "Probe_" + std::to_string(probeMap.size() + 1) + "_Map.csv");

            // Time
            tempMap.t = {probes[i]["t"][0].asDouble(), probes[i]["t"][1].asDouble()};

            // Position
            tempMap.xPos = {static_cast<size_t>(std::lower_bound(Msh.Nodes[0].begin(), Msh.Nodes[0].end(), probes[i]["x0"][0].asDouble()) - Msh.Nodes[0].begin()), static_cast<size_t>(std::lower_bound(Msh.Nodes[0].begin(), Msh.Nodes[0].end(), probes[i]["x1"][0].asDouble()) - Msh.Nodes[0].begin())};
            tempMap.yPos = {static_cast<size_t>(std::lower_bound(Msh.Nodes[1].begin(), Msh.Nodes[1].end(), probes[i]["x0"][1].asDouble()) - Msh.Nodes[1].begin()), static_cast<size_t>(std::lower_bound(Msh.Nodes[1].begin(), Msh.Nodes[1].end(), probes[i]["x1"][1].asDouble()) - Msh.Nodes[1].begin())};
            
            // Header
            for (int j = tempMap.xPos[0]; j <= tempMap.xPos[1]; j++){
                for (int k = tempMap.yPos[0]; k <= tempMap.yPos[1]; k++){
                    tempMap.file << ",[" << Msh.Nodes[0][j] << " " << Msh.Nodes[1][k] << "]";
                }
            } tempMap.file << "\n";

            // Control
            probeMap.push_back(std::move(tempMap));
            tempMap = {};

        } else {
            std::cerr << "Probe type not recognized. \n";
        }
        
    }

    // Control
    probePoint.file << "\n";

}

void Probe::checkProbes(Mesh Msh, double t){
    
    // Point Probe
    std::vector<bool> bSave{}; bSave.resize(probePoint.xPos.size(), false);
    for (size_t i = 0; i < probePoint.xPos.size(); i++){
        if (t >= probePoint.t0[i] && t <= probePoint.t1[i]){
            bSave[i] = true;
        }
    }

    // Save Values
    if (std::find(bSave.begin(), bSave.end(), true) != bSave.end()){
        probePoint.file << t;
        for (size_t i = 0; i < probePoint.xPos.size(); i++){
            probePoint.file << ",";
            if (t >= probePoint.t0[i] && t <= probePoint.t1[i]){
                probePoint.file << Msh.nT[probePoint.xPos[i]][probePoint.yPos[i]];
            }
        } probePoint.file << "\n";
    }

    // Map Probe
    bSave = {}; bSave.resize(probeMap.size(), false);
    for (size_t i = 0; i < probeMap.size(); i++){
        if (t >= probeMap[i].t[0] && t <= probeMap[i].t[1]){
            bSave[i] = true;
        }
    }
    
    // Save Values
    for (size_t i = 0; i < probeMap.size(); i++){
        if (!bSave[i]){continue;}

        probeMap[i].file << t;
        for (size_t j = probeMap[i].xPos[0]; j <= probeMap[i].xPos[1]; j++){
            for (size_t k = probeMap[i].yPos[0]; k <= probeMap[i].yPos[1]; k++){
                probeMap[i].file << "," << Msh.nT[j][k];
            }
        } probeMap[i].file << "\n";
    }

}

Probe::~Probe(){
    
    // Point Probe
    if (!probePoint.bFile){
        probePoint.file.close();
    }
    
    // Map Probes
    if (!probeMap.empty()){
        for (int i = 0; i < probeMap.size(); i++){
            probeMap[i].file.close();
        }
    }

}