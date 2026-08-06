// Imports
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <json/json.h>
#include <fstream>
#include <ctime>
#include <filesystem>

#include "o05_Probe.h"
#include "o02_Mesh.h"
#include "o04_Solver.h"

std::string createFolder(std::string tempScheme, std::string spatScheme, std::string fName, std::string& dirName){

    // Control
    if (tempScheme.find('-') != std::string::npos){tempScheme.erase(tempScheme.find('-'), 1);}

    // Directory
    std::filesystem::path pBase = std::filesystem::current_path();
    pBase /= "ioRes";

    // Timestamp
    time_t timeStamp = std::time(nullptr);
    struct tm datetime = *localtime(&timeStamp);
    char oName[35]; strftime(oName, sizeof(oName), "%Y%m%d%H%M%S_", &datetime);
    
    // Folder Name
    int iPos = fName.find(".json"); dirName = oName + fName.substr(0, iPos) + "_" + tempScheme + "_" + spatScheme;
    pBase /= dirName;

    // Create Folder
    std::filesystem::create_directories(pBase);

    return pBase.string();

}

std::ofstream createFile(std::filesystem::path fName){

    // Open File 
    std::ofstream file(fName);
    if (!file.is_open()){
        std::cerr << "Failed to open file. \n";
    }

    // File Header
    file << "Time";

    return file;

}

Probe::Probe(Mesh Msh, Json::Value probes, std::string tempScheme, std::string spatScheme, std::string fName){
    
    // Create Folder
    std::string tempString{};
    std::filesystem::path newPath(fName);
    pathBase = createFolder(tempScheme, spatScheme, newPath.filename().string(), dirName);
    newPath = pathBase;
    
    // Add Probes and Create Files
    pMap tempMap{}; pBug tempBug{}; pFld tempFld{}; pBnd tempBnd{};
    for (Json::Value::ArrayIndex i = 0; i < probes.size(); i++){

        if (probes[i]["type"].asString() == "Point"){

            // Create File
            if (probePoint.bFile){
                probePoint.file = createFile(newPath / "Probe_0_Point.csv"); probePoint.bFile = false;
            }

            // Header
            probePoint.file << "," << probes[i]["x0"][0].asDouble() << " " << probes[i]["x0"][1].asDouble();
            
            // Time
            probePoint.t0.push_back(probes[i]["t"][0].asDouble()); probePoint.t1.push_back(probes[i]["t"][1].asDouble());

            // Position
            probePoint.xPos.push_back(std::lower_bound(Msh.Nodes[0].begin(), Msh.Nodes[0].end(), probes[i]["x0"][0].asDouble()) - Msh.Nodes[0].begin());
            probePoint.yPos.push_back(std::lower_bound(Msh.Nodes[1].begin(), Msh.Nodes[1].end(), probes[i]["x0"][1].asDouble()) - Msh.Nodes[1].begin());

        } else if (probes[i]["type"].asString() == "Map"){
            
            // Create File
            tempString = "Probe_" + std::to_string(probeMap.size() + 1) + "_Map.csv";
            tempMap.file = createFile(newPath / tempString);

            // Time
            tempMap.t = {probes[i]["t"][0].asDouble(), probes[i]["t"][1].asDouble()};

            // Position
            tempMap.xPos = {static_cast<size_t>(std::lower_bound(Msh.Nodes[0].begin(), Msh.Nodes[0].end(), probes[i]["x0"][0].asDouble()) - Msh.Nodes[0].begin()), static_cast<size_t>(std::lower_bound(Msh.Nodes[0].begin(), Msh.Nodes[0].end(), probes[i]["x1"][0].asDouble()) - Msh.Nodes[0].begin())};
            tempMap.yPos = {static_cast<size_t>(std::lower_bound(Msh.Nodes[1].begin(), Msh.Nodes[1].end(), probes[i]["x0"][1].asDouble()) - Msh.Nodes[1].begin()), static_cast<size_t>(std::lower_bound(Msh.Nodes[1].begin(), Msh.Nodes[1].end(), probes[i]["x1"][1].asDouble()) - Msh.Nodes[1].begin())};
            
            // Header
            for (int j = tempMap.xPos[0]; j < tempMap.xPos[1]; j++){
                for (int k = tempMap.yPos[0]; k < tempMap.yPos[1]; k++){
                    tempMap.file << "," << Msh.Nodes[0][j] << " " << Msh.Nodes[1][k];
                }
            } tempMap.file << "\n";

            // Control
            probeMap.push_back(std::move(tempMap));
            tempMap = {};
        
        } else if (probes[i]["type"].asString() == "Field"){

            // Not using this one for now since I am using a static field and I can calculate it again later but will fully add it later

            /* // Create File */
            /* tempString = "Probe_" + std::to_string(probeMap.size() + probeFld.size() + 1) + "_Field.csv"; */
            /* tempFld.file = createFile(newPath / tempString); */

            /* // Time */
            /* tempFld.t = {probes[i]["t"][0].asDouble(), probes[i]["t"][1].asDouble()}; */

            /* // Position */

            /* // Header */
            /* for (int j = tempFld.xPos[0]; j <= tempFld.xPos[1]; j++){ */
            /*     for (int k = tempFld.yPos[0]; k <= tempFld.yPos[1]; k++){ */
            /*         tempFld.file << "," << Msh.Faces[0][j] << " " << Msh.Faces[1][k]; */
            /*     } */
            /* } */

            /* // Control */
            /* probeFld.push_back(std::move(tempFld)); */
            /* tempFld = {}; */
            

        } else if (probes[i]["type"].asString() == "Boundary"){

            // Create File
            tempString = "Probe_" + std::to_string(probeMap.size() + probeBnd.size() + 1) + "_Boundary.csv";
            tempBnd.file = createFile(newPath / tempString);

            // Data
            tempBnd.side = probes[i]["side"].asInt(); double epsFind = 1e-5;
            tempBnd.t = {probes[i]["t"][0].asDouble(), probes[i]["t"][1].asDouble()};

            // Position
            tempBnd.xPos = {probes[i]["x0"][0].asDouble(), probes[i]["x1"][0].asDouble()};
            tempBnd.yPos = {probes[i]["x0"][1].asDouble(), probes[i]["x1"][1].asDouble()};


            // Boundaries
            for (size_t k = 0; k < Msh.boundaryConditions.size(); k++){

                // Headers
                if (Msh.boundaryConditions[k].x0[0] >= tempBnd.xPos[0] && Msh.boundaryConditions[k].x1[0] <= tempBnd.xPos[1] && Msh.boundaryConditions[k].x0[1] >= tempBnd.yPos[0] && Msh.boundaryConditions[k].x1[1] <= tempBnd.yPos[1]){

                    // Control
                    tempBnd.iBC.push_back(k); int l{}, m{};

                    if (tempBnd.xPos[0] == tempBnd.xPos[1]){
                        // xBoundary
                        if (Msh.boundaryConditions[k].side == 0){l = 0;} else if (Msh.boundaryConditions[k].side == 1){l = Msh.Faces[0].size()-1;}
                        /* tempBnd.iPos = {static_cast<size_t>(std::lower_bound(Msh.Nodes[1].begin(), Msh.Nodes[1].end(), Msh.boundaryConditions[k].x0[1] - epsFind) - Msh.Nodes[1].begin()), static_cast<size_t>(std::lower_bound(Msh.Nodes[1].begin(), Msh.Nodes[1].end(), Msh.boundaryConditions[k].x1[1] - epsFind) - Msh.Nodes[1].begin())}; */
                        tempBnd.iPos.push_back(static_cast<size_t>(std::lower_bound(Msh.Nodes[1].begin(), Msh.Nodes[1].end(), Msh.boundaryConditions[k].x0[1] - epsFind) - Msh.Nodes[1].begin()));
                        tempBnd.iPos.push_back(static_cast<size_t>(std::lower_bound(Msh.Nodes[1].begin(), Msh.Nodes[1].end(), Msh.boundaryConditions[k].x1[1] - epsFind) - Msh.Nodes[1].begin()));
                        for (int m = tempBnd.iPos[tempBnd.iPos.size()-2]; m < tempBnd.iPos.back(); m++){tempBnd.file << "," << Msh.Faces[0][l] << " " << Msh.Nodes[1][m];}

                    } else if (tempBnd.yPos[0] == tempBnd.yPos[1]){
                        // yBoundary
                        if (Msh.boundaryConditions[k].side == 0){m = 0;} else if (Msh.boundaryConditions[k].side == 1){m = Msh.Faces[1].size()-1;}
                        /* tempBnd.iPos = {static_cast<size_t>(std::lower_bound(Msh.Nodes[0].begin(), Msh.Nodes[0].end(), Msh.boundaryConditions[k].x0[0] - epsFind) - Msh.Nodes[0].begin()), static_cast<size_t>(std::lower_bound(Msh.Nodes[0].begin(), Msh.Nodes[0].end(), Msh.boundaryConditions[k].x1[0] - epsFind) - Msh.Nodes[0].begin())}; */
                        tempBnd.iPos.push_back(static_cast<size_t>(std::lower_bound(Msh.Nodes[0].begin(), Msh.Nodes[0].end(), Msh.boundaryConditions[k].x0[0] - epsFind) - Msh.Nodes[0].begin()));
                        tempBnd.iPos.push_back(static_cast<size_t>(std::lower_bound(Msh.Nodes[0].begin(), Msh.Nodes[0].end(), Msh.boundaryConditions[k].x1[0] - epsFind) - Msh.Nodes[0].begin()));
                        for (int l = tempBnd.iPos[tempBnd.iPos.size()-2]; l < tempBnd.iPos.back(); l++){tempBnd.file << "," << Msh.Nodes[0][l] << " " << Msh.Faces[1][m];}
                    } else {std::cerr << "Probe range not specified correctly\n";}

                }


            } tempBnd.file << "\n";

            // Control
            probeBnd.push_back(std::move(tempBnd));
            tempBnd = {};

        } else if (probes[i]["type"].asString() == "Debug"){

            // Create File
            tempString = "Probe_" + std::to_string(probeMap.size() + probeBnd.size() + probeBug.size() + 1) + "_Bug.csv";
            tempBug.file = createFile(newPath / tempString);

            // Time
            tempBug.t = {probes[i]["t"][0].asDouble(), probes[i]["t"][1].asDouble()};

            // Position
            tempBug.xPos = {static_cast<size_t>(std::lower_bound(Msh.Nodes[0].begin(), Msh.Nodes[0].end(), probes[i]["x0"][0].asDouble()) - Msh.Nodes[0].begin()), static_cast<size_t>(std::lower_bound(Msh.Nodes[0].begin(), Msh.Nodes[0].end(), probes[i]["x1"][0].asDouble()) - Msh.Nodes[0].begin())};
            tempBug.yPos = {static_cast<size_t>(std::lower_bound(Msh.Nodes[1].begin(), Msh.Nodes[1].end(), probes[i]["x0"][1].asDouble()) - Msh.Nodes[1].begin()), static_cast<size_t>(std::lower_bound(Msh.Nodes[1].begin(), Msh.Nodes[1].end(), probes[i]["x1"][1].asDouble()) - Msh.Nodes[1].begin())};
            
            // Header
            tempBug.file << ",lastIter,lastRes\n";

            // Control
            probeBug.push_back(std::move(tempBug));
            tempBug = {};

        } else {
            std::cerr << "Probe type not recognized. \n";
        }
        
    }

    // Control
    probePoint.file << "\n";

}

void Probe::checkProbes(Mesh Msh, Solver* Sol, double t){

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
                probePoint.file << Msh.vPhi[probePoint.xPos[i]][probePoint.yPos[i]];
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
        for (size_t j = probeMap[i].xPos[0]; j < probeMap[i].xPos[1]; j++){
            for (size_t k = probeMap[i].yPos[0]; k < probeMap[i].yPos[1]; k++){
                probeMap[i].file << "," << Msh.vPhi[j][k];
            }
        } probeMap[i].file << "\n";
    }

    // Boundary Probe
    bSave = {}; bSave.resize(probeBnd.size(), false);
    for (size_t i = 0; i < probeBnd.size(); i++){
        if (t >= probeBnd[i].t[0] && t <= probeBnd[i].t[1]){
            bSave[i] = true;
        }
    }

    // Save Values
    for (size_t i = 0; i < probeBnd.size(); i++){
        if (!bSave[i]){continue;}

        probeBnd[i].file << t;
        for (size_t j = 0; j < probeBnd[i].iBC.size(); j++){
            for (size_t k = probeBnd[i].iPos[2*j]; k < probeBnd[i].iPos[2*j+1]; k++){
                probeBnd[i].file << "," << Msh.boundaryConditions[probeBnd[i].iBC[j]].Phi[k];
            }
        } probeBnd[i].file << "\n";
    }

    // Bug Probe
    bSave = {}; bSave.resize(probeBug.size(), false);
    for (size_t i = 0; i < probeBug.size(); i++){
        if (t >= probeBug[i].t[0] && t <= probeBug[i].t[1]){
            bSave[i] = true;
        }
    }

    // Save Values
    for (size_t i = 0; i < probeBug.size(); i++){
        if (!bSave[i]){continue;}
        probeBug[i].file << t << "," << Sol->lastIter << "," << Sol->lastRes << "\n";
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

    // Boundary Probes
    if (!probeBnd.empty()){
        for (int i = 0; i < probeBnd.size(); i++){
            probeBnd[i].file.close();
        }
    }

    // Bug Probes
    if (!probeBug.empty()){
	    for (int i = 0; i < probeBug.size(); i++){
		    probeBug[i].file.close();
	    }
    }

}
