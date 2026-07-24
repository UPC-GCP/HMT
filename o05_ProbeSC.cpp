// Imports
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>
#include <json/json.h>
#include <fstream>
#include <ctime>
#include <filesystem>

#include "o05_ProbeSC.h"

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
    pMap tempMap{}; pBug tempBug{}; pFld tempFld{};
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
            probePoint.xPos.push_back(std::lower_bound(Msh.p.Nodes[0].begin(), Msh.p.Nodes[0].end(), probes[i]["x0"][0].asDouble()) - Msh.p.Nodes[0].begin());
            probePoint.yPos.push_back(std::lower_bound(Msh.p.Nodes[1].begin(), Msh.p.Nodes[1].end(), probes[i]["x0"][1].asDouble()) - Msh.p.Nodes[1].begin());

        } else if (probes[i]["type"].asString() == "Map"){

            ///// Energy  /////
            if (!Msh.T.Phi.empty()){
                // Create File
                tempString = "Probe_" + std::to_string(probeMap.size() + 1) + "_Map.csv";
                tempMap.file = createFile(newPath / tempString);

                // Time
                tempMap.t = {probes[i]["t"][0].asDouble(), probes[i]["t"][1].asDouble()};
                tempMap.nWrite = probes[i]["nWrite"].isNull() ? 1 : probes[i]["nWrite"].asInt();

                // Position
                tempMap.xPos = {static_cast<size_t>(std::lower_bound(Msh.T.Nodes[0].begin(), Msh.T.Nodes[0].end(), probes[i]["x0"][0].asDouble()) - Msh.T.Nodes[0].begin()), static_cast<size_t>(std::lower_bound(Msh.T.Nodes[0].begin(), Msh.T.Nodes[0].end(), probes[i]["x1"][0].asDouble()) - Msh.T.Nodes[0].begin())};
                tempMap.yPos = {static_cast<size_t>(std::lower_bound(Msh.T.Nodes[1].begin(), Msh.T.Nodes[1].end(), probes[i]["x0"][1].asDouble()) - Msh.T.Nodes[1].begin()), static_cast<size_t>(std::lower_bound(Msh.T.Nodes[1].begin(), Msh.T.Nodes[1].end(), probes[i]["x1"][1].asDouble()) - Msh.T.Nodes[1].begin())};

                // Header
                for (int j = tempMap.xPos[0]; j < tempMap.xPos[1]; j++){
                    for (int k = tempMap.yPos[0]; k < tempMap.yPos[1]; k++){
                        tempMap.file << "," << Msh.T.Nodes[0][j] << " " << Msh.T.Nodes[1][k];
                    }
                } tempMap.file << "\n";

                // Control
                probeMap.push_back(std::move(tempMap));
                tempMap = {};
            }

            ///// Pressure /////
            // Create File
            tempString = "Probe_" + std::to_string(probeMap.size() + 1) + "_Map.csv";
            tempMap.file = createFile(newPath / tempString);

            // Time
            tempMap.t = {probes[i]["t"][0].asDouble(), probes[i]["t"][1].asDouble()};
            tempMap.nWrite = probes[i]["nWrite"].isNull() ? 1 : probes[i]["nWrite"].asInt();

            // Position
            tempMap.xPos = {static_cast<size_t>(std::lower_bound(Msh.p.Nodes[0].begin(), Msh.p.Nodes[0].end(), probes[i]["x0"][0].asDouble()) - Msh.p.Nodes[0].begin()), static_cast<size_t>(std::lower_bound(Msh.p.Nodes[0].begin(), Msh.p.Nodes[0].end(), probes[i]["x1"][0].asDouble()) - Msh.p.Nodes[0].begin())};
            tempMap.yPos = {static_cast<size_t>(std::lower_bound(Msh.p.Nodes[1].begin(), Msh.p.Nodes[1].end(), probes[i]["x0"][1].asDouble()) - Msh.p.Nodes[1].begin()), static_cast<size_t>(std::lower_bound(Msh.p.Nodes[1].begin(), Msh.p.Nodes[1].end(), probes[i]["x1"][1].asDouble()) - Msh.p.Nodes[1].begin())};
            
            // Header
            for (int j = tempMap.xPos[0]; j < tempMap.xPos[1]; j++){
                for (int k = tempMap.yPos[0]; k < tempMap.yPos[1]; k++){
                    tempMap.file << "," << Msh.p.Nodes[0][j] << " " << Msh.p.Nodes[1][k];
                }
            } tempMap.file << "\n";

            // Control
            probeMap.push_back(std::move(tempMap));
            tempMap = {};
        
        } else if (probes[i]["type"].asString() == "Field"){

            // Create File - uMesh
            tempString = "Probe_" + std::to_string(probeMap.size() + 1) + "u_Field.csv";
            tempFld.file = createFile(newPath / tempString);

            // Time - uMesh
            tempFld.t = {probes[i]["t"][0].asDouble(), probes[i]["t"][1].asDouble()};
            tempFld.nWrite = probes[i]["nWrite"].isNull() ? 1 : probes[i]["nWrite"].asInt();

            // Position - uMesh
            tempFld.xPos = {static_cast<size_t>(std::lower_bound(Msh.u.Nodes[0].begin(), Msh.u.Nodes[0].end(), probes[i]["x0"][0].asDouble()) - Msh.u.Nodes[0].begin()), static_cast<size_t>(std::lower_bound(Msh.u.Nodes[0].begin(), Msh.u.Nodes[0].end(), probes[i]["x1"][0].asDouble()) - Msh.u.Nodes[0].begin())};
            tempFld.yPos = {static_cast<size_t>(std::lower_bound(Msh.u.Nodes[1].begin(), Msh.u.Nodes[1].end(), probes[i]["x0"][1].asDouble()) - Msh.u.Nodes[1].begin()), static_cast<size_t>(std::lower_bound(Msh.u.Nodes[1].begin(), Msh.u.Nodes[1].end(), probes[i]["x1"][1].asDouble()) - Msh.u.Nodes[1].begin())};

            // Header - uMesh
            for (int j = tempFld.xPos[0]; j < tempFld.xPos[1]; j++){
                for (int k = tempFld.yPos[0]; k < tempFld.yPos[1]; k++){
                    tempFld.file << "," << Msh.u.Nodes[0][j] << " " << Msh.u.Nodes[1][k];
                }
            } tempFld.file << "\n";

            // Control
            probeFld.push_back(std::move(tempFld));
            tempFld = {};

            // Create File - vMesh
            tempString = "Probe_" + std::to_string(probeMap.size() + 1) + "v_Field.csv";
            tempFld.file = createFile(newPath / tempString);

            // Time - vMesh
            tempFld.t = {probes[i]["t"][0].asDouble(), probes[i]["t"][1].asDouble()};
            tempFld.nWrite = probes[i]["nWrite"].isNull() ? 1 : probes[i]["nWrite"].asInt();

            // Position - vMesh
            tempFld.xPos = {static_cast<size_t>(std::lower_bound(Msh.v.Nodes[0].begin(), Msh.v.Nodes[0].end(), probes[i]["x0"][0].asDouble()) - Msh.v.Nodes[0].begin()), static_cast<size_t>(std::lower_bound(Msh.v.Nodes[0].begin(), Msh.v.Nodes[0].end(), probes[i]["x1"][0].asDouble()) - Msh.v.Nodes[0].begin())};
            tempFld.yPos = {static_cast<size_t>(std::lower_bound(Msh.v.Nodes[1].begin(), Msh.v.Nodes[1].end(), probes[i]["x0"][1].asDouble()) - Msh.v.Nodes[1].begin()), static_cast<size_t>(std::lower_bound(Msh.v.Nodes[1].begin(), Msh.v.Nodes[1].end(), probes[i]["x1"][1].asDouble()) - Msh.v.Nodes[1].begin())};

            // Header - vMesh
            for (int j = tempFld.xPos[0]; j < tempFld.xPos[1]; j++){
                for (int k = tempFld.yPos[0]; k < tempFld.yPos[1]; k++){
                    tempFld.file << "," << Msh.v.Nodes[0][j] << " " << Msh.v.Nodes[1][k];
                }
            } tempFld.file << "\n";

            // Control
            probeFld.push_back(std::move(tempFld));
            tempFld = {};

        } else if (probes[i]["type"].asString() == "Debug"){

            // Create File
            tempString = "Probe_" + std::to_string(probeMap.size() + probeBug.size() + 1) + "_Bug.csv";
            tempBug.file = createFile(newPath / tempString);

            // Time
            tempBug.t = {probes[i]["t"][0].asDouble(), probes[i]["t"][1].asDouble()};

            // Position
            tempBug.xPos = {static_cast<size_t>(std::lower_bound(Msh.p.Nodes[0].begin(), Msh.p.Nodes[0].end(), probes[i]["x0"][0].asDouble()) - Msh.p.Nodes[0].begin()), static_cast<size_t>(std::lower_bound(Msh.p.Nodes[0].begin(), Msh.p.Nodes[0].end(), probes[i]["x1"][0].asDouble()) - Msh.p.Nodes[0].begin())};
            tempBug.yPos = {static_cast<size_t>(std::lower_bound(Msh.p.Nodes[1].begin(), Msh.p.Nodes[1].end(), probes[i]["x0"][1].asDouble()) - Msh.p.Nodes[1].begin()), static_cast<size_t>(std::lower_bound(Msh.p.Nodes[1].begin(), Msh.p.Nodes[1].end(), probes[i]["x1"][1].asDouble()) - Msh.p.Nodes[1].begin())};
            
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
                probePoint.file << Msh.p.Phi[probePoint.xPos[i]][probePoint.yPos[i]];
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

    // Save Values (even index = T-map, odd index = p-map)
    for (size_t i = 0; i < probeMap.size(); i++){
        if (!bSave[i]){continue;}
        if (probeMap[i].nCount++ % probeMap[i].nWrite != 0){continue;}

        probeMap[i].file << t;
        if (i % 2 == 0 && !Msh.T.Phi.empty()){
            for (size_t j = probeMap[i].xPos[0]; j < probeMap[i].xPos[1]; j++){
                for (size_t k = probeMap[i].yPos[0]; k < probeMap[i].yPos[1]; k++){
                    probeMap[i].file << "," << Msh.T.Phi[j][k];
                }
            }
        } else {
            for (size_t j = probeMap[i].xPos[0]; j < probeMap[i].xPos[1]; j++){
                for (size_t k = probeMap[i].yPos[0]; k < probeMap[i].yPos[1]; k++){
                    probeMap[i].file << "," << Msh.p.Phi[j][k];
                }
            }
        }
        probeMap[i].file << "\n";
    }


    // Field Probe
    bSave = {}; bSave.resize(probeFld.size(), false);
    for (size_t i = 0; i < probeFld.size(); i++){
        if (t >= probeFld[i].t[0] && t <= probeFld[i].t[1]){
            bSave[i] = true;
        }
    }

    // Save uMesh
    if (probeFld.size() > 0 && bSave[0]){
        if (probeFld[0].nCount++ % probeFld[0].nWrite == 0){
            probeFld[0].file << t;
            for (size_t i = probeFld[0].xPos[0]; i < probeFld[0].xPos[1]; i++){
                for (size_t j = probeFld[0].yPos[0]; j < probeFld[0].yPos[1]; j++){
                    probeFld[0].file << "," << Msh.u.Phi[i][j];
                }
            } probeFld[0].file << "\n";
        }
    }

    // Save vMesh
    if (probeFld.size() > 1 && bSave[1]){
        if (probeFld[1].nCount++ % probeFld[1].nWrite == 0){
            probeFld[1].file << t;
            for (size_t i = probeFld[1].xPos[0]; i < probeFld[1].xPos[1]; i++){
                for (size_t j = probeFld[1].yPos[0]; j < probeFld[1].yPos[1]; j++){
                    probeFld[1].file << "," << Msh.v.Phi[i][j];
                }
            } probeFld[1].file << "\n";
        }
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

    // Field Probes
    if (!probeFld.empty()){
        for (int i = 0; i < probeFld.size(); i++){
            probeFld[i].file.close();
        }
    }

    // Bug Probes
    if (!probeBug.empty()){
	    for (int i = 0; i < probeBug.size(); i++){
		    probeBug[i].file.close();
	    }
    }

}

