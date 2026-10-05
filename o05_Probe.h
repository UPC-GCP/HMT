#ifndef PROBE_H_
#define PROBE_H_

// Imports
#include <ctime>
#include <array>
#include <vector>
#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <filesystem>
#include <json/json.h>

template <size_t Dim> struct probeBase {
    std::ofstream file{};
    size_t nWrite{1}, nCount{};
};

template <size_t Dim> struct probePoint : probeBase<Dim> {
    // pPoint Definition -- Single object for all probes
    // Targets specific coordinates and stores the value within a range of time
    // Each point can be stored in different intervals -- t = vector<double>
    // Each point can have N dimensions -- i0 = array<Dim> 
    // Can have multiple poins -- i0 = vector<array<Dim>>
    std::vector<std::array<size_t, Dim>> i0{};
    std::vector<std::array<double, 2>> t{};
};

template <size_t Dim> struct probeRange : probeBase<Dim> {
    // pRange Definition - One object for each probe
    // Targets range of coordinates and stores values within a range of time
    // Each range can have N dimensions -- i0, i1 = array<Dim>
    // Each range is defined with an interval -- t = array<double, 2>
    std::array<size_t, Dim> i0{}, i1{};
    std::array<double, 2> t{};
};

template <size_t Dim> class Probe {
private:

public:
    // Variables
    std::string pathBase{}, dirName{}, uName{};
    probePoint<Dim> pPoint{}; probeRange<Dim> pRange{}; // Probe obects -- not sure if leave here or move to specific files
    /* pPoint probePoint{}; */
    /* std::vector<pMap> probeMap{}; */
    /* std::vector<pFld> probeFld{}; */
    /* std::vector<pBug> probeBug{}; */

    // Constructor
    Probe(const Json::Value& probes, std::string fName);
    /* Probe(Mesh Msh, Json::Value probes, std::string tempScheme, std::string spatScheme, std::string fName); */
    
    // Destructor
    ~Probe();

    // Headers
    /* void checkProbes(Mesh Msh, Solver* Sol, double t=0); */
};

// Functions
inline std::string createFolder(std::string fName, std::string& dirName) {
    // Timestamp
    time_t timeStamp = std::time(nullptr);
    struct tm datetime = *localtime(&timeStamp);
    char oName[35]; strftime(oName, sizeof(oName), "%Y%m%d%H%M%S_", &datetime);
    
    // Folder Name
    size_t iPos = fName.find(".json"); dirName = oName + fName.substr(0, iPos);
    std::filesystem::path pBase = std::filesystem::current_path(); pBase /= dirName;

    // Create Folder
    std::filesystem::create_directories(pBase);

    return pBase.string();
}

inline std::ofstream createFile(std::filesystem::path fName){
    // Create File
    std::ofstream file(fName); if (!file.is_open()){ throw std::runtime_error("File could not be created: " + fName.string()); }
    file << "Time";
    
    return file;
}

// Constuctor
template <size_t Dim> Probe<Dim>::Probe(const Json::Value& probes, std::string fName) {
    // Create Folder
    std::filesystem::path newPath(fName);
    pathBase = createFolder(newPath.filename().string(), dirName);

    // CHECK INITIALIZER LIST TO PROPERLY CONFIGURE THIS WITHIN EACH ONE AND LEAVE CREATE FOLDER HERE
    /* // Add Probes -- This should probably be done within each one */
    /* for (Json::Value::ArrayIndex k = 0; k < probes.size(); k++) { */

    /*     if (probes[k]["type"].asString() == "Point") { */

    /*     } else if (probes[k]["type"].asString() == "Map") { */

    /*     } else if (probes[k]["type"].asString() == "Field") { */

    /*     } */

    /* } */
}

// Destructor
template <size_t Dim> Probe<Dim>::~Probe() {
    // Point Probe

}


#endif
