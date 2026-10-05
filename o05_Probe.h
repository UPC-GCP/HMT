#ifndef PROBE_H_
#define PROBE_H_

// Imports
#include <array>
#include <cstddef>
#include <fstream>
#include <ctime>
#include <filesystem>
#include <stdexcept>

template <size_t Dim> struct pBase {
    std::ofstream file{};
    size_t nWrite{1}, nCount{};
    std::array<size_t, Dim> i0{}, i1{};
    std::array<double, 2> t{};
};

template <size_t Dim> class Probe {
private:

public:
    // Variables
    std::string pathBase{}, dirName{}, uName{};
    /* pPoint probePoint{}; */
    /* std::vector<pMap> probeMap{}; */
    /* std::vector<pFld> probeFld{}; */
    /* std::vector<pBug> probeBug{}; */

    // Constructor
    /* Probe(Mesh Msh, Json::Value probes, std::string tempScheme, std::string spatScheme, std::string fName); */
    
    // Destructor
    ~Probe();

    // Headers
    /* void checkProbes(Mesh Msh, Solver* Sol, double t=0); */
};

inline std::string createFolder(std::string fName, std::string& dirName) {
    // Timestamp
    time_t timeStamp = std::time(nullptr);
    struct tm datetime = *localtime(&timeStamp);
    char oName[35]; strftime(oName, sizeof(oName), "%Y%m%d%H%M%S_", &datetime);
    
    // Folder Name
    std::filesystem::path pBase = std::filesystem::current_path();
    size_t iPos = fName.find(".json"); dirName = oName + fName.substr(0, iPos);
    pBase /= dirName;

    // Create Folder
    std::filesystem::create_directories(pBase);

    return pBase.string();
}

inline std::ofstream createFile(std::filesystem::path fName){
    // Open File 
    std::ofstream file(fName);
    if (!file.is_open()){ throw std::runtime_error("File could not be created: " + fName.string()); }

    // File Header
    file << "Time";

    return file;
}

#endif
