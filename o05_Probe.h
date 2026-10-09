#ifndef PROBE_H_
#define PROBE_H_

// Imports
/* #include <iostream> */
#include <ctime>
#include <array>
#include <vector>
#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <filesystem>
#include <json/json.h>

template <size_t Dim> struct probeBase {
    std::ofstream file{}; // File
    size_t nWrite{1}, nCount{}; // Counters
};

template <size_t Dim> struct probePoint : probeBase<Dim> { // Each coordinate i0[k] is stored within timestamps t[k]
    std::vector<std::array<double, 2>> t{}; // Time interval
    std::vector<std::array<size_t, Dim>> i0{}; // Coordinates
};

template <size_t Dim> struct probeRange : probeBase<Dim> { // Each range i0, i1 is stored within timestamps t[k]
    /* size_t type{}; // Need to make it detect P, T, V for FVM, maybe this goes within */
    std::array<double, 2> t{}; // Time interval
    std::array<size_t, Dim> i0{}, i1{}; // Coordinates
};

template <size_t Dim> class Probe {
private:

public:
    // Variables
    size_t iProbes{};
    std::string pathBase{};
    std::filesystem::path fPath;

    // Constructor
    Probe(const Json::Value& probes, std::string fName);
};

// Functions
inline std::string createFolder(std::string fName) {
    // Timestamp
    time_t timeStamp = std::time(nullptr); struct tm datetime = *localtime(&timeStamp);
    char oName[35]; strftime(oName, sizeof(oName), "%Y%m%d%H%M%S_", &datetime);
    
    // User
    const char* user = std::getenv("USER");
    (user == nullptr) ? user = std::getenv("LOGNAME") : nullptr;
    std::string uName = std::string(user);

    // Directory
    std::filesystem::path pBase = std::filesystem::current_path(); 
    if (uName == "upc_gcp") { pBase /= "ioRes"; } else if (uName == "gonzalo") { pBase /= "../../../home_nobck/ioRes"; pBase = std::filesystem::weakly_canonical(pBase); }

    // Create Folder
    size_t iPos = fName.find(".json"); std::string dirName = oName + fName.substr(0, iPos); pBase /= dirName;
    std::filesystem::create_directories(pBase);

    return pBase.string();
}

inline std::ofstream createFile(const std::filesystem::path& fName){
    // Create File
    std::ofstream file(fName); if (!file.is_open()){ throw std::runtime_error("File could not be created: " + fName.string()); } file << "Time";
    
    return file;
}

// Constuctor
template <size_t Dim> Probe<Dim>::Probe(const Json::Value& probes, std::string fName) {
    // Control
    fPath = fName; if (probes.size() == 0) { return; }

    // Create Folder
    pathBase = createFolder(fPath.filename().string()); fPath = pathBase;
}

#endif
