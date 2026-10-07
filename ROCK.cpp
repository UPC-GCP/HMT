// Imports
#include <chrono>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <json/json.h>

////////// JSON PARSER ///////////
Json::Value getParsedData(const std::string& fileName) {
    ///// Open File /////
    std::ifstream file(fileName, std::ifstream::binary);

    ///// Filter /////
    if (!file.is_open()) { throw std::runtime_error("Configuration file could not be read."); }

    ///// Parsing /////
    Json::Value data; Json::CharReaderBuilder readerBuilder; std::string errs;
    Json::parseFromStream(readerBuilder, file, &data, &errs);

    ///// Close File /////
    file.close();

    return data;
}

////////// SOLVER TYPE //////////

///// Enumerator /////
enum class enumSolver { FVMScalar, FVMNS, Spectral };
enumSolver getSolverType(const std::string& configSolver) {
    if (configSolver == "FVMScalar") { return enumSolver::FVMScalar; }
    else if (configSolver == "FVMNS") { return enumSolver::FVMNS; }
    else if (configSolver == "Spectral") { return enumSolver::Spectral; }
    else throw std::logic_error("Unknown solver type: " + configSolver);
}

///// Forward Declaration /////
void runFVMScalar(const Json::Value& data, std::string configName);
void runFVMNS(const Json::Value& data, std::string configName);
void runSpectral(const Json::Value& data, std::string configName);

////////// MAIN //////////
int main(int argc, char* argv[]){ // ROCK: Research Oriented Computational Kernel
    ///// Setup /////
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Initializing model ... \n" << std::fixed << std::setprecision(2); // std::fixed, std::scientific, std::hexfloat, std::defaultfloat

    try {
        ///// Data /////
        std::cout << "Reading data ... \n";
        Json::Value data = getParsedData(argv[1]); std::cout << "Data parsed successfully. \n";

        ///// SOLVER SELECTION /////
        enumSolver configSolver = getSolverType(data["configSolver"].asString());
        switch (configSolver) {
            case enumSolver::FVMScalar: runFVMScalar(data, argv[1]); break;
            case enumSolver::FVMNS: runFVMNS(data, argv[1]); break;
            case enumSolver::Spectral: runSpectral(data, argv[1]); break;
        }
    } catch (const std::exception& e) {std::cerr << "Program shutdown...\n"; return EXIT_FAILURE;}

    ///// Control /////
    auto t2 = std::chrono::high_resolution_clock::now(); std::chrono::duration<double, std::milli> msDoub = t2 - t1; double tTime = msDoub.count()/1000/60; std::cout << "Time elapsed: " << int(tTime) << " minutes and " << (tTime - int(tTime))*60 << " seconds.\n";
}
