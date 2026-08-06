#ifndef PROBE_H_
#define PROBE_H_

#include <cstddef>
#include <fstream>
#include <vector>
#include <json/json.h>

#include "o02_Mesh.h"
#include "o04_Solver.h"

struct Prb{
    int type{};
    std::vector<double> x0{}, x1{};
    double t{};
    std::ofstream file{};
};

struct pBug{
    std::ofstream file{};
    std::vector<size_t> xPos{}, yPos{};
    std::vector<double> t;
};

struct pFld{
    std::ofstream file{};
    std::vector<size_t> xPos{}, yPos{};
    std::vector<double> t;
};

struct pMap{
    std::ofstream file{};
    std::vector<size_t> xPos{}, yPos{};
    std::vector<double> t;
};

struct pBnd{
    std::ofstream file{};
    int side{}, nPos{};
    std::vector<size_t> iBC{}, iPos{};
    std::vector<double> t;
    std::vector<double> xPos{}, yPos{};
};

struct pPoint{
    bool bFile = true;
    std::ofstream file{};
    std::vector<int> xPos{}, yPos{};
    std::vector<double> t0{}, t1{}; 
};

class Probe
{
private:

public:

    // Variables
    std::string pathBase{};
    std::string dirName{};
    pPoint probePoint{};

    // Vectors
    std::vector<pMap> probeMap{};
    std::vector<pFld> probeFld{};
    std::vector<pBug> probeBug{};
    std::vector<pBnd> probeBnd{};

    // Constructor
    Probe(Mesh Msh, Json::Value probes, std::string tempScheme, std::string spatScheme, std::string fName);
    
    // Destructor
    ~Probe();

    // Functions
    void checkProbes(Mesh Msh, Solver* Sol, double t=0);
    
};

#endif
