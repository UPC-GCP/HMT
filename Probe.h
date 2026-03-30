#ifndef PROBE_H_
#define PROBE_H_

#include <iostream>
#include <vector>
#include <json/json.h>

#include "Mesh.h"

struct Prb{
    int type{};
    std::vector<double> x0{}, x1{};
    double t{};
    std::ofstream file{};
};

struct pMap{
    std::ofstream file{};
    std::vector<size_t> xPos{}, yPos{};
    std::vector<double> t;
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
    pPoint probePoint{};

    // Vectors
    std::vector<pMap> probeMap{};

    // Constructor
    Probe(Mesh Msh, Json::Value probes, std::string scheme, std::string fName);
    
    // Destructor
    ~Probe();

    // Functions
    void checkProbes(Mesh Msh, double t=0);
    
};

#endif