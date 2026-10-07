#ifndef PROBESPECTRAL_H_
#define PROBESPECTRAL_H_

// Imports

// Self-Imports
#include "o05_Probe.h"
#include "o02_MeshSpectral.h"

// Class
class ProbeSpectral : Probe<1> {
private:

public:
    // Variables
    std::vector<probeRange<1>> pModal{};

    // Constructor
    ProbeSpectral(const MeshBurgers& Burg, const Json::Value& probes, std::string fName);

    // Destructor
    ~ProbeSpectral();
};


#endif
