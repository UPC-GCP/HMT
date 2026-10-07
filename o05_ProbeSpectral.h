#ifndef PROBESPECTRAL_H_
#define PROBESPECTRAL_H_

// Imports

// Self-Imports
#include "o05_Probe.h"

// Class
class ProbeSpectral : Probe<1> {
private:

public:
    // Variables
    std::vector<probeRange<1>> pModal{};

    // Constructor
    ProbeSpectral(const Json::Value& probes, std::string fName);

    // Destructor
    ~ProbeSpectral();
};


#endif
