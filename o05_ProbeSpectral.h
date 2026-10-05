#ifndef PROBESPECTRAL_H_
#define PROBESPECTRAL_H_

// Imports

// Self-Imports
#include "o05_Probe.h"

// Types
struct probeModes : probeBase<1> {
    // Not sure if I need this or if I just defined the object as pModes = probeRange<1>
};

// Class
class ProbeSpectral : Probe<1> {
private:

public:
    // Variables

    // Constructor
    ProbeSpectral(const Json::Value& probes, std::string fName);
};


#endif
