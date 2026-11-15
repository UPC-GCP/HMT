#ifndef PROBESPECTRAL_H_
#define PROBESPECTRAL_H_

// Imports

// Self-Imports
#include "o05_Probe.h"
#include "o02_MeshSpectral.h"

// Enumerator
enum class prbSpectral { E, u };

// Types 
struct probeRangeSpectral : probeRange<1> {
    prbSpectral type{};
};

// Class
class ProbeSpectral : public Probe<1> {
private:

public:
    // Variables
    std::vector<probeRangeSpectral> pModal{};
    
    // Constructor
    ProbeSpectral(const MeshBurgers& Burg, const Json::Value& probes, std::string fName);

    // Destructor
    ~ProbeSpectral();

    // Headers
    void checkProbes(const MeshBurgers& Burg, double t);
};

#endif
