#ifndef PROBEFVM_H_
#define PROBEFVM_H_

// Imports
/* #include <cstddef> */

// Self-Imports
#include "o05_Probe.h"
#include "o02_MeshFVM.h"
#include <cstddef>

// Enumerator
enum class prbFVM{ Point, Field };

// Types 
template <std::size_t Dim> struct probeFVM : probeRange<Dim> {
    size_t type{};
};

// Class
template <std::size_t Dim> class ProbeFVM : Probe<Dim> {
private:

public:
    // Variables
    /* std::vector<probeFVM<Dim>> pModal{}; */ // Whatever holds the probes and should be deleted in destructor
    
    // Constructor
    ProbeFVM(const MeshFVM<Dim>& Msh, const Json::Value& probes, std::string fName);

    // Destructor
    ~ProbeFVM();

    // Headers
    void checkProbes(const MeshFVM<Dim>& Msh, double t);
};

#endif
