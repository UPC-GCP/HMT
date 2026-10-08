#ifndef MEDIC_H_
#define MEDIC_H_

// Imports
#include <cstddef>
#include <filesystem>

// Self-Imports
#include "o05_Probe.h"

// Class
template <std::size_t Dim> class Medic {
private:

public:
    // Variables
    std::filesystem::path fPath{};

    // Medic
    Medic(const Probe<Dim>& Prb);
};

// Constructor
template <std::size_t Dim> Medic<Dim>::Medic(const Probe<Dim>& Prb) {
    // Folder
    fPath = Prb.fPath;
}

#endif
