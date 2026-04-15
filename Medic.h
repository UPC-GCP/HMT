#ifndef MEDIC_H_
#define MEDIC_H_

#include <vector>
#include <string>
#include <iostream>

#include "Material.h"
#include "Mesh.h"
#include "Discretizer.h"
#include "Probe.h"

class Medic
{
private:

public:
    // Variables
    std::string pathBase{};
    std::ofstream file{}, fileR{};
    
    // Constructor
    Medic(Mesh Msh, Probe& Prb);

    // Functions
    void getDiagnostic(Material Mat, Mesh Msh, Discretizer Dsc, std::vector<std::vector<double>> oldTemp, double t);
    void getGlobalBalance(Material Mat, Mesh Msh, Discretizer Dsc);
    void getSystemResidual(Material Mat, Mesh Msh, Discretizer Dsc);
};

#endif