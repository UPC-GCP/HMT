#ifndef MEDIC_H_
#define MEDIC_H_

#include <vector>
#include <string>

#include "o01_Material.h"
#include "o02_Mesh.h"
#include "o03_Discretizer.h"
#include "o05_Probe.h"

class Medic
{
private:

public:
    // Variables
    std::string pathBase{};
    std::ofstream file{}, fileR{};
    
    // Constructor
    Medic(Mesh Msh, Probe& Prb, bool bExit);

    // Destructor
    ~Medic();

    // Functions
    void getDiagnostic(Material Mat, Mesh Msh, Discretizer Dsc, std::vector<std::vector<double>> oldPhi, double t);
    void getGlobalBalance(Material Mat, Mesh Msh, Discretizer Dsc);
    void getSystemResidual(Material Mat, Mesh Msh, Discretizer Dsc, double t);
};

#endif
