#ifndef MEDICNS_H_
#define MEDICNS_H_

#include <string>

#include "o01_Material.h"
#include "o02_MeshSC.h"
#include "o03_DiscretizerSC.h"
#include "o05_ProbeSC.h"

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
    void getDiagnostic(Material Mat, Mesh Msh, Discretizer Dsc, double t);
    void getGlobalBalance(Material Mat, Mesh Msh, Discretizer Dsc);
    void getSystemResidual(Material Mat, Mesh Msh, Discretizer Dsc, double t);
};

#endif
