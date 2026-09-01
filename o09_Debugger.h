#ifndef DEBUGGER_H_
#define DEBUGGER_H_

// Self-Imports
#include "o02_Mesh.h"

struct debugOptions {
    bool bGeneral=false, bMat=false, bSurf=false, bVol=false, bObs=false; // General Mesh
    bool bBoundaries=false, bPhi=false; // General Boundary
};

void print1D(MeshSolver<1> Msh, debugOptions dOps);
void print2D(MeshSolver<2> Msh, debugOptions dOps);
void print3D(MeshSolver<3> Msh, debugOptions dOps);

#endif
