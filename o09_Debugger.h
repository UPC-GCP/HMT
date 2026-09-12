#ifndef DEBUGGER_H_
#define DEBUGGER_H_

// Imports
#include <iostream>
#include <cstddef>

// Self-Imports
/* #include "o02_Mesh.h" */

// Types
struct debugOptions {
    bool bGeneral=false, bMat=false, bSurf=false, bVol=false, bObs=false; // General
    bool bBoundaries=false, bPhi=false, bPhiBC=false; // Boundary
};

// Headers
/* void print1D(MeshSolver<1> Msh, debugOptions dOps); */
/* void print2D(MeshSolver<2> Msh, debugOptions dOps); */
template <template <size_t> class MeshClass> void print3D(const MeshClass<3>& Msh, debugOptions dOps);
template <template <size_t> class MeshClass> void print2D(const MeshClass<2>& Msh, debugOptions dOps);
template <template <size_t> class MeshClass> void print1D(const MeshClass<1>& Msh, debugOptions dOps);

// Functions
template <template <size_t> class MeshClass, size_t Dim> void printDebug(const MeshClass<Dim> Msh, debugOptions dOps){

    // General Mesh
    if (dOps.bGeneral) {
        for (size_t i = 0; i < Dim; i++) {
            std::cout << "\nAxis: " << i << "\n";
            std::cout << "Faces: "; for (double val : Msh.Faces[i]) {std::cout << val << " ";} std::cout << "\n";
            std::cout << "Nodes: "; for (double val : Msh.Nodes[i]) {std::cout << val << " ";} std::cout << "\n";
            std::cout << "DeltaX: "; for (double val : Msh.deltaX[i]) {std::cout << val << " ";} std::cout << "\n";
            std::cout << "dX: "; for (double val : Msh.dX[i]) {std::cout << val << " ";} std::cout << "\n";
        } std::cout << "\n";
    }
    
    // General Boundaries
    if (dOps.bBoundaries) {
        for (size_t i = 0; i < Msh.BC.size(); i++) {
            std::cout << "\nBC: " << i << "\n";
            std::cout << "Type - Side: " << Msh.BC[i].type << " " << Msh.BC[i].side << "\n";
            std::cout << "Expression: " << Msh.BC[i].bUpdate << " - " << Msh.BC[i].iExpr << " " << Msh.BC[i].iEq << " " << Msh.BC[i].expression << "\n";
            std::cout << "ExpressionA: " << Msh.BC[i].bA << " - " << Msh.BC[i].iExprA << " " << Msh.BC[i].iEqA << " " << Msh.BC[i].expressionA << "\n";
            std::cout << "Pos0: "; for (size_t val : Msh.BC[i].i0) { std::cout << val << " "; } std::cout << "\n";
            std::cout << "Pos1: "; for (size_t val : Msh.BC[i].i1) { std::cout << val << " "; } std::cout << "\n";
        } std::cout << "\n";
    }

    // Dimensional
    if constexpr (Dim == 1) {print1D(Msh, dOps);}
    else if constexpr (Dim == 2) {print2D(Msh, dOps);}
    else if constexpr (Dim == 3) {print3D(Msh, dOps);}

}

#endif
