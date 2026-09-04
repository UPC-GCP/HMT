// Imports
#include <cstddef>
#include <iostream>

// Self-Imports
#include "o09_Debugger.h"
#include "o02_Mesh.h"

// PENDING CHANGES:
// Detect type of variable and block out corresponding sections (MeshBase/MeshSolver)

void print3D(MeshSolver<3> Msh, debugOptions dOps) {
    // Geometry
    if (dOps.bSurf) { // Surfaces
        for (size_t nD = 0; nD < 3; nD++) {
            std::cout << "S " << nD << ":\n";
            for (size_t k = 0; k < Msh.N[2]; k++) {
                for (size_t i = 0; i < Msh.N[0]; i++) {
                    for (size_t j = 0; j < Msh.N[1]; j++) {
                        std::cout << Msh.S[nD][calcIndex(i, j, Msh.N[1], k, Msh.N[0])] << " ";
                    } std::cout << "\n";
                } std::cout << "\n";
            } std::cout << "\n";
        } std::cout << "\n";
    }

    if (dOps.bVol) { // Volume
        std::cout << "V:\n";
        for (size_t k = 0; k < Msh.N[2]; k++) {
            for (size_t i = 0; i < Msh.N[0]; i++) {
                for (size_t j = 0; j < Msh.N[1]; j++) {
                    std::cout << Msh.Vp[calcIndex(i, j, Msh.N[1], k, Msh.N[0])] << " ";
                } std::cout << "\n";
            } std::cout << "\n";
        } std::cout << "\n";
    }
    
    // MeshSolver
    /* if constexpr (std::is_same_v<T, MeshSolver<3>>) { */
        // Material
        if (dOps.bMat) {
            std::cout << "nMat:\n"; 
            for (size_t k = 0; k < Msh.N[2]; k++) {
                for (size_t i = 0; i < Msh.N[0]; i++) {
                    for (size_t j = 0; j < Msh.N[1]; j++) {
                        std::cout << Msh.nMat[calcIndex(i, j, Msh.N[1], k, Msh.N[2])] << " ";
                    } std::cout << "\n";
                } std::cout << "\n";
            } std::cout << "\n";
        }

        // Obstacles
        if (dOps.bObs) {
            std::cout << "Obs:\n";
            for (size_t k = 0; k < Msh.N[2]; k++) {
                for (size_t i = 0; i < Msh.N[0]; i++) {
                    for (size_t j = 0; j < Msh.N[1]; j++) {
                        std::cout << Msh.bObs[calcIndex(i, j, Msh.N[1], k, Msh.N[0])] << " ";
                    } std::cout << "\n";
                } std::cout << "\n";
            } std::cout << "\n";
        }
    /* } */

    // Value
    if (dOps.bPhi) {
        std::cout << "Phi:\n";
        for (size_t k = 0; k < Msh.N[2]; k++) {
            for (size_t i = 0; i < Msh.N[0]; i++) {
                for (size_t j = 0; j < Msh.N[1]; j++) {
                    std::cout << Msh.Phi[calcIndex(i, j, Msh.N[1], k, Msh.N[0])] << " ";
                } std::cout << "\n";
            } std::cout << "\n";
        } std::cout << "\n";
    }
        /* runLoopMesh<3>(Msh.N, [&](size_t i, size_t j, size_t k) { std::cout << Msh.Phi[calcIndex(i, j, k)] << "\n";}); */
        // Maybe won't work for this but should work in calculations
}

void print2D(MeshSolver<2> Msh, debugOptions dOps) {
    // Geometry
    if (dOps.bSurf) { // Surfaces
        for (size_t nD = 0; nD < 2; nD++) {
            std::cout << "S " << nD << ":\n";
            for (size_t i = 0; i < Msh.N[0]; i++) {
                for (size_t j = 0; j < Msh.N[1]; j++) {
                    std::cout << Msh.S[nD][calcIndex(i, j, Msh.N[1])] << " ";
                } std::cout << "\n";
            } std::cout << "\n";
        } std::cout << "\n";
    }

    if (dOps.bVol) { // Volume
        std::cout << "V:\n";
        for (size_t i = 0; i < Msh.N[0]; i++) {
            for (size_t j = 0; j < Msh.N[1]; j++) {
                std::cout << Msh.Vp[calcIndex(i, j, Msh.N[1])] << " ";
            } std::cout << "\n";
        }
    }

    // MeshSolver
    /* if constexpr (std::is_same_v<T, MeshSolver<2>>) { */
        // Material
        if (dOps.bMat) {
            std::cout << "nMat:\n"; 
            for (size_t i = 0; i < Msh.N[0]; i++) {
                for (size_t j = 0; j < Msh.N[1]; j++) {
                    std::cout << Msh.nMat[calcIndex(i, j, Msh.N[1])] << " ";
                } std::cout << "\n";
            } std::cout << "\n";
        }

        // Obstacles
        if (dOps.bObs) {
            std::cout << "Obs:\n";
            for (size_t i = 0; i < Msh.N[0]; i++) {
                for (size_t j = 0; j < Msh.N[1]; j++) {
                    std::cout << Msh.bObs[calcIndex(i, j, Msh.N[1])] << " ";
                } std::cout << "\n";
            } std::cout << "\n";
        }
    /* } */

    // Value
    if (dOps.bPhi) {
        std::cout << "Phi:\n";
        runLoopMesh<2>(Msh.N,  [&](size_t i, size_t j, size_t k) { std::cout << Msh.Phi[calcIndex(i, j)] << "\n";});
    }
}

void print1D(MeshSolver<1> Msh, debugOptions dOps) {
    // Geometry
    if (dOps.bSurf) { // Surfaces
        for (size_t nD = 0; nD < 1; nD++) {
            std::cout << "S " << nD << ":\n";
            for (size_t i = 0; i < Msh.N[0]; i++) {
                std::cout << Msh.S[nD][calcIndex(i)] << " ";
            } std::cout << "\n";
        } std::cout << "\n";
    }

    if (dOps.bVol) { // Volume
        std::cout << "V:\n";
        for (size_t i = 0; i < Msh.N[0]; i++) {
            std::cout << Msh.Vp[calcIndex(i)] << " ";
        } std::cout << "\n";
    }

    // MeshSolver
    /* if constexpr (std::is_same_v<T, MeshSolver<1>>) { */
        // Material
        if (dOps.bMat) {
            std::cout << "nMat:\n";
            for (double val : Msh.nMat) {
                std::cout << val << " ";
            } std::cout << "\n"; std::cout << "\n";
        }

        // Obstacles
        if (dOps.bObs) {
            std::cout << "Obs:\n";
            for (size_t i = 0; i < Msh.N[0]; i++) {
                std::cout << Msh.bObs[calcIndex(i)] << " ";
            } std::cout << "\n";
        }
    /* } */

    // Value
    if (dOps.bPhi) {
        std::cout << "Phi\n";
        runLoopMesh<1>(Msh.N, [&](size_t i, size_t j, size_t k) { std::cout << Msh.Phi[calcIndex(i)] << "\n";});
    }
}
