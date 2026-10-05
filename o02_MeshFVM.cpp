// Self-Imports
#include "o02_MeshFVM.h"
#include <cstddef>


template <size_t Dim> void MeshFVM<Dim>::generateMeshSolver(MeshSolver<Dim>& Msh, Json::Value qNode, Json::Value sections, Json::Value refinement, Json::Value obstacles) {

}


template <size_t Dim> void MeshFVM<Dim>::addBoundariesSolver(MeshSolver<Dim>& Msh, Material Mat, Parser& Prs, Json::Value boundaries, double dInit, std::string sInit) {

}

template <size_t Dim> void MeshFVM<Dim>::deriveMeshBase(MeshSolver<Dim> p, std::array<MeshBase<Dim>, Dim>& V) {

}


    // Was it worth it to pass the std::array<MeshBase<Dim>, Dim>& V ? Or should I just pass MeshBase<Dim>
    // Test with addBoundariesBase() -- The code doesn't know what index to access if I send them independently, would need to pass that as well -- in main() -> for (size_t i = 0; i < Dim; i++) { Msh.addBoundariesBase(i, V, Mat, Prs, boundaries, bInit, sInit) }


template <size_t Dim> void MeshFVM<Dim>::addBoundariesBase(size_t i, std::array<MeshBase<Dim>, Dim>& V, Material Mat, Parser& Prs, Json::Value boundaries, std::vector<double> dInit, std::vector<std::string> sInit) {

}


// Compiler Instances
template class MeshFVM<1>;
template class MeshFVM<2>;
template class MeshFVM<3>;
