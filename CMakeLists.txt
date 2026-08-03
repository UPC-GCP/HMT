cmake_minimum_required(VERSION 3.10)
project(HMT)

# 1. C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# 2. Packages 
find_package(PkgConfig REQUIRED)
pkg_check_modules(JSONCPP REQUIRED IMPORTED_TARGET jsoncpp)

# 3. Executable
add_executable(MainSolver
    o00_MainSolver.cpp
    o01_Material.cpp
    o02_Mesh.cpp
    o03_Discretizer.cpp
    o04_Solver.cpp
    o04_CG.cpp
    o04_BCG.cpp
    o04_GS.cpp
    o05_Probe.cpp
    o09_ExpressionParser.cpp
    o09_Medic.cpp
)

add_executable(LDCSolver
    o00_NSSolver.cpp
    o01_Material.cpp
    o02_MeshNS.cpp
    o03_DiscretizerNS.cpp
    o04_SolverNS.cpp
    o04_CGNS.cpp
    o04_BCG.cpp
    o05_ProbeNS.cpp
    o09_ExpressionParser.cpp
    o09_MedicNS.cpp
)

add_executable(DHCSolver
    o00_DHCSolver.cpp
    o01_Material.cpp
    o02_MeshNS.cpp
    o03_DiscretizerNS.cpp
    o04_SolverNS.cpp
    o04_CGNS.cpp
    o04_BCG.cpp
    o05_ProbeNS.cpp
    o09_ExpressionParser.cpp
    o09_MedicNS.cpp
)

add_executable(SCSolver
    o00_SCSolver.cpp
    o01_Material.cpp
    o02_MeshSC.cpp
    o03_DiscretizerSC.cpp
    o04_SolverSC.cpp
    o04_CGSC.cpp
    # o04_BCG.cpp
    o05_ProbeSC.cpp
    o09_ExpressionParser.cpp
    o09_MedicSC.cpp
)

add_executable(3DSolver
    o00_3DSolver.cpp
    o01_Material.cpp
    o02_Mesh3D.cpp
#     o03_Discretizer3D.cpp
#     o04_Solver3D.cpp
#     o04_CG3D.cpp
#     o04_BCG3D.cpp
#     o05_Probe3D.cpp
#     o09_Medic3D.cpp
    o09_ExpressionParser.cpp
)

# 4. Linking
foreach (target MainSolver LDCSolver DHCSolver SCSolver 3DSolver)
    target_link_libraries(${target} PRIVATE PkgConfig::JSONCPP)
    target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
endforeach()



