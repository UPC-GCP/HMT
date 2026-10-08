cmake_minimum_required(VERSION 3.10)
project(HMT)

# 1. C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# 2. Packages 
find_package(OpenMP REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(JSONCPP REQUIRED IMPORTED_TARGET jsoncpp)

# 3. Executable
add_executable(ROCK
    ROCK.cpp
    o00_FVMNS.cpp
    o00_FVMScalar.cpp
    o00_Spectral.cpp
    o01_Material.cpp
    # o01_MaterialFVM.cpp
    # o01_MaterialSpectral.cpp - Can split this later but honestly doesn't cost too much
    o02_MeshFVM.cpp
    o02_MeshSpectral.cpp
    # o03_DiscretizerFVM.cpp
    o03_DiscretizerSpectral.cpp
    # o04_SolverFVM.cpp
    o04_SolverSpectral.cpp
    # o05_ProbeFVM.cpp
    o05_ProbeSpectral.cpp
    # o06_MedicFVM.cpp
    o06_MedicSpectral.cpp
    o09_Parser.cpp
    # o09_Debugger.cpp
    # o09_DebuggerFVM.cpp
    # o09_DebuggerSpectral.cpp
)

# 4. Linking
foreach (target ROCK)
    target_compile_options(${target} PRIVATE -O3 -march=native) # gcc
    target_link_libraries(${target} PRIVATE PkgConfig::JSONCPP) # jsoncpp
    target_link_libraries(${target} PRIVATE OpenMP::OpenMP_CXX) # OpenMP
    target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
endforeach()
