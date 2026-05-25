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
    # o04_Solver.cpp
    # o04_CG.cpp
    # GS.cpp
    # o05_Probe.cpp
    o09_ExpressionParser.cpp
    # o09_Medic.cpp
)

add_executable(TestSolver
    MainSolver.cpp
    Material.cpp
    Mesh.cpp
    Discretizer.cpp
    Solver.cpp
    CG.cpp
    Probe.cpp
    ExpressionParser.cpp
    Medic.cpp
)

# # 4. Directories
# target_include_directories(MainSolver PRIVATE
# 	${JSONCPP_INCLUDE_DIRS}
# 	${CMAKE_CURRENT_SOURCE_DIR}
# )

# # 4. Link Libraries
# target_link_libraries(MainSolver PRIVATE ${JSONCPP_LIBRARIES})

# 4. Linking
foreach (target MainSolver TestSolver)
    target_link_libraries(${target} PRIVATE PkgConfig::JSONCPP)
    target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
endforeach()



