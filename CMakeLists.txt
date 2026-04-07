cmake_minimum_required(VERSION 3.10)
project(HMT)

# 1. Use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 2. Find JsonCpp (vcpkg makes this work)
add_definitions(-DJSON_DLL)
find_package(jsoncpp CONFIG REQUIRED)

# 3. Define your executable
add_executable(MainSolver
    MainSolver.cpp
    Material.cpp
    Mesh.cpp
    Discretizer.cpp
    Solver.cpp
    CG.cpp
    GS.cpp
    ExpressionParser.cpp
    Probe.cpp
)

# Add BOTH the root include and the specific json folder just in case
target_include_directories(MainSolver PRIVATE 
    "C:/Users/gonce/vcpkg/installed/x64-mingw-dynamic/include"
    "C:/Users/gonce/vcpkg/installed/x64-mingw-dynamic/include/json"
    "C:\\Users\\gonce\\vcpkg\\installed\\x64-mingw-dynamic\\include\\json"
)

target_compile_definitions(MainSolver PRIVATE JSON_DLL)

# This flag allows the object files to handle the massive ExprTk templates
target_compile_options(MainSolver PRIVATE "-Wa,-mbig-obj")

# 4. Link the library (This handles the .h and .dll links automatically)
# target_link_libraries(MainSolver PRIVATE JsonCpp::JsonCpp)

# Manual path fix (ensure this path is correct for your PC)
target_link_libraries(MainSolver PRIVATE "C:/Users/gonce/vcpkg/installed/x64-mingw-dynamic/lib/libjsoncpp.dll.a")

target_include_directories(MainSolver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
