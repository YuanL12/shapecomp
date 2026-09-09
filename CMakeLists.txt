cmake_minimum_required(VERSION 3.22)
project(shapecomp LANGUAGES CXX)

option(SHAPECOMP_BUILD_TESTS "Build ShapeComp native smoke/test executables" OFF)
option(SHAPECOMP_BUILD_EXAMPLES "Build ShapeComp example executables" OFF)

# Disable Taskflow demos when embedding the library.
set(TF_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(TF_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fdiagnostics-color=always")

# Set default build type to Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Prefer Ninja generator if available
if(NOT CMAKE_GENERATOR)
    find_program(NINJA_EXECUTABLE ninja)
    if(NINJA_EXECUTABLE)
        set(CMAKE_GENERATOR "Ninja" CACHE STRING "CMake generator" FORCE)
    endif()
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    # Debug flags: no optimization, debug symbols, comprehensive warnings
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wuninitialized")
endif()

# Set all libraries to position independent code
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# C++ flags
set(CMAKE_CXX_STANDARD 17) 
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# CMake module path
list(PREPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Add submodules
add_subdirectory(deps/eigen)
add_subdirectory(deps/geometry-central)
add_subdirectory(deps/taskflow)
add_subdirectory(deps/ReebHanTun)

# ----------------- Build the library -----------------
# Add Taskflow from deps/
add_library(taskflow INTERFACE)
target_include_directories(taskflow INTERFACE ${CMAKE_SOURCE_DIR}/deps/taskflow)

# create shapecomp static library
list(APPEND GP_INCLUDE_DIRS "include")
# Use explicit file lists instead of GLOB to prevent CMake from re-running unnecessarily
set(GP_SOURCES
    # mesh sources
    src/mesh/Corner.cpp
    src/mesh/Edge.cpp
    src/mesh/Face.cpp
    src/mesh/HalfEdge.cpp
    src/mesh/Mesh.cpp
    src/mesh/MeshIO.cpp
    src/mesh/PolygonSoup.cpp
    src/mesh/Vertex.cpp
    # generators sources
    src/generators/Dijkstra.cpp
    src/generators/Generators.cpp
    src/generators/ReebLoop.cpp
    # matrix sources
    src/matrix/Laplacian.cpp
    # embedding sources
    src/embedding/torus/Tutte.cpp
    # torus_comparison sources
    src/torus_comparison/GeodesicPath.cpp
    src/torus_comparison/MeshSurgery.cpp
    src/torus_comparison/Planar.cpp
    # extern/geometry_central sources
    src/extern/geometry_central/geodesic.cpp
)
add_library(gp STATIC ${GP_SOURCES})
set_target_properties(gp PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(gp PUBLIC ${GP_INCLUDE_DIRS})
target_link_libraries(gp PUBLIC Eigen3::Eigen taskflow geometry-central ReebHanTun::ReebHanTunLib)

# ----------------- Native tests / examples -----------------
if(SHAPECOMP_BUILD_TESTS OR SHAPECOMP_BUILD_EXAMPLES)
    enable_testing()

    add_executable(main src/test/main.cpp)
    target_link_libraries(main PRIVATE gp)
    target_link_libraries(main PRIVATE eigen)
    set_target_properties(main PROPERTIES OUTPUT_NAME main_out)

    add_executable(load_mesh src/test/load_mesh.cpp)
    target_link_libraries(load_mesh PRIVATE gp)
    set_target_properties(load_mesh PROPERTIES OUTPUT_NAME load_mesh.out)

    add_executable(test_flip_geodesic src/test/test_gc_flip_geodesic.cpp)
    set_target_properties(test_flip_geodesic PROPERTIES OUTPUT_NAME test_flip_geodesic.out)
    target_link_libraries(test_flip_geodesic PRIVATE Eigen3::Eigen gp)

    add_executable(test_planar_locator src/test/test_planar_locator.cpp)
    set_target_properties(test_planar_locator PROPERTIES OUTPUT_NAME test_planar_locator.out)
    target_link_libraries(test_planar_locator PRIVATE Eigen3::Eigen gp)

    add_executable(test_CircularLinkedList src/test/test_CircularLinkedList.cpp)
    set_target_properties(test_CircularLinkedList PROPERTIES OUTPUT_NAME test_CircularLinkedList.out)
    target_link_libraries(test_CircularLinkedList PRIVATE Eigen3::Eigen gp)

    add_executable(test_ReebHanTen src/test/test_ReebHanTen.cpp)
    set_target_properties(test_ReebHanTen PROPERTIES OUTPUT_NAME test_ReebHanTen.out)
    target_link_libraries(test_ReebHanTen PRIVATE Eigen3::Eigen gp)

    # Smoke: load tracked torus fixture and exit 0.
    add_test(
        NAME load_mesh_torus_smoke
        COMMAND load_mesh "${CMAKE_SOURCE_DIR}/test/input/torus/torus.obj"
    )
endif()

# ----------------- Python binding -----------------
# scikit-build-core provides pybind11 for the Python package build. Pure C++
# configures (tests/examples) can skip the binding when pybind11 is absent.
option(SHAPECOMP_BUILD_PYTHON "Build the pybind11 Python extension" ON)
if(SHAPECOMP_BUILD_PYTHON)
    find_package(pybind11 CONFIG QUIET)
    if(pybind11_FOUND)
        add_subdirectory(binding)
    else()
        message(STATUS "pybind11 not found; skipping Python binding (set SHAPECOMP_BUILD_PYTHON=OFF to silence)")
    endif()
endif()
