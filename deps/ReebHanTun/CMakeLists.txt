#set(CMAKE_CXX_FLAGS "-lboost_program_options")
#include_directories(${CMAKE_SOURCE_DIR}/src "${CMAKE_SOURCE_DIR}/src/include" "${CMAKE_SOURCE_DIR}/src/CycleOptimization" "${CMAKE_SOURCE_DIR}/src/CycleOptimization/Annotation")

cmake_minimum_required(VERSION 3.22)

# Suppress CMP0167 warning about FindBoost being removed
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 OLD)
endif()

# Force use of FindBoost instead of BoostConfig (for compatibility with Boost 1.89)
# Must be set BEFORE find_package and unset any Boost_DIR that might force BoostConfig
unset(Boost_DIR CACHE)
set(Boost_NO_BOOST_CMAKE ON)
# Find Boost - in Boost 1.89, system is header-only, so we only require program_options
find_package(Boost 1.40 COMPONENTS program_options REQUIRED)
# boost_system is header-only in Boost 1.89, use Boost::boost for header-only components
set(BOOST_SYSTEM_LIB Boost::boost)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

# CGAL
find_package(CGAL REQUIRED)
if(CGAL_USE_FILE)
    include(${CGAL_USE_FILE})
else()
    # Newer CGAL versions don't use CGAL_USE_FILE
    find_package(CGAL COMPONENTS REQUIRED)
endif()

# ANN - build from submodule automatically, then fallback to system
# CMAKE_SOURCE_DIR here is the project root (set by parent CMakeLists.txt)
get_filename_component(PROJECT_ROOT "${CMAKE_CURRENT_SOURCE_DIR}" DIRECTORY)
get_filename_component(PROJECT_ROOT "${PROJECT_ROOT}" DIRECTORY)
set(ANN_LOCAL_DIR "${PROJECT_ROOT}/deps/annlib")
set(ANN_LOCAL_LIBRARY "${ANN_LOCAL_DIR}/lib/libANN.a")
set(ANN_BUILD_TARGET "")
if(EXISTS "${ANN_LOCAL_DIR}/Makefile" AND EXISTS "${ANN_LOCAL_DIR}/include/ANN/ANN.h")
    find_program(ANN_MAKE_EXECUTABLE NAMES gmake make)
    if(NOT ANN_MAKE_EXECUTABLE)
        message(FATAL_ERROR "Cannot build ANN from submodule: 'make' was not found.")
    endif()

    add_custom_command(
        OUTPUT "${ANN_LOCAL_LIBRARY}"
        COMMAND "${ANN_MAKE_EXECUTABLE}" clean
        COMMAND "${ANN_MAKE_EXECUTABLE}" targets
                "ANNLIB = libANN.a"
                "C++ = g++"
                "CFLAGS = -O3 -fPIC"
                "MAKELIB = ar ruv"
                "RANLIB = true"
        WORKING_DIRECTORY "${ANN_LOCAL_DIR}/src"
        COMMENT "Building ANN static library from submodule with -fPIC"
        VERBATIM
    )
    add_custom_target(ReebHanTunAnnBuild DEPENDS "${ANN_LOCAL_LIBRARY}")

    set(ANN_BUILD_TARGET ReebHanTunAnnBuild)
    set(ANN_LIBRARY "${ANN_LOCAL_LIBRARY}")
    set(ANN_INCLUDE_DIR "${ANN_LOCAL_DIR}/include")
    message(STATUS "Using ANN submodule at: ${ANN_LOCAL_DIR}")
    include_directories(${ANN_INCLUDE_DIR})
else()
    # Fallback to system installation
    find_library(ANN_LIBRARY NAMES ann ANN PATHS /usr/lib /usr/lib/x86_64-linux-gnu /usr/local/lib)
    find_path(ANN_INCLUDE_DIR ANN/ANN.h PATHS /usr/include /usr/local/include)
    if(ANN_LIBRARY AND ANN_INCLUDE_DIR)
        message(STATUS "Found ANN library (system): ${ANN_LIBRARY}")
        message(STATUS "Found ANN headers (system): ${ANN_INCLUDE_DIR}")
        include_directories(${ANN_INCLUDE_DIR})
    else()
        message(FATAL_ERROR "ANN library not found. Install the system/Conda ANN package (e.g. libann-dev or conda-forge::ann). Optional local source at ${ANN_LOCAL_DIR} is no longer shipped.")
    endif()
endif()

# Create static library (without main.cpp)
add_library(ReebHanTunLib STATIC
        ComputeReebGraph.cpp
        CycleOptimization/Annotation/AnnotationComputation.cpp
        CycleOptimization/ShortCycles.cpp
        CycleOptimization/canonical_loops.cpp
        CycleOptimization/DijkstraAlgorithm.cpp
        CycleOptimization/edge_annotations_gauss.cpp
        CycleOptimization/FibonacciHeap.cpp
        AdjListGraph.cpp
        ANNSearch.cpp
        FilesOutputForOptimalCycles.cpp
        IntersectionOnMesh.cpp
        InverseLinkNumberMatrix.cpp
        LinkNumberComputing.cpp
        MapLoopsBackToMeshLevelSetAndArc.cpp
        NonOverlappedLevelCycleAndArc.cpp
        NonOverlappingCycles.cpp
        psbmReebGraph.cpp
        ReebGraphPairing.cpp
        SegmentIntersection.cpp
        SegmentPairIntersection.cpp
        SimpleMesh.cpp
        UniqueDirection.cpp
        )

if(ANN_BUILD_TARGET)
    add_dependencies(ReebHanTunLib ${ANN_BUILD_TARGET})
endif()

# Link dependencies to library
target_link_libraries(ReebHanTunLib PUBLIC
        ${Boost_LIBRARIES}
        Boost::program_options
        ${BOOST_SYSTEM_LIB}
        CGAL::CGAL
        ${ANN_LIBRARY})

# Make include directories PUBLIC so other projects can use the headers
# Set to parent directory so <ReebHanTun/SimpleMesh.h> resolves correctly
# Also include current directory for library's own source files
get_filename_component(REEBHANTUN_PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)
target_include_directories(ReebHanTunLib PUBLIC 
    ${REEBHANTUN_PARENT_DIR}      # For <ReebHanTun/SimpleMesh.h> includes
    ${CMAKE_CURRENT_SOURCE_DIR}   # For direct "SimpleMesh.h" includes in library sources
)


# Compile definitions
target_compile_definitions(ReebHanTunLib PUBLIC BOOST_TIMER_ENABLE_DEPRECATED)

# Create an alias for cleaner usage in parent projects
add_library(ReebHanTun::ReebHanTunLib ALIAS ReebHanTunLib)

# Create executable that uses the library
add_executable(ReebHanTun main.cpp)
target_link_libraries(ReebHanTun PRIVATE ReebHanTunLib)