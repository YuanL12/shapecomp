# Find the module development requirements (requires FindPython from 3.17 or
# scikit-build-core's built-in backport)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# Copy the shapecomp package into the build directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/shapecomp DESTINATION ${CMAKE_BINARY_DIR})

# Add a library using FindPython's tooling (pybind11 also provides a helper like this)
python_add_library(_shapecomp MODULE 
    core.cpp geomp.cpp 
    parallel_taskflow.cpp 
    geometry_central/geometry_central.cpp
    WITH_SOABI)
target_link_libraries(_shapecomp PRIVATE pybind11::headers Eigen3::Eigen gp geometry-central taskflow)

# Set the output directory for the .so file
set_target_properties(_shapecomp PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/shapecomp)

# The install directory is the output (wheel) directory
install(TARGETS _shapecomp DESTINATION shapecomp)
install(FILES ${CMAKE_BINARY_DIR}/shapecomp/__init__.py DESTINATION shapecomp)
