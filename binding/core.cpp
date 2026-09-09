#include "core.hpp"

PYBIND11_MODULE(_shapecomp, m) {
    m.doc() = "ShapeComp: tools for comparing shapes and surfaces";

    // simple test function to check if the binding is working add
    m.def("add", [](int a, int b) { return a + b; });

    // Mesh
    init_mesh(m);

    // PlanarLocator
    init_planar_locator(m);
    
    // Compute geodesic by Geometry-Central
    init_geodesic(m);

    // Parallel computation supprted by Taskflow
    init_parallel_taskflow(m);

    // MeshIO
    init_meshio(m);
    
    // // Mapping
    // m.def("tutte_embedding", &find_tutte_embedding);
    // init_normalize(m);
}