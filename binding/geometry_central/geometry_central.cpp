/*
Python Binding from geometrycentral to Python. It supports
    - Geodesic Path
    - Flip Geodesic Path
    - Convert Vector3 to numpy array(dense matrix)
*/

#include "geometry_central.hpp"

void init_geodesic(py::module &m) {
    py::class_<geomp::GeodesicsManager>(m, "GeodesicsManager")
        .def(py::init<GC::DenseMatrix<double>, GC::DenseMatrix<int64_t>>(),
            "Load mesh from vertices and faces matrix",
            py::arg("verts"),
            py::arg("faces"))
        .def(py::init<std::string>(), "Load mesh from file path", py::arg("mesh_path"))

        .def("find_exact_geodesic_path",
            py::overload_cast<int64_t, int64_t, const std::vector<double> &, const std::vector<double> &>(
                &geomp::GeodesicsManager::find_exact_geodesic_path),
            py::arg("startFaceIdx"),
            py::arg("endFaceIdx"),
            py::arg("startFaceCoords"),
            py::arg("endFaceCoords"),
            "Find the geodesic path between two points")

        .def("find_exact_geodesic_path",
            py::overload_cast<std::pair<int64_t, const std::vector<double> &>,
                std::pair<int64_t, const std::vector<double> &>>(&geomp::GeodesicsManager::find_exact_geodesic_path),
            py::arg("start_face_idx_coords_pair"),
            py::arg("end_face_idx_coords_pair"),
            "Find the geodesic path between two points")

        .def("find_surface_points_on_exact_geodesic_path",
            &geomp::GeodesicsManager::find_surface_points_on_exact_geodesic_path,
            "Find the surface points on the exact geodesic path",
            py::arg("startFaceIdx"),
            py::arg("endFaceIdx"),
            py::arg("startFaceCoords"),
            py::arg("endFaceCoords"))

        .def("find_flip_geodesic_path",
            py::overload_cast<int64_t, int64_t>(&geomp::GeodesicsManager::find_flip_geodesic_path),
            "Find the flip geodesic path between two vertices",
            py::arg("startVertexIdx"),
            py::arg("endVertexIdx"))

        .def("find_flip_geodesic_path",
            py::overload_cast<int64_t, int64_t, const std::vector<double> &, const std::vector<double> &>(
                &geomp::GeodesicsManager::find_flip_geodesic_path),
            "Find the flip geodesic path between two face points",
            py::arg("startFaceIdx"),
            py::arg("endFaceIdx"),
            py::arg("startFaceCoords"),
            py::arg("endFaceCoords"));

    // Surface points on a flip-geodesic path is no longer the surface points on the original mesh (due to the mesh
    // modification), so we do not bind this function. .def("find_surface_points_on_flip_geodesic_path",
    //     &GeodesicsManager::find_surface_points_on_flip_geodesic_path,
    //     "Find the surface points on the flip geodesic path",
    //     py::arg("startFaceIdx"),
    //     py::arg("endFaceIdx"),
    //     py::arg("startFaceCoords"),
    //     py::arg("endFaceCoords"));
}
