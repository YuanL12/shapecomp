#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "geomp/mesh/Types.h"
#include "geomp/mesh/MeshIO.h"
// #include <basic.hpp>
#include <load.hpp>
// #include <Tutte.h>
// #include <Normalize.h>
#include <pybind11/stl.h>
#include "geomp/torus_comparison/Planar.h"
#include "geomp/torus_comparison/GeodesicPath.h"
#include "geomp/mesh/Mesh.h"
namespace py = pybind11;


// Template function to convert Eigen matrices to a row-major NumPy array
template <typename T>
pybind11::array_t<T> eigen_to_row_major_numpy(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& col_major_matrix) {
    // Convert to row-major layout
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> row_major_matrix = col_major_matrix;

    // Return as NumPy array
    return pybind11::array_t<T>(
        {row_major_matrix.rows(), row_major_matrix.cols()}, // Shape
        {sizeof(T) * row_major_matrix.cols(), sizeof(T)},   // Strides
        row_major_matrix.data()                             // Pointer to data
    );
}

using namespace geomp;

std::shared_ptr<Mesh> loadMesh(const std::string& inputPath, std::optional<bool> normalize = std::nullopt) {
    Model model;
    std::vector<uint8_t> isSurfaceClosed;
    // Load the model from the input path
    loadModel(inputPath, model, isSurfaceClosed, normalize);
    // Validate that there's only one connected component
    if (model.size() != 1) {
        throw std::runtime_error("There are " + std::to_string(model.size()) + " connected components in the input path.");
    }
    // Create a shared pointer to the first mesh in the model
    return std::make_shared<Mesh>(std::move(model[0]));
}

std::shared_ptr<Mesh> loadMesh(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces, std::optional<bool> runNormalize = std::nullopt) {
    Model model;
    loadModel(vertices, faces, model, runNormalize);
    // Validate that there's only one connected component
    if (model.size() != 1) {
        throw std::runtime_error("There are " + std::to_string(model.size()) + " connected components in the input path.");
    }
    // Create a shared pointer to the first mesh in the model
    return std::make_shared<Mesh>(std::move(model[0]));
}

void init_mesh(py::module& m){
    m.def("load_mesh", py::overload_cast<const std::string&, std::optional<bool>>(&loadMesh), "Load mesh from file", 
        py::arg("inputPath"), py::arg("normalize") = std::nullopt);
    m.def("load_mesh", py::overload_cast<const Eigen::MatrixXd&, const Eigen::MatrixXi&, std::optional<bool>>(&loadMesh), 
        "Load mesh from numpy arrays", 
        py::arg("vertices"), py::arg("faces"), py::arg("normalize") = std::nullopt);

    py::class_<Mesh, std::shared_ptr<Mesh>>(m, "Mesh")
        // .def(py::init<const Mesh &>())
        .def(py::init<>())
        .def_property_readonly("eulerCharacteristic", &Mesh::eulerCharacteristic)
        .def("get_vertices_shared", &Mesh::getVertices, "Get the vertices of the mesh, with the data shared with cpp")
        .def("get_faces_shared", &Mesh::getFaces, "Get the faces of the mesh, with the data shared with cpp")
        .def("get_edges", &Mesh::getEdges, "Get the edges of the mesh")
        .def("get_vertices", [](const Mesh& mesh) {
            
            // Call the original method to get the column-major matrix
            Eigen::MatrixXd col_major_matrix = mesh.getVertices();

            // Convert to row-major NumPy array using the utility function
            return eigen_to_row_major_numpy(col_major_matrix);

        }, "Get the vertices of the mesh in row-major order")
        
        .def("get_faces", [](const Mesh& mesh) {
            // Call the original method to get the column-major matrix
            Eigen::MatrixXi col_major_matrix = mesh.getFaces();

            // Convert to row-major NumPy array using the utility function
            return eigen_to_row_major_numpy(col_major_matrix);

        }, "Get the faces of the mesh in row-major order")
        // .def("view_vf", &Mesh::view_vf, "get vertices and faces", py::arg("new_poisition") = false)
        // .def("__str__", &Mesh::sizes_info)
        // .def("cross_edge_ratios", &Mesh::cross_edge_ratios, "cross edge ratios on all edges")
        ;
}

void init_planar_locator(py::module& m){
    py::class_<PlanarLocator>(m, "PlanarLocator")
        .def(py::init<>(), "Construct a default Null PlanarLocator")
        .def(py::init<Mesh&, std::optional<std::string>>(), 
            "Construct a PlanarLocator by passing torus mesh (before cutting) and the Tutte embedding type \ncurrent supported types are: Uniform, coTan, MeanValue, Authalic", 
            py::arg("mesh"), py::arg("Tutte_embedding_type") = std::nullopt)
        .def("constructPlanarLocatorShortest", &PlanarLocator::constructPlanarLocatorShortest, 
            "Construct a PlanarLocator by passing torus mesh (before cutting) and the Tutte embedding type \nshortest path will be used to find the generators", 
            py::arg("mesh"), py::arg("Tutte_embedding_type") = std::nullopt)
        .def("constructPlanarLocatorReebGraph", &PlanarLocator::constructPlanarLocatorReebGraph, 
            "Construct a PlanarLocator by passing torus mesh (before cutting) and the Tutte embedding type \nReeb graph algorithm will be used to find the generators", 
            py::arg("mesh"), py::arg("Tutte_embedding_type") = std::nullopt,
            py::arg("swap_generators") = false)
        .def("constructPlanarLocatorReebGraphOriented", &PlanarLocator::constructPlanarLocatorReebGraphOriented, 
            "Construct a PlanarLocator by passing torus mesh (before cutting) and the Tutte embedding type \nReeb graph algorithm will be used to find the generators, the orientation of the generators will be used to construct the planar locator", 
            py::arg("mesh"), 
            py::arg("Tutte_embedding_type") = std::nullopt,
            py::arg("distinctDirection") = std::nullopt,
            py::arg("swap_generators") = false)
        .def(py::init<Mesh&, const Eigen::MatrixXd&>(), "Construct a PlanarLocator by passing cutted mesh and uv positions")
        .def("get_mesh", &PlanarLocator::get_mesh_ptr, "Get the shared pointer to the subdivided mesh")
        .def("get_identification_map", &PlanarLocator::getIdentificationMap, "Get the identification map")
        .def("get_uv_positions", [](const PlanarLocator& locator) {
            // Call the original method to get the column-major matrix
            Eigen::MatrixXd col_major_matrix = locator.get_uv_positions();

            // Convert to row-major NumPy array using the utility function
            return eigen_to_row_major_numpy(col_major_matrix);

        }, "Get the UV positions of the vertices in row-major order")
        // get the cutted mesh vertices, faces, edges
        .def("get_cutted_mesh_vertices", 
            [](const PlanarLocator& locator){
                auto mesh_ptr = locator.get_cutted_mesh_ptr();
                if (!mesh_ptr) throw std::runtime_error("No cutted mesh in PlanarLocator!");
                return eigen_to_row_major_numpy(mesh_ptr->getVertices());
            }, 
            "Get the cutted mesh vertices")
        .def("get_cutted_mesh_faces", 
            [](const PlanarLocator& locator){
                auto mesh_ptr = locator.get_cutted_mesh_ptr();
                if (!mesh_ptr) throw std::runtime_error("No cutted mesh in PlanarLocator!");
                return eigen_to_row_major_numpy(mesh_ptr->getFaces());
            }, 
            "Get the cutted mesh faces")
        .def("get_cutted_mesh_edges", 
            [](const PlanarLocator& locator){
                auto mesh_ptr = locator.get_cutted_mesh_ptr();
                if (!mesh_ptr) throw std::runtime_error("No cutted mesh in PlanarLocator!");
                return eigen_to_row_major_numpy(mesh_ptr->getEdges());
            }, 
            "Get the cutted mesh edges")
        // get the subdivided mesh vertices, faces, edges
        .def("get_subdivided_mesh_ptr", &PlanarLocator::get_subdivided_mesh_ptr, "Get the shared pointer to the subdivided mesh")
        .def("get_subdivided_mesh_vertices", 
            [](const PlanarLocator& locator){
                auto mesh_ptr = locator.get_subdivided_mesh_ptr();
                if (!mesh_ptr) throw std::runtime_error("No subdivided mesh in PlanarLocator!");
                return eigen_to_row_major_numpy(mesh_ptr->getVertices());
            }, 
            "Get the subdivided mesh vertices")
        .def("get_subdivided_mesh_faces", 
            [](const PlanarLocator& locator){
                auto mesh_ptr = locator.get_subdivided_mesh_ptr();
                if (!mesh_ptr) throw std::runtime_error("No subdivided mesh in PlanarLocator!");
                return eigen_to_row_major_numpy(mesh_ptr->getFaces());
            }, 
            "Get the subdivided mesh faces")
        .def("get_subdivided_mesh_edges", 
            [](const PlanarLocator& locator){
                auto mesh_ptr = locator.get_subdivided_mesh_ptr();
                if (!mesh_ptr) throw std::runtime_error("No subdivided mesh in PlanarLocator!");
                return eigen_to_row_major_numpy(mesh_ptr->getEdges());
            }, 
            "Get the subdivided mesh edges")
        // find the barycentric coordinates of the triangle containing the point p
        .def("find_bary_coords", 
            // (std::pair<int, Eigen::Vector3d> (geomp::PlanarLocator::*)(const std::pair<double, double>&) const) &geomp::PlanarLocator::findBaryCoords, 
            py::overload_cast<const std::pair<double, double>&>(&PlanarLocator::findBaryCoords, py::const_),
            "Find barycentric coordinates of a point", 
            py::arg("P_uv"))

        .def("find_bary_coords_serial", 
            [](const PlanarLocator& locator, const Eigen::Vector2d& point) {
                return locator.findBaryCoords(point);
            }, 
            "Find barycentric coordinates of a point in serial", 
            py::arg("point"))


        // (Parallel) find the barycentric coordinates of triangles containing given points
        .def("find_bary_coords_parallel", 
            [](const PlanarLocator& locator, const Eigen::MatrixXd& points) {
                return locator.findBaryCoordsParallel(points);
            }, 
            "Find barycentric coordinates of multiple points in parallel", 
            py::arg("points"))
        
        // Compute geodesic path between two UV points
        .def("compute_geodesic_path_with_homotopy_class", 
            [](const PlanarLocator& locator, 
               const std::pair<double, double>& a,
               const std::pair<double, double>& b) {
                Eigen::MatrixXd pathMatrix = computeGeodesicPathFromPlanarLocator(
                    locator, a, b);
                return eigen_to_row_major_numpy(pathMatrix);
            },
            "Compute geodesic path on torus T between two UV points a and b. "
            "Uses the subdivided mesh from the PlanarLocator.",
            py::arg("a"),          // Start point in UV coordinates (can be in [-1,2]^2)
            py::arg("b"))          // End point in UV coordinates (can be in [-1,2]^2)
        .def("compute_dijkstra_f3_path_3d", 
            [](const PlanarLocator& locator, 
               const std::pair<double, double>& a,
               const std::pair<double, double>& b) {
                Eigen::MatrixXd pathMatrix = computeDijkstraF3Path3D(
                    locator, a, b);
                return eigen_to_row_major_numpy(pathMatrix);
            },
            "Compute the dijkstraF3 path in 3D coordinates (debug: the input halfedge path (before shortening) as 3D coordinates). "
            "Returns the vertex path along the halfedges as a 3D numpy array.",
            py::arg("a"),          // Start point in UV coordinates (can be in [-1,2]^2)
            py::arg("b"))          // End point in UV coordinates (can be in [-1,2]^2)
        .def("compute_dijkstra_f3_path_2d", 
            [](const PlanarLocator& locator, 
               const std::pair<double, double>& a,
               const std::pair<double, double>& b) {
                Eigen::MatrixXd path2D = computeDijkstraF3Path2D(
                    locator, a, b);
                return eigen_to_row_major_numpy(path2D);
            },
            "Debug version: Compute the dijkstraF3 path in 2D UV coordinates. "
            "Returns the F3 vertex path projected to 2D UV space as a numpy array (Nx2).",
            py::arg("a"),          // Start point in UV coordinates (can be in [-1,2]^2)
            py::arg("b"))          // End point in UV coordinates (can be in [-1,2]^2)
            
        .def("get_generator1", &PlanarLocator::getGenerator1, "Get the vertex indices (in cutted mesh) of the first generator path")
        .def("get_generator2", &PlanarLocator::getGenerator2, "Get the vertex indices (in cutted mesh) of the second generator path")

        .def("get_generator_paths_3d", 
            [](const PlanarLocator& locator) {
                auto [path1, path2] = locator.getGeneratorPaths3D();
                return std::make_pair(
                    eigen_to_row_major_numpy(path1),
                    eigen_to_row_major_numpy(path2)
                );
            },
            "Get the 3D points along the two generator paths. "
            "Returns a tuple of two numpy arrays, each Nx3 representing the 3D points along a generator path.")

    ;
}


void init_meshio(py::module& m){
    m.def("write_mesh_into_obj", &MeshIO::writeMeshOBJ, "Write mesh to obj file", py::arg("fileName"), py::arg("mesh"));
    m.def("subdivide_mesh",  py::overload_cast<const Mesh&>(&MeshIO::BarycenterDivideMesh), "Subdivide mesh", py::arg("mesh"));
}

// void init_normalize(py::module& m){
//     py::class_<Normalize>(m, "Normalize")
//         .def_static("normalize", &Normalize::normalize, "Normalzie through `Moebius Registration` by A. Baden, K. Crane, and M. Kazhdan",  py::arg("mesh"), py::arg("verbose") = false)
//         .def_static("inversion2zero", &Normalize::inversion2zero, "Normalzie through applying double inversion")
//         .def_static("find_center_of_mass", &Normalize::find_center_of_mass)
//         ;
// }