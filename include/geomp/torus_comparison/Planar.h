#pragma once
#include <Eigen/Core>
#include <memory>
#include <Eigen/Sparse>
#include "geomp/mesh/Types.h"
#include "geomp/mesh/MeshData.h"
#include "geomp/generators/Generators.h"
#include "geomp/generators/MeshSurgery.h"
#include "geomp/embedding/torus/Tutte.h"
#include "geomp/mesh/MeshIO.h"
// Include Taskflow for parallelization of findBaryCoords
#include <taskflow/taskflow.hpp>  
#include <taskflow/core/flow_builder.hpp>
#include <taskflow/algorithm/sort.hpp>
#include <taskflow/algorithm/for_each.hpp>
// #include <torch/extension.h>

namespace geomp
{
/*
    PlanarLocator class
	Given a planar mesh (a triangulation the unit square) and a point p in R^2
	find the bary coordinates of the traingle containing p. 
	Specifically, the mesh connectivity is given in Mesh subdivided_mesh instance
	the uv positions of the vertices are given in VertexData<std::pair<double, double>> uv_postiosns
	the point p is given as a pair<double, double> p
	the ouput out should contain 1. the face index 2. bary_coords as a tuple
*/
class PlanarLocator
{
    using UV_TYPE = VertexData<std::pair<double, double>>;
private:
    Mesh* mesh_ptr; // the original mesh (PlanarLocator does not own it)
    std::shared_ptr<Mesh> subdivided_mesh_ptr; // the subdivided mesh
    std::shared_ptr<Mesh> cutted_mesh_ptr; // the cutted mesh
    std::shared_ptr<UV_TYPE> uv_positions_ptr; // the uv positions of the vertices of the cutted mesh
    std::vector<int> generator1; // vertex indices of the first generator path (in cutted mesh)
    std::vector<int> generator2; // vertex indices of the second generator path (in cutted mesh)
    int get_Tutte_embedding_type(std::string type) const; // get the type number of the Tutte embedding type

public:
    //Default constructor
    PlanarLocator():
        mesh_ptr(nullptr),
        subdivided_mesh_ptr(nullptr),
        cutted_mesh_ptr(nullptr),
        uv_positions_ptr(nullptr),
        generator1(std::vector<int>()),
        generator2(std::vector<int>()) {};

    // constructor by passing torus mesh, a subdivided mesh, its cutted version and uv_positions will be stored 
    PlanarLocator(Mesh& mesh, std::optional<std::string> Tutte_embedding_type = std::nullopt);

    // constructor by passing torus mesh, a subdivided mesh, its cutted version and uv_positions will be stored 
    void constructPlanarLocatorShortest(Mesh& mesh, std::optional<std::string> Tutte_embedding_type = std::nullopt);

    // constructor by passing torus mesh, a subdivided mesh, its cutted version and uv_positions will be stored 
    void constructPlanarLocatorReebGraph(Mesh& mesh,
        std::optional<std::string> Tutte_embedding_type = std::nullopt,
        bool swap_generators = false);

    // constructor by passing torus mesh, a subdivided mesh, its cutted version and uv_positions will be stored 
    void constructPlanarLocatorReebGraphOriented(Mesh& mesh, 
        std::optional<std::string> Tutte_embedding_type = std::nullopt,
        std::optional<std::vector<double>> distinctDirection = std::nullopt,
        bool swap_generators = false);
    
    // Getter for generator1 and generator2
    inline std::vector<int> getGenerator1() const {return generator1;}
    inline std::vector<int> getGenerator2() const {return generator2;}

    // constructor by passing a cutted mesh and uv_positions (subdvision will not be done here)
    PlanarLocator(Mesh& cutted_mesh_, const UV_TYPE& uv_positions_):
        cutted_mesh_ptr(std::make_shared<Mesh>(cutted_mesh_)), 
        uv_positions_ptr(std::make_shared<UV_TYPE>(uv_positions_)),
        generator1(),
        generator2() {
    };

    // constructor by passing a cutted mesh and uv_positions in Eigen::MatrixXd
    PlanarLocator(Mesh& cutted_mesh_, const Eigen::MatrixXd& uv_positions_matrix):
        cutted_mesh_ptr(std::make_shared<Mesh>(cutted_mesh_)),
        generator1(),
        generator2() {
        uv_positions_ptr = std::make_shared<UV_TYPE>(cutted_mesh_, std::make_pair(0.0, 0.0));
        for (int i = 0; i < uv_positions_matrix.rows(); i++) {
            VertexCIter v = cutted_mesh_.vertices.begin() + i;
            (*uv_positions_ptr)[v] = std::make_pair(uv_positions_matrix(i, 0), uv_positions_matrix(i, 1));
        }
    };

    // copy constructor
    PlanarLocator(const PlanarLocator& other): 
        mesh_ptr(other.mesh_ptr),
        cutted_mesh_ptr(other.cutted_mesh_ptr), 
        subdivided_mesh_ptr(other.subdivided_mesh_ptr),
        uv_positions_ptr(other.uv_positions_ptr),
        generator1(other.generator1),
        generator2(other.generator2) {};

    // find the barycentric coordinates of the triangle containing the point p
    std::pair<int, Eigen::Vector3d> findBaryCoords(const geomp::Vector& p) const;
    std::pair<int, Eigen::Vector3d> findBaryCoords(const std::pair<double, double>& P_uv) const;
    std::pair<int, Eigen::Vector3d> findBaryCoords(const Eigen::Vector2d& point) const;
    
    // (Parallel) find the barycentric coordinates of triangles containing given points
    // each row of points is a point in R^2
    std::vector<std::pair<int, Eigen::Vector3d>> findBaryCoordsParallel(const Eigen::MatrixXd& points) const;
    // each element of points is a point in R^2
    std::vector<std::pair<int, Eigen::Vector3d>> findBaryCoordsParallel(const std::vector<std::pair<double, double>>& points) const;


    // Get the identification map of the cutted mesh, i.e. vertex index -> its copy vertex indices in the cutted mesh
    std::unordered_map<int, std::vector<int>> getIdentificationMap() const;

    // Getter for mesh_ptr
    inline Mesh* get_mesh_ptr() const {
        return mesh_ptr;
    }

    // Getter for subdivided_mesh_ptr
    inline std::shared_ptr<Mesh> get_subdivided_mesh_ptr() const {
        return subdivided_mesh_ptr;
    }

    // Getter for cutted_mesh_ptr
    inline std::shared_ptr<Mesh> get_cutted_mesh_ptr() const {
        return cutted_mesh_ptr;
    }

    // Getter for uv_positions_ptr
    inline std::shared_ptr<UV_TYPE> get_uv_positions_ptr() const {
        return uv_positions_ptr;
    }    

    // Getter for uv_positions_ptr, return Eigen::MatrixXd
    inline Eigen::MatrixXd get_uv_positions() const {
        // pointer safe check
        if (!cutted_mesh_ptr || !uv_positions_ptr) {
            throw std::runtime_error("cutted_mesh or uv_positions is null.");
        }
        // get the uv positions of the vertices
        Eigen::MatrixXd uvs_vect(cutted_mesh_ptr->vertices.size(), 2);
        for (VertexCIter v = cutted_mesh_ptr->vertices.begin(); v != cutted_mesh_ptr->vertices.end(); v++) {
            // get the uv position of the vertex
            const auto& uv_positions = *uv_positions_ptr;
            uvs_vect(v->index, 0) = uv_positions[v].first;
            uvs_vect(v->index, 1) = uv_positions[v].second;
        }
        return uvs_vect;
    }

    // Debug functions
    // Get 3D points along the two generator paths
    // Returns a pair of matrices, each Nx3 representing the 3D points along a generator path
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> getGeneratorPaths3D() const;

    // Get the two generators from the circular linked list of vertices on the boundary
    // reverse the linked list if the starting edge is on the second generator
    void getGeneratorsFromBoundary(CircularLinkedList& linked_vertices_bd, bool swap_generators = false);
    
};

} // namespace geomp


