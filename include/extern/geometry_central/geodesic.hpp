#include <cmath>
#include <iostream>
#include <iterator>
#include <optional>
#include <vector>

// Eigen
#include "Eigen/Dense"

// geometrycentral
#include "geometrycentral/surface/exact_geodesics.h"
#include "geometrycentral/surface/flip_geodesics.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

// name space alias for geometrycentral
namespace GC = geometrycentral;
namespace GCS = geometrycentral::surface;

namespace geomp {

class GeodesicsManager {
   public:
    GeodesicsManager(GC::DenseMatrix<double> verts, GC::DenseMatrix<int64_t> faces);
    GeodesicsManager(std::string mesh_path);

    // Compute the geodesic path between two points, return a n*3 matrix
    GC::DenseMatrix<double> find_exact_geodesic_path(int64_t startFaceIdx,
        int64_t endFaceIdx,
        const std::vector<double> &startFaceCoords,
        const std::vector<double> &endFaceCoords);

    GC::DenseMatrix<double> find_exact_geodesic_path(
        std::pair<int64_t, const std::vector<double> &> start_face_idx_coords_pair,
        std::pair<int64_t, const std::vector<double> &> end_face_idx_coords_pair);

    // Compute the geodesic path between two points, return the path as a vector
    // of surface points
    std::vector<std::pair<int64_t, std::vector<double>>> find_surface_points_on_exact_geodesic_path(
        int64_t startFaceIdx,
        int64_t endFaceIdx,
        const std::vector<double> &startFaceCoords,
        const std::vector<double> &endFaceCoords);

    // Compute the geodesic path between two vertices, return the path as a
    // vector of surface points
    GC::DenseMatrix<double> find_flip_geodesic_path(int64_t startVertexIdx, int64_t endVertexIdx);

    // Compute the geodesic path between two points, return a n*3 matrix
    GC::DenseMatrix<double> find_flip_geodesic_path(int64_t startFaceIdx,
        int64_t endFaceIdx,
        const std::vector<double> &startFaceCoords,
        const std::vector<double> &endFaceCoords);

    // Compute the geodesic path between two points, return the path as a vector
    // of surface points
    std::vector<std::pair<int64_t, std::vector<double>>> find_surface_points_on_flip_geodesic_path(int64_t startFaceIdx,
        int64_t endFaceIdx,
        const std::vector<double> &startFaceCoords,
        const std::vector<double> &endFaceCoords);

    static GCS::SurfacePoint convert_face_point_to_specific(const GCS::SurfacePoint &surface_point);

    // getter functions (return raw pointers is fine, since unique_ptr will hanlde the memory management)
    inline GCS::ManifoldSurfaceMesh *get_mesh() { return mesh.get(); };
    inline GCS::VertexPositionGeometry *get_geometry() { return geometry.get(); };
    
    // Compute the geodesic path between two points.
    // It cannot handle edge case (two points on the same edge or the same face) directly.
    // @param start_surface_point: the start surface point
    // @param end_surface_point: the end surface point
    // @return the path as a vector of GCS::SurfacePoint
    std::vector<GCS::SurfacePoint> find_flip_geodesic_path_GC(GCS::SurfacePoint start_surface_point,
        GCS::SurfacePoint end_surface_point);

    // Compute the geodesic path between two points.
    // It can handle edge case (two points on the same edge or the same face) directly.
    // @param startFaceIdx: the start face index
    // @param endFaceIdx: the end face index
    // @param startFaceCoords: the start face coordinates
    // @param endFaceCoords: the end face coordinates
    // @return the path as a vector of GCS::SurfacePoint
    std::vector<GCS::SurfacePoint> find_flip_geodesic_path_GC(int64_t startFaceIdx,
        int64_t endFaceIdx,
        const std::vector<double> &startFaceCoords,
        const std::vector<double> &endFaceCoords);
    
    // Check if two surface points are on the same face
    // @param surface_point_1: the first surface point
    // @param surface_point_2: the second surface point
    // @return true if they are on the same face, false otherwise
    static bool is_on_same_face(const GCS::SurfacePoint &surface_point_1, const GCS::SurfacePoint &surface_point_2);
   private:
    std::unique_ptr<GCS::ManifoldSurfaceMesh> mesh;
    std::unique_ptr<GCS::VertexPositionGeometry> geometry;

    // Helper function to convert a surface point to a face index and
    // barycentric coordinates
    std::pair<int64_t, std::vector<double>> surface_point_to_face_idx_coords(const GCS::SurfacePoint &surface_point);

    // Compute the geodesic path between two points, return the path as a vector
    // of GCS::SurfacePoint
    std::vector<GCS::SurfacePoint> find_exact_geodesic_path_GC(int64_t startFaceIdx,
        int64_t endFaceIdx,
        const std::vector<double> &startFaceCoords,
        const std::vector<double> &endFaceCoords);

    // Shared handling for degenerate configurations that flip / exact geodesic backends
    // handle poorly (same triangle, same edge, face–edge on one triangle). Returns a
    // 2-point path when resolved here; otherwise fills converted endpoints for the
    // full algorithm and returns std::nullopt.
    std::optional<std::vector<GCS::SurfacePoint>> preprocess_surface_points(int64_t startFaceIdx,
        int64_t endFaceIdx,
        const std::vector<double> &startFaceCoords,
        const std::vector<double> &endFaceCoords,
        GCS::SurfacePoint &out_start_converted,
        GCS::SurfacePoint &out_end_converted);

    // Convert a surface point to a vertex
    // note: this function will modify the mesh and geometry if the surface
    // point is an edge or face
    GCS::Vertex convert_surface_point_to_vertex(GCS::SurfacePoint surface_point);

    // Convert a vector of surface points to a matrix of surface points
    GC::DenseMatrix<double> surface_points_to_dense_matrix(const std::vector<GCS::SurfacePoint> &surface_points);

    // Convert a vector of surface points to a vector of face index and
    // barycentric coordinates
    std::vector<std::pair<int64_t, std::vector<double>>> surface_points_to_face_idx_coords_list(
        const std::vector<GCS::SurfacePoint> &surface_points);
};

}  // namespace geomp