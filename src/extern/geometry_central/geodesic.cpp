#include "extern/geometry_central/geodesic.hpp"

#include <stdexcept>

namespace geomp {

GeodesicsManager::GeodesicsManager(GC::DenseMatrix<double> verts, GC::DenseMatrix<int64_t> faces) {
    // Construct the internal mesh and geometry
    mesh = std::make_unique<GCS::ManifoldSurfaceMesh>(faces);
    geometry = std::make_unique<GCS::VertexPositionGeometry>(*mesh);

    // Load the vertex positions
    for (size_t i = 0; i < mesh->nVertices(); i++) {
        for (size_t j = 0; j < 3; j++) {
            geometry->inputVertexPositions[i][j] = verts(i, j);
        }
    }
};

GeodesicsManager::GeodesicsManager(std::string mesh_path) {
    // Load mesh from file path
    std::tie(mesh, geometry) = GCS::readManifoldSurfaceMesh(mesh_path);
};

// Convert a face-type SurfacePoint to a more specific type (vertex/edge) when possible,
// following the same corner/edge ordering used by SurfacePoint::inSomeFace().
// - If one barycentric coord is ~1, returns a Vertex-type point at that corner.
// - If one barycentric coord is ~0 (and others sum to ~1), returns an Edge-type point
//   on the opposite edge with t consistent with inSomeFace() mappings.
// - Otherwise, returns the original Face-type point.
GCS::SurfacePoint GeodesicsManager::convert_face_point_to_specific(const GCS::SurfacePoint &facePoint) {
    if (facePoint.type != GCS::SurfacePointType::Face) return facePoint;  // already a specific type

    const double eps = 1e-12;

    // Normalize barycentric if very slightly off
    GC::Vector3 bc = facePoint.faceCoords;
    double s = bc.x + bc.y + bc.z;
    if (std::abs(s - 1.0) > eps && s > eps) {
        bc = {bc.x / s, bc.y / s, bc.z / s};
    }

    // Clamp tiny negatives/overflows
    auto clamp01 = [](double v) { return v < 0 ? 0. : (v > 1 ? 1. : v); };
    bc = {clamp01(bc.x), clamp01(bc.y), clamp01(bc.z)};

    GCS::Halfedge he0 = facePoint.face.halfedge();
    GCS::Halfedge he1 = he0.next();
    GCS::Halfedge he2 = he1.next();

    // Vertex cases: (1,0,0), (0,1,0), (0,0,1)
    if (bc.x >= 1.0 - eps) return GCS::SurfacePoint(he0.vertex());
    if (bc.y >= 1.0 - eps) return GCS::SurfacePoint(he1.vertex());
    if (bc.z >= 1.0 - eps) return GCS::SurfacePoint(he2.vertex());

    // Edge cases (using inverse of inSomeFace() mappings):
    if (bc.z <= eps) {
        // Edge corresponding to he0: {1 - t, t, 0}
        GCS::Edge e = he0.edge();
        // he0 is a halfedge of e. Check if he0 is e.halfedge() or its twin
        // SurfacePoint(Edge, t) uses t relative to e.halfedge()
        if (he0 == e.halfedge()) {
            return GCS::SurfacePoint(e, bc.y);
        } else {
            // he0 is the twin, so we need to flip t
            return GCS::SurfacePoint(e, 1.0 - bc.y);
        }
    }
    if (bc.y <= eps) {
        // Edge corresponding to he2: {t, 0, 1 - t}
        GCS::Edge e = he2.edge();
        // he2 is a halfedge of e. Check if he2 is e.halfedge() or its twin
        if (he2 == e.halfedge()) {
            return GCS::SurfacePoint(e, bc.x);
        } else {
            // he2 is the twin, so we need to flip t
            return GCS::SurfacePoint(e, 1.0 - bc.x);
        }
    }
    if (bc.x <= eps) {
        // Edge corresponding to he1: {0, 1 - t, t}
        GCS::Edge e = he1.edge();
        // he1 is a halfedge of e. Check if he1 is e.halfedge() or its twin
        if (he1 == e.halfedge()) {
            return GCS::SurfacePoint(e, bc.z);
        } else {
            // he1 is the twin, so we need to flip t
            return GCS::SurfacePoint(e, 1.0 - bc.z);
        }
    }

    // Remain a face point
    return facePoint;
}

// Helper function to convert a surface point to a face index and barycentric coordinates
std::pair<int64_t, std::vector<double>> GeodesicsManager::surface_point_to_face_idx_coords(
    const GCS::SurfacePoint &surface_point) {
    // get all face indices
    const auto &face_indices = mesh->getFaceIndices();
    auto in_some_face_point = surface_point.inSomeFace();
    const auto &face_coords = in_some_face_point.faceCoords;
    std::vector<double> face_coords_vec{face_coords.x, face_coords.y, face_coords.z};
    return std::make_pair(face_indices[in_some_face_point.face], face_coords_vec);
}

std::vector<std::pair<int64_t, std::vector<double>>> GeodesicsManager::surface_points_to_face_idx_coords_list(
    const std::vector<GCS::SurfacePoint> &surface_points) {
    std::vector<std::pair<int64_t, std::vector<double>>> out;
    for (const auto &surface_point : surface_points) {
        auto face_idx_coords = surface_point_to_face_idx_coords(surface_point);
        out.push_back(face_idx_coords);
    }
    return out;
}

GC::DenseMatrix<double> GeodesicsManager::surface_points_to_dense_matrix(
    const std::vector<GCS::SurfacePoint> &surface_points) {
    GC::DenseMatrix<double> out(surface_points.size(), 3);
    geometry->requireVertexPositions();
    for (size_t i = 0; i < surface_points.size(); i++) {
        // Use vertexPositions (current positions) instead of inputVertexPositions
        // to handle cases where vertices have been inserted/modified
        auto pt_post = surface_points[i].interpolate(geometry->vertexPositions);
        for (size_t j = 0; j < 3; j++) {
            out(i, j) = pt_post[j];
        }
    }
    return out;
}

GC::DenseMatrix<double> GeodesicsManager::find_flip_geodesic_path(int64_t startVertexIdx, int64_t endVertexIdx) {
    // Create a path network as a Dijkstra path between endpoints
    std::unique_ptr<GCS::FlipEdgeNetwork> edgeNetwork;
    GCS::Vertex vStart = mesh->vertex(startVertexIdx);
    GCS::Vertex vEnd = mesh->vertex(endVertexIdx);
    edgeNetwork = GCS::FlipEdgeNetwork::constructFromDijkstraPath(*mesh, *geometry, vStart, vEnd);

    // Make the path a geodesic
    edgeNetwork->iterativeShorten();

    // Extract the result as a polyline along the surface
    edgeNetwork->posGeom = geometry.get();
    std::vector<std::vector<GCS::SurfacePoint>> polyline = edgeNetwork->getPathPolyline();

    // Assert that there is one polyline
    if (polyline.size() != 1) {
        throw std::runtime_error("Expected one polyline, but got " + std::to_string(polyline.size()));
    }

    // Convert the polyline to a matrix of surface points
    return surface_points_to_dense_matrix(polyline[0]);
}

GCS::Vertex GeodesicsManager::convert_surface_point_to_vertex(GCS::SurfacePoint surface_point) {
    // depending on the type of the surface point, convert it to a vertex index
    if (surface_point.type == GCS::SurfacePointType::Vertex) {
        // If the surface point is already a vertex, just return its index
        return surface_point.vertex;
    } else if (surface_point.type == GCS::SurfacePointType::Edge) {
        // compute the 3D position of the surface point and assign it to the new vertex
        geometry->requireVertexPositions();
        GC::Vector3 newPosition = surface_point.interpolate(geometry->vertexPositions);

        // insert a vertex along the edge at the specified position
        GCS::Halfedge newHe = mesh->splitEdgeTriangular(surface_point.edge);
        GCS::Vertex newVertex = newHe.vertex();

        // assign it to the new vertex
        geometry->vertexPositions[newVertex] = newPosition;
        return newVertex;
    } else if (surface_point.type == GCS::SurfacePointType::Face) {
        // compute the 3D position of the surface point and assign it to the new vertex
        geometry->requireVertexPositions();
        GC::Vector3 newPosition = surface_point.interpolate(geometry->vertexPositions);

        // insert a vertex in the face (this automatically triangulates the face)
        GCS::Vertex newVertex = mesh->insertVertex(surface_point.face);

        // assign the 3D position to the new vertex
        geometry->vertexPositions[newVertex] = newPosition;
        return newVertex;
    }
    throw std::runtime_error("Unsupported surface point type or no type is provided");
}

std::vector<GCS::SurfacePoint> GeodesicsManager::find_flip_geodesic_path_GC(GCS::SurfacePoint start_surface_point,
    GCS::SurfacePoint end_surface_point) {
    // convert the surface points to vertices
    GCS::Vertex start_vertex = convert_surface_point_to_vertex(start_surface_point);
    GCS::Vertex end_vertex = convert_surface_point_to_vertex(end_surface_point);

    // compress the mesh
    mesh->compress();

    // Create a path network as a Dijkstra path between endpoints
    std::unique_ptr<GCS::FlipEdgeNetwork> edgeNetwork;
    edgeNetwork = GCS::FlipEdgeNetwork::constructFromDijkstraPath(*mesh, *geometry, start_vertex, end_vertex);

    // Make the path a geodesic
    edgeNetwork->iterativeShorten();

    // Extract the result as a polyline along the surface
    edgeNetwork->posGeom = geometry.get();
    std::vector<std::vector<GCS::SurfacePoint>> polyline = edgeNetwork->getPathPolyline();

    // Assert that there is one polyline
    if (polyline.size() != 1) {
        throw std::runtime_error("Expected one polyline, but got " + std::to_string(polyline.size()));
    }
    return polyline[0];
}

bool GeodesicsManager::is_on_same_face(const GCS::SurfacePoint &surface_point_1, const GCS::SurfacePoint &surface_point_2) {
    // check if the two surface points are on the same face
    return surface_point_1.face == surface_point_2.face;
}

std::optional<std::vector<GCS::SurfacePoint>> GeodesicsManager::preprocess_surface_points(int64_t startFaceIdx,
    int64_t endFaceIdx,
    const std::vector<double> &startFaceCoords,
    const std::vector<double> &endFaceCoords,
    GCS::SurfacePoint &out_start_converted,
    GCS::SurfacePoint &out_end_converted) {
    // check if the input is valid
    if (startFaceCoords.size() != 3 || endFaceCoords.size() != 3) {
        throw std::invalid_argument("Face barycentric coordinates must be a vector of 3 elements");
    }
    GCS::SurfacePoint start_surface_point = GCS::SurfacePoint(mesh->face(startFaceIdx),
        GC::Vector3{startFaceCoords[0], startFaceCoords[1], startFaceCoords[2]});
    GCS::SurfacePoint end_surface_point =
        GCS::SurfacePoint(mesh->face(endFaceIdx), GC::Vector3{endFaceCoords[0], endFaceCoords[1], endFaceCoords[2]});

    // edge cases
    // 1. Same triangle: straight segment in the same face.
    if (startFaceIdx == endFaceIdx) {
        return std::vector<GCS::SurfacePoint>{start_surface_point, end_surface_point};
    }
    // 2. if the start and end face are on the same edge, return the path on the edge
    auto start_converted = convert_face_point_to_specific(start_surface_point);
    auto end_converted = convert_face_point_to_specific(end_surface_point);
    out_start_converted = start_converted;
    out_end_converted = end_converted;
    if (start_converted.type == GCS::SurfacePointType::Edge && end_converted.type == GCS::SurfacePointType::Edge) {
        if (start_converted.edge == end_converted.edge) {
            return std::vector<GCS::SurfacePoint>{start_surface_point, end_surface_point};
        }
    }

    // 3. One point in face interior, the other on an edge of that face (or its twin face).
    bool on_same_face = false;
    if (start_converted.type == GCS::SurfacePointType::Edge && end_converted.type == GCS::SurfacePointType::Face) {
        // start point is on an edge, check its face and its twin face is the same as the end point
        auto start_half_edge = start_converted.edge.halfedge();
        on_same_face = (start_half_edge.face() == end_converted.face);
        on_same_face = on_same_face || (start_half_edge.twin().face() == end_converted.face);
    }
    if (start_converted.type == GCS::SurfacePointType::Face && end_converted.type == GCS::SurfacePointType::Edge) {
        // end point is on an edge, check its face and its twin face is the same as the start point
        auto end_half_edge = end_converted.edge.halfedge();
        on_same_face = (end_half_edge.face() == start_converted.face);
        on_same_face = on_same_face || (end_half_edge.twin().face() == start_converted.face);
    }
    if (on_same_face) {
        return std::vector<GCS::SurfacePoint>{start_surface_point, end_surface_point};
    }

    return std::nullopt;
}

std::vector<GCS::SurfacePoint> GeodesicsManager::find_flip_geodesic_path_GC(int64_t startFaceIdx,
    int64_t endFaceIdx,
    const std::vector<double> &startFaceCoords,
    const std::vector<double> &endFaceCoords) {
    GCS::SurfacePoint start_converted, end_converted;
    auto early = preprocess_surface_points(
        startFaceIdx, endFaceIdx, startFaceCoords, endFaceCoords, start_converted, end_converted);
    if (early.has_value()) {
        return *early;
    }
    return find_flip_geodesic_path_GC(start_converted, end_converted);
}

GC::DenseMatrix<double> GeodesicsManager::find_flip_geodesic_path(int64_t startFaceIdx,
    int64_t endFaceIdx,
    const std::vector<double> &startFaceCoords,
    const std::vector<double> &endFaceCoords) {
    std::vector<GCS::SurfacePoint> path =
        find_flip_geodesic_path_GC(startFaceIdx, endFaceIdx, startFaceCoords, endFaceCoords);
    return surface_points_to_dense_matrix(path);
}

std::vector<std::pair<int64_t, std::vector<double>>> GeodesicsManager::find_surface_points_on_flip_geodesic_path(
    int64_t startFaceIdx,
    int64_t endFaceIdx,
    const std::vector<double> &startFaceCoords,
    const std::vector<double> &endFaceCoords) {
    std::vector<GCS::SurfacePoint> path =
        find_flip_geodesic_path_GC(startFaceIdx, endFaceIdx, startFaceCoords, endFaceCoords);

    // convert the path to a vector of face index and barycentric coordinates
    return surface_points_to_face_idx_coords_list(path);
}

// Compute the geodesic path between two points, return a n*3 matrix
// where n is the number of points in the path
std::vector<GCS::SurfacePoint> GeodesicsManager::find_exact_geodesic_path_GC(int64_t startFaceIdx,
    int64_t endFaceIdx,
    const std::vector<double> &startFaceCoords,
    const std::vector<double> &endFaceCoords) {
    GCS::SurfacePoint start_converted, end_converted;
    auto early = preprocess_surface_points(
        startFaceIdx, endFaceIdx, startFaceCoords, endFaceCoords, start_converted, end_converted);
    if (early.has_value()) {
        return *early;
    }
    // now the two surface points are not the edge cases(sharing the same faces)

    // use the exact geodesic algorithm to find the path
    std::unique_ptr<GCS::GeodesicAlgorithmExact> mmp = std::make_unique<GCS::GeodesicAlgorithmExact>(*mesh, *geometry);

    double max_propagation_distance = GC::surface::GEODESIC_INF;
    mmp->propagate({start_converted}, max_propagation_distance, {end_converted});

    // Get the geodesic path from the end point to the start point
    double pathLength;
    std::vector<GCS::SurfacePoint> path = mmp->traceBack(end_converted, pathLength);
    
    // Reverse the path to return source -> target instead of target -> source
    std::reverse(path.begin(), path.end());
    return path;
};

// Compute the geodesic path between two points, return a n*3 matrix
// where n is the number of points in the path
GC::DenseMatrix<double> GeodesicsManager::find_exact_geodesic_path(int64_t startFaceIdx,
    int64_t endFaceIdx,
    const std::vector<double> &startFaceCoords,
    const std::vector<double> &endFaceCoords) {
    std::vector<GCS::SurfacePoint> path =
        find_exact_geodesic_path_GC(startFaceIdx, endFaceIdx, startFaceCoords, endFaceCoords);

    // convert the path to a matrix
    return surface_points_to_dense_matrix(path);
};

// Compute the geodesic path between two points, return a n*3 matrix
// where n is the number of points in the path
GC::DenseMatrix<double> GeodesicsManager::find_exact_geodesic_path(
    std::pair<int64_t, const std::vector<double> &> start_face_idx_coords_pair,
    std::pair<int64_t, const std::vector<double> &> end_face_idx_coords_pair) {
    // unpack the input
    int64_t endFaceIdx = end_face_idx_coords_pair.first;
    int64_t startFaceIdx = start_face_idx_coords_pair.first;
    std::vector<double> endFaceCoords = end_face_idx_coords_pair.second;
    std::vector<double> startFaceCoords = start_face_idx_coords_pair.second;

    return find_exact_geodesic_path(startFaceIdx, endFaceIdx, startFaceCoords, endFaceCoords);
};

// Compute the geodesic path between two points, return the path as a vector of surface points
std::vector<std::pair<int64_t, std::vector<double>>> GeodesicsManager::find_surface_points_on_exact_geodesic_path(
    int64_t startFaceIdx,
    int64_t endFaceIdx,
    const std::vector<double> &startFaceCoords,
    const std::vector<double> &endFaceCoords) {
    std::vector<GCS::SurfacePoint> path =
        find_exact_geodesic_path_GC(startFaceIdx, endFaceIdx, startFaceCoords, endFaceCoords);

    // convert the path to a vector of face index and barycentric coordinates
    return surface_points_to_face_idx_coords_list(path);
}

}  // namespace geomp