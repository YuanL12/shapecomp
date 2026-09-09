#include <cmath>
#include <iostream>
#include <iterator>

#include "extern/geometry_central/geodesic.hpp"

std::string to_string(GCS::SurfacePointType surface_point_type) {
    switch (surface_point_type) {
        case GCS::SurfacePointType::Vertex:
            return "vertex";
        case GCS::SurfacePointType::Edge:
            return "edge";
        case GCS::SurfacePointType::Face:
            return "face";
    }
    return "unknown";
}

int test_surface_points_conversion(std::string mesh_path) {
    geomp::GeodesicsManager geodesic_manager(mesh_path);

    int64_t face_idx = 2278;  // the start face index

    std::vector<std::pair<GC::Vector3, GCS::SurfacePointType>> face_coords_type_list = {
        {GC::Vector3{0.6, 0.0, 0.4}, GCS::SurfacePointType::Edge},                      // edge
        {GC::Vector3{0.5, 0.2, 0.3}, GCS::SurfacePointType::Face},                      // face
        {GC::Vector3{0.0, 0.0, 1.0}, GCS::SurfacePointType::Vertex},                    // vertex
        {GC::Vector3{0.0, 0.53846, 0.46153999999999995}, GCS::SurfacePointType::Edge},  // edge
    };

    for (const auto &face_coords_type : face_coords_type_list) {
        GC::Vector3 face_coords = face_coords_type.first;
        GCS::SurfacePointType target_surface_point_type = face_coords_type.second;
        try {
            GCS::ManifoldSurfaceMesh *mesh = geodesic_manager.get_mesh();
            GCS::VertexPositionGeometry *geometry = geodesic_manager.get_geometry();
            GCS::SurfacePoint surface_point = GCS::SurfacePoint(mesh->face(face_idx), face_coords);
            auto converted_surface_point = geodesic_manager.convert_face_point_to_specific(surface_point);

            // assert the converted surface points are the same as the original surface points in 3D
            auto surface_3d_position_original = surface_point.interpolate(geometry->vertexPositions);
            auto surface_3d_position_converted = converted_surface_point.interpolate(geometry->vertexPositions);
            GC::Vector3 diff = surface_3d_position_original - surface_3d_position_converted;
            if (diff.norm() > 1e-6) {
                throw std::runtime_error("diff is not zero, norm: " + std::to_string(diff.norm()));
            }
            if (converted_surface_point.type != target_surface_point_type) {
                std::string input_surface_point_type_name = to_string(surface_point.type);
                std::string target_surface_point_type_name = to_string(target_surface_point_type);
                throw std::runtime_error(
                    "converted surface point type(" + input_surface_point_type_name +
                    ") is not the target surface point type(" + target_surface_point_type_name + ")");
            }
        } catch (const std::exception &e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return 1;
        }
    }
    std::cout << "All tests passed" << std::endl;
    return 0;
}

int test_flip_geodesic_on_two_surface_points(std::string mesh_path) {
    geomp::GeodesicsManager geodesic_manager(mesh_path);
    int64_t start_face_idx = 2278;
    int64_t end_face_idx = 112;
    try {
        std::vector<double> start_face_coords = {0.0, 1.0, 0.0};                    // the start edge coordinates
        std::vector<double> end_face_coords = {0.0, 0.53846, 0.46153999999999995};  // the end face coordinates
        GC::DenseMatrix<double> path =
            geodesic_manager.find_flip_geodesic_path(start_face_idx, end_face_idx, start_face_coords, end_face_coords);
        std::cout << "path matrix: \n" << path << std::endl;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

int test_flip_geodesic_on_two_surface_points_2(std::string mesh_path) {
    try {
        geomp::GeodesicsManager geodesic_manager(mesh_path);
        int64_t start_face_idx = 6;
        int64_t end_face_idx = 50;
        std::vector<double> start_face_coords = {0.08691, 0.9130900000000001, 0.0};  
        std::vector<double> end_face_coords = {0.3, 0.4, 0.3};  

        GCS::SurfacePoint start_surface_point = GCS::SurfacePoint(geodesic_manager.get_mesh()->face(start_face_idx), GC::Vector3{start_face_coords[0], start_face_coords[1], start_face_coords[2]});
        GCS::SurfacePoint end_surface_point = GCS::SurfacePoint(geodesic_manager.get_mesh()->face(end_face_idx), GC::Vector3{end_face_coords[0], end_face_coords[1], end_face_coords[2]});
        GC::Vector3 start_vertex_position = start_surface_point.interpolate(geodesic_manager.get_geometry()->vertexPositions);
        GC::Vector3 end_vertex_position = end_surface_point.interpolate(geodesic_manager.get_geometry()->vertexPositions);
        std::vector<GCS::SurfacePoint> path =
            geodesic_manager.find_flip_geodesic_path_GC(start_face_idx, end_face_idx, start_face_coords, end_face_coords);
        
        std::cout << "start_vertex_position: " << start_vertex_position << std::endl;
        std::cout << "end_vertex_position: " << end_vertex_position << std::endl;
        std::cout << "path matrix: \n";
        for (const auto& point : path) {
            std::cout << point << " " << point.interpolate(geodesic_manager.get_geometry()->vertexPositions) << std::endl;
        }

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

int test_flip_geodesic_on_two_surface_points_3(std::string mesh_path) {
    try {
        geomp::GeodesicsManager geodesic_manager(mesh_path);
        int64_t start_face_idx = 6;
        int64_t end_face_idx = 113;
        std::vector<double> start_face_coords = {0.003999999999999954, 0.996, 0.0};  
        std::vector<double> end_face_coords = {0.08691, 0.9130900000000001, 0.0};  

        GCS::SurfacePoint start_surface_point = GCS::SurfacePoint(geodesic_manager.get_mesh()->face(start_face_idx), GC::Vector3{start_face_coords[0], start_face_coords[1], start_face_coords[2]});
        GCS::SurfacePoint end_surface_point = GCS::SurfacePoint(geodesic_manager.get_mesh()->face(end_face_idx), GC::Vector3{end_face_coords[0], end_face_coords[1], end_face_coords[2]});
        
        // Convert to specific types (Vertex/Edge/Face) to check if they're on the same edge
        GCS::SurfacePoint start_converted = geomp::GeodesicsManager::convert_face_point_to_specific(start_surface_point);
        GCS::SurfacePoint end_converted = geomp::GeodesicsManager::convert_face_point_to_specific(end_surface_point);
        
        std::cout << "Start converted type: " << to_string(start_converted.type) << std::endl;
        std::cout << "End converted type: " << to_string(end_converted.type) << std::endl;
        
        // Check if they're on the same edge
        bool on_same_edge = false;
        if (start_converted.type == GCS::SurfacePointType::Edge && end_converted.type == GCS::SurfacePointType::Edge) {
            on_same_edge = (start_converted.edge == end_converted.edge);
            std::cout << "On same edge: " << (on_same_edge ? "YES" : "NO") << std::endl;
            if (on_same_edge) {
                std::cout << "  Start tEdge: " << start_converted.tEdge << std::endl;
                std::cout << "  End tEdge: " << end_converted.tEdge << std::endl;
            }
        }
        
        GC::Vector3 start_vertex_position = start_surface_point.interpolate(geodesic_manager.get_geometry()->vertexPositions);
        GC::Vector3 end_vertex_position = end_surface_point.interpolate(geodesic_manager.get_geometry()->vertexPositions);
        std::vector<GCS::SurfacePoint> path =
            geodesic_manager.find_flip_geodesic_path_GC(start_face_idx, end_face_idx, start_face_coords, end_face_coords);
        
        std::cout << "start_vertex_position: " << start_vertex_position << std::endl;
        std::cout << "end_vertex_position: " << end_vertex_position << std::endl;
        std::cout << "path matrix: \n";
        for (const auto& point : path) {
            std::cout << point << " " << point.interpolate(geodesic_manager.get_geometry()->vertexPositions) << std::endl;
        }

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

int test_flip_geodesic_on_two_surface_points_4(std::string mesh_path) {
    try {
        geomp::GeodesicsManager geodesic_manager(mesh_path);
        int64_t start_face_idx = 52;
        int64_t end_face_idx = 53;
        std::vector<double> start_face_coords = {6.55112222065421e-17, 0.3, 0.7};  
        std::vector<double> end_face_coords = {0.93828, 0.00743, 0.054290000000000005};  

        GCS::SurfacePoint start_surface_point = GCS::SurfacePoint(geodesic_manager.get_mesh()->face(start_face_idx), GC::Vector3{start_face_coords[0], start_face_coords[1], start_face_coords[2]});
        GCS::SurfacePoint end_surface_point = GCS::SurfacePoint(geodesic_manager.get_mesh()->face(end_face_idx), GC::Vector3{end_face_coords[0], end_face_coords[1], end_face_coords[2]});
        
        // Convert to specific types (Vertex/Edge/Face) to check if they're on the same edge
        GCS::SurfacePoint start_converted = geomp::GeodesicsManager::convert_face_point_to_specific(start_surface_point);
        GCS::SurfacePoint end_converted = geomp::GeodesicsManager::convert_face_point_to_specific(end_surface_point);
        
        std::cout << "Start converted type: " << to_string(start_converted.type) << std::endl;
        std::cout << "End converted type: " << to_string(end_converted.type) << std::endl;
        
        // Check if they're on the same face
        bool on_same_face = false;
        if (start_converted.type == GCS::SurfacePointType::Edge && end_converted.type == GCS::SurfacePointType::Face) {
            // start point is on an edge, check its face and its twin face is the same as the end point
            auto start_half_edge = start_converted.edge.halfedge();
            on_same_face = (start_half_edge.face() == end_converted.face);
            on_same_face = on_same_face || (start_half_edge.twin().face() == end_converted.face);
            std::cout << "On same face: " << (on_same_face ? "YES" : "NO") << std::endl;
        }
        if (start_converted.type == GCS::SurfacePointType::Face && end_converted.type == GCS::SurfacePointType::Edge) {
            // end point is on an edge, check its face and its twin face is the same as the start point
            auto end_half_edge = end_converted.edge.halfedge();
            on_same_face = (end_half_edge.face() == start_converted.face);
            on_same_face = on_same_face || (end_half_edge.twin().face() == start_converted.face);
            std::cout << "On same face: " << (on_same_face ? "YES" : "NO") << std::endl;
        }
        
        GC::Vector3 start_vertex_position = start_surface_point.interpolate(geodesic_manager.get_geometry()->vertexPositions);
        GC::Vector3 end_vertex_position = end_surface_point.interpolate(geodesic_manager.get_geometry()->vertexPositions);
        std::vector<GCS::SurfacePoint> path =
            geodesic_manager.find_flip_geodesic_path_GC(start_face_idx, end_face_idx, start_face_coords, end_face_coords);
        
        std::cout << "start_vertex_position: " << start_vertex_position << std::endl;
        std::cout << "end_vertex_position: " << end_vertex_position << std::endl;
        std::cout << "path matrix: \n";
        for (const auto& point : path) {
            std::cout << point << " " << point.interpolate(geodesic_manager.get_geometry()->vertexPositions) << std::endl;
        }

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

int main() {
    // Prefer argv mesh path; default to the tracked torus fixture (run from repo root).
    std::string b2_path = "test/input/torus/torus.obj";
    // test_surface_points_conversion(mesh_path);
    // std::cout << "--------------------------------" << std::endl;
    // test_flip_geodesic_on_two_surface_points_2(b2_path);
    // std::cout << "--------------------------------" << std::endl;
    // std::cout << "test_flip_geodesic_on two surface points on the same edge" << std::endl;
    // test_flip_geodesic_on_two_surface_points_3(b2_path);
    // std::cout << "--------------------------------" << std::endl;
    std::cout << "test_flip_geodesic_on two surface points: one on the edge, one on the face" << std::endl;
    test_flip_geodesic_on_two_surface_points_4(b2_path);
    std::cout << "--------------------------------" << std::endl;
    return 0;
}
