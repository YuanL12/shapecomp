#pragma once

#include "geomp/torus_comparison/Planar.h"
#include "geomp/mesh/Mesh.h"
#include <vector>
#include <unordered_map>
#include <queue>
#include <limits>
#include <memory>
#include <Eigen/Core>
#include "geometrycentral/surface/exact_geodesics.h"
#include "geometrycentral/surface/flip_geodesics.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geomp {

namespace GCS = geometrycentral::surface;
namespace GC = geometrycentral;

/**
 * Represents a vertex in the F3 tiled mesh.
 * Each vertex knows its tile offset and the base vertex index in the fundamental domain F.
 */
struct F3Vertex {
    int tileU;  // Tile offset in U direction: -1, 0, or 1
    int tileV;  // Tile offset in V direction: -1, 0, or 1
    int baseVertexIndex;  // Index of the vertex in the cutted mesh (fundamental domain F)
    int f3VertexIndex;    // Unique index in the F3 mesh
    
    F3Vertex() : tileU(0), tileV(0), baseVertexIndex(-1), f3VertexIndex(-1) {}
    F3Vertex(int u, int v, int baseIdx, int f3Idx) 
        : tileU(u), tileV(v), baseVertexIndex(baseIdx), f3VertexIndex(f3Idx) {}
    
    bool operator==(const F3Vertex& other) const {
        return f3VertexIndex == other.f3VertexIndex;
    }
    
    bool operator<(const F3Vertex& other) const {
        return f3VertexIndex < other.f3VertexIndex;
    }
};

/**
 * Represents an edge in the F3 tiled mesh.
 */
struct F3Edge {
    int v1;  // F3 vertex index
    int v2;  // F3 vertex index
    int baseEdgeIndex;  // Index of the base edge in the cutted mesh (or -1 if seam edge)
    
    F3Edge() : v1(-1), v2(-1), baseEdgeIndex(-1) {}
    F3Edge(int v1_, int v2_, int baseEdgeIdx) 
        : v1(v1_), v2(v2_), baseEdgeIndex(baseEdgeIdx) {}
};

/**
 * F3TileMesh: A 3×3 tiled mesh of the fundamental domain.
 * This creates 9 copies of the fundamental domain F, each translated by UV shifts in {-1,0,1}^2.
 */
class F3TileMesh {
public:
    // Build F3 from a PlanarLocator (which contains the fundamental domain F)
    F3TileMesh(const PlanarLocator& planarLocator);
    
    // Find the nearest vertex in F3 to a given UV point
    F3Vertex nearestVertex(const std::pair<double, double>& pointUV) const;
    
    // Get the F3 vertex at a given index
    const F3Vertex& getVertex(int f3VertexIndex) const;
    
    // Get neighbor F3 vertex indices for a given vertex (for Dijkstra)
    std::vector<int> getNeighborIndices(int f3VertexIndex) const;
    
    // Get the base vertex index in the original torus T for an F3 vertex
    int getBaseVertexIndexOnT(const F3Vertex& f3Vertex, const Mesh& cuttedMesh) const;
    
    // Get the number of vertices in F3
    size_t numVertices() const { return vertices.size(); }
    
private:
    std::vector<F3Vertex> vertices;
    std::vector<F3Edge> edges;
    std::unordered_map<int, std::vector<int>> vertexToNeighbors;  // F3 vertex index -> neighbor F3 vertex indices
    
    const Mesh* cuttedMeshPtr;  // Pointer to the fundamental domain mesh F
    const VertexData<std::pair<double, double>>* uvPositionsPtr;  // UV coordinates
    
    // Helper: Build edges for a single tile
    void buildTile(int tileU, int tileV);
    
    // Helper: Wire seam edges between adjacent tiles
    void wireSeamEdges();
    
    // Helper: Get tile offset for a UV coordinate
    std::pair<int, int> getTileOffset(const std::pair<double, double>& uv) const;
    
    // Helper: Convert tile UV to base domain UV [0,1)^2
    std::pair<double, double> tileUVToBaseUV(const std::pair<double, double>& tileUV, int tileU, int tileV) const;
};

// Forward declarations - implementations are in .cpp file which has geometry-central includes
// These functions are only used in code that links to geometry-central

/**
 * Compute geodesic path using the subdivided mesh from PlanarLocator.
 * The PlanarLocator already contains the subdivided mesh, so no additional mesh data is needed.
 * 
 * @param planarLocator: PlanarLocator containing the fundamental domain F and subdivided mesh
 * @param a: Start point in UV coordinates (can be in [-1,2]^2)
 * @param b: End point in UV coordinates (can be in [-1,2]^2)
 * @return: Nx3 matrix of 3D points representing the geodesic path
 */
Eigen::MatrixXd computeGeodesicPathFromPlanarLocator(
    const PlanarLocator& planarLocator,
    const std::pair<double, double>& a,
    const std::pair<double, double>& b);

/**
 * Debug version: Compute the input halfedge path (before shortening) as 3D coordinates.
 * Returns the vertex path along the halfedges as a 3D numpy array.
 * 
 * @param planarLocator: PlanarLocator containing the fundamental domain F and subdivided mesh
 * @param a: Start point in UV coordinates (can be in [-1,2]^2)
 * @param b: End point in UV coordinates (can be in [-1,2]^2)
 * @return: Nx3 matrix of 3D points representing the halfedge path vertices
 */
Eigen::MatrixXd computeDijkstraF3Path3D(
    const PlanarLocator& planarLocator,
    const std::pair<double, double>& a,
    const std::pair<double, double>& b);

/**
 * Debug version: Compute the dijkstraF3 path in 2D UV coordinates.
 * Returns the F3 vertex path projected to 2D UV space as a numpy array.
 * 
 * @param planarLocator: PlanarLocator containing the fundamental domain F and subdivided mesh
 * @param a: Start point in UV coordinates (can be in [-1,2]^2)
 * @param b: End point in UV coordinates (can be in [-1,2]^2)
 * @return: Nx2 matrix of UV coordinates representing the F3 vertex path
 */
Eigen::MatrixXd computeDijkstraF3Path2D(
    const PlanarLocator& planarLocator,
    const std::pair<double, double>& a,
    const std::pair<double, double>& b);

} // namespace geomp

