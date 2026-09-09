#include "geomp/torus_comparison/GeodesicPath.h"
#include "geomp/mesh/GeometryUtils.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/flip_geodesics.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <tuple>
#include <Eigen/Dense>

namespace GC = geometrycentral;

namespace geomp {

// Helper: Get tile offset for a UV coordinate
std::pair<int, int> F3TileMesh::getTileOffset(const std::pair<double, double>& uv) const {
    // Convert UV to tile coordinates
    // UV space is [-1, 2]^2, tiles are at shifts {-1, 0, 1}^2
    int tileU = static_cast<int>(std::floor(uv.first));
    int tileV = static_cast<int>(std::floor(uv.second));
    
    // Clamp to valid tile range
    tileU = std::max(-1, std::min(1, tileU));
    tileV = std::max(-1, std::min(1, tileV));
    
    return {tileU, tileV};
}

// Helper: Convert tile UV to base domain UV [0,1)^2
std::pair<double, double> F3TileMesh::tileUVToBaseUV(
    const std::pair<double, double>& tileUV, int tileU, int tileV) const {
    return {
        tileUV.first - tileU,
        tileUV.second - tileV
    };
}

F3TileMesh::F3TileMesh(const PlanarLocator& planarLocator) {
    cuttedMeshPtr = planarLocator.get_cutted_mesh_ptr().get();
    uvPositionsPtr = planarLocator.get_uv_positions_ptr().get();
    
    if (!cuttedMeshPtr || !uvPositionsPtr) {
        throw std::runtime_error("F3TileMesh: PlanarLocator must have valid cutted mesh and UV positions");
    }
    
    // const Mesh& F = *cuttedMeshPtr;
    // const auto& uvPositions = *uvPositionsPtr;
    
    // Build 9 tiles
    for (int tileU = -1; tileU <= 1; ++tileU) {
        for (int tileV = -1; tileV <= 1; ++tileV) {
            buildTile(tileU, tileV);
        }
    }
    
    // Wire seam edges between adjacent tiles
    wireSeamEdges();
    
    // Build adjacency list for Dijkstra (avoid duplicates)
    vertexToNeighbors.clear();
    for (const auto& edge : edges) {
        // Check if this neighbor is already added to avoid duplicates
        auto& neighbors1 = vertexToNeighbors[edge.v1];
        if (std::find(neighbors1.begin(), neighbors1.end(), edge.v2) == neighbors1.end()) {
            neighbors1.push_back(edge.v2);
        }
        
        auto& neighbors2 = vertexToNeighbors[edge.v2];
        if (std::find(neighbors2.begin(), neighbors2.end(), edge.v1) == neighbors2.end()) {
            neighbors2.push_back(edge.v1);
        }
    }
}

void F3TileMesh::buildTile(int tileU, int tileV) {
    const Mesh& cuttedMesh = *cuttedMeshPtr;
    const auto& uvPositions = *uvPositionsPtr;
    
    // Create vertices for this tile
    std::unordered_map<int, int> baseVertexToF3Index;  // Maps base vertex index to F3 vertex index in this tile
    
    for (VertexCIter v = cuttedMesh.vertices.begin(); v != cuttedMesh.vertices.end(); ++v) {
        int baseVertexIndex = v->index;
        int f3VertexIndex = static_cast<int>(vertices.size());
        
        F3Vertex f3v(tileU, tileV, baseVertexIndex, f3VertexIndex);
        vertices.push_back(f3v);
        baseVertexToF3Index[baseVertexIndex] = f3VertexIndex;
    }
    
    // Create edges for this tile (interior edges only - seam edges will be wired separately)
    for (EdgeCIter e = cuttedMesh.edges.begin(); e != cuttedMesh.edges.end(); ++e) {
        auto [v1, v2] = e->twoEndpoints();
        int baseV1Idx = v1->index;
        int baseV2Idx = v2->index;
        
        // Find F3 vertex indices for this tile
        auto it1 = baseVertexToF3Index.find(baseV1Idx);
        auto it2 = baseVertexToF3Index.find(baseV2Idx);
        
        if (it1 != baseVertexToF3Index.end() && it2 != baseVertexToF3Index.end()) {
            int f3v1 = it1->second;
            int f3v2 = it2->second;
            
            // Add all interior edges (seam edges will be wired separately in wireSeamEdges)
            edges.emplace_back(f3v1, f3v2, e->index);
        }
    }
}

void F3TileMesh::wireSeamEdges() {
    const Mesh& cuttedMesh = *cuttedMeshPtr;
    const auto& uvPositions = *uvPositionsPtr;
    
    // Build maps for efficient lookup
    // Map: (tileU, tileV, baseVertexIndex) -> F3 vertex index
    std::map<std::tuple<int, int, int>, int> tileVertexMap;
    // Map: (tileU, tileV, referenceIndex) -> F3 vertex indices (for boundary vertices)
    std::map<std::tuple<int, int, int>, std::vector<int>> referenceMap;
    
    for (size_t i = 0; i < vertices.size(); ++i) {
        const auto& f3v = vertices[i];
        tileVertexMap[{f3v.tileU, f3v.tileV, f3v.baseVertexIndex}] = static_cast<int>(i);
        
        // Also build referenceIndex map for boundary vertices
        VertexCIter baseV = cuttedMesh.vertices.begin() + f3v.baseVertexIndex;
        if (baseV->referenceIndex != -1) {
            referenceMap[{f3v.tileU, f3v.tileV, baseV->referenceIndex}].push_back(static_cast<int>(i));
        }
    }
    
    // Use a set to track edges we've already added to avoid duplicates
    std::set<std::pair<int, int>> addedEdges;  // (min(v1, v2), max(v1, v2))
    
    const double EPS = 1e-10;  // Small epsilon for floating point comparison
    
    // New approach: Iterate over tiles and connect boundary vertices directly
    // For each tile, check its boundary vertices and connect them to adjacent tiles
    for (int tileU = -1; tileU <= 1; ++tileU) {
        for (int tileV = -1; tileV <= 1; ++tileV) {
            // For each vertex in this tile, check if it's on a boundary
            for (const auto& f3v : vertices) {
                if (f3v.tileU != tileU || f3v.tileV != tileV) {
                    continue;  // Not in this tile
                }
                
                VertexCIter baseV = cuttedMesh.vertices.begin() + f3v.baseVertexIndex;
                const auto& uv = uvPositions[baseV];
                int baseVertexIdx = baseV->index;
                int refIdx = baseV->referenceIndex;
                
                // Check which boundaries this vertex is on
                bool onLeftBoundary = (std::abs(uv.first - 0.0) < EPS);
                bool onRightBoundary = (std::abs(uv.first - 1.0) < EPS);
                bool onBottomBoundary = (std::abs(uv.second - 0.0) < EPS);
                bool onTopBoundary = (std::abs(uv.second - 1.0) < EPS);
                
                // Connect left boundary (u == 0) to right boundary (u == 1) in left neighbor tile
                if (onLeftBoundary && tileU > -1) {
                    int neighborTileU = tileU - 1;
                    
                    // Find corresponding vertex in the left neighbor tile
                    std::vector<int> targetVertices;
                    if (refIdx != -1) {
                        // Boundary vertex - match by referenceIndex
                        auto refIt = referenceMap.find({neighborTileU, tileV, refIdx});
                        if (refIt != referenceMap.end()) {
                            targetVertices = refIt->second;
                        }
                    } else {
                        // Non-boundary vertex on boundary - this shouldn't happen, but handle it
                        auto it = tileVertexMap.find({neighborTileU, tileV, baseVertexIdx});
                        if (it != tileVertexMap.end()) {
                            targetVertices.push_back(it->second);
                        }
                    }
                    
                    // Connect to matching vertices
                    for (int targetIdx : targetVertices) {
                        if (f3v.f3VertexIndex == targetIdx) continue;  // Skip self-loops
                        std::pair<int, int> edgeKey = {std::min(f3v.f3VertexIndex, targetIdx), 
                                                       std::max(f3v.f3VertexIndex, targetIdx)};
                        if (addedEdges.find(edgeKey) == addedEdges.end()) {
                            edges.emplace_back(f3v.f3VertexIndex, targetIdx, -1);
                            addedEdges.insert(edgeKey);
                        }
                    }
                }
                
                // Connect right boundary (u == 1) to left boundary (u == 0) in right neighbor tile
                if (onRightBoundary && tileU < 1) {
                    int neighborTileU = tileU + 1;
                    
                    std::vector<int> targetVertices;
                    if (refIdx != -1) {
                        auto refIt = referenceMap.find({neighborTileU, tileV, refIdx});
                        if (refIt != referenceMap.end()) {
                            targetVertices = refIt->second;
                        }
                    } else {
                        auto it = tileVertexMap.find({neighborTileU, tileV, baseVertexIdx});
                        if (it != tileVertexMap.end()) {
                            targetVertices.push_back(it->second);
                        }
                    }
                    
                    for (int targetIdx : targetVertices) {
                        if (f3v.f3VertexIndex == targetIdx) continue;
                        std::pair<int, int> edgeKey = {std::min(f3v.f3VertexIndex, targetIdx), 
                                                       std::max(f3v.f3VertexIndex, targetIdx)};
                        if (addedEdges.find(edgeKey) == addedEdges.end()) {
                            edges.emplace_back(f3v.f3VertexIndex, targetIdx, -1);
                            addedEdges.insert(edgeKey);
                        }
                    }
                }
                
                // Connect bottom boundary (v == 0) to top boundary (v == 1) in bottom neighbor tile
                if (onBottomBoundary && tileV > -1) {
                    int neighborTileV = tileV - 1;
                    
                    std::vector<int> targetVertices;
                    if (refIdx != -1) {
                        auto refIt = referenceMap.find({tileU, neighborTileV, refIdx});
                        if (refIt != referenceMap.end()) {
                            targetVertices = refIt->second;
                        }
                    } else {
                        auto it = tileVertexMap.find({tileU, neighborTileV, baseVertexIdx});
                        if (it != tileVertexMap.end()) {
                            targetVertices.push_back(it->second);
                        }
                    }
                    
                    for (int targetIdx : targetVertices) {
                        if (f3v.f3VertexIndex == targetIdx) continue;
                        std::pair<int, int> edgeKey = {std::min(f3v.f3VertexIndex, targetIdx), 
                                                       std::max(f3v.f3VertexIndex, targetIdx)};
                        if (addedEdges.find(edgeKey) == addedEdges.end()) {
                            edges.emplace_back(f3v.f3VertexIndex, targetIdx, -1);
                            addedEdges.insert(edgeKey);
                        }
                    }
                }
                
                // Connect top boundary (v == 1) to bottom boundary (v == 0) in top neighbor tile
                if (onTopBoundary && tileV < 1) {
                    int neighborTileV = tileV + 1;
                    
                    std::vector<int> targetVertices;
                    if (refIdx != -1) {
                        auto refIt = referenceMap.find({tileU, neighborTileV, refIdx});
                        if (refIt != referenceMap.end()) {
                            targetVertices = refIt->second;
                        }
                    } else {
                        auto it = tileVertexMap.find({tileU, neighborTileV, baseVertexIdx});
                        if (it != tileVertexMap.end()) {
                            targetVertices.push_back(it->second);
                        }
                    }
                    
                    for (int targetIdx : targetVertices) {
                        if (f3v.f3VertexIndex == targetIdx) continue;
                        std::pair<int, int> edgeKey = {std::min(f3v.f3VertexIndex, targetIdx), 
                                                       std::max(f3v.f3VertexIndex, targetIdx)};
                        if (addedEdges.find(edgeKey) == addedEdges.end()) {
                            edges.emplace_back(f3v.f3VertexIndex, targetIdx, -1);
                            addedEdges.insert(edgeKey);
                        }
                    }
                }
            }
        }
    }
}

F3Vertex F3TileMesh::nearestVertex(const std::pair<double, double>& pointUV) const {
    const auto& uvPositions = *uvPositionsPtr;
    
    // First, determine which tile(s) this point belongs to
    auto [tileU, tileV] = getTileOffset(pointUV);
    
    // Convert pointUV to base domain UV [0,1]^2
    auto basePointUV = tileUVToBaseUV(pointUV, tileU, tileV);
    
    // Search in the primary tile and adjacent tiles (for points near boundaries)
    double minDist = std::numeric_limits<double>::max();
    F3Vertex nearest;
    
    // Search in relevant tiles: primary tile and possibly adjacent ones
    for (int searchTileU = std::max(-1, tileU - 1); searchTileU <= std::min(1, tileU + 1); ++searchTileU) {
        for (int searchTileV = std::max(-1, tileV - 1); searchTileV <= std::min(1, tileV + 1); ++searchTileV) {
            // Search all vertices in this tile
            for (const auto& f3v : vertices) {
                if (f3v.tileU != searchTileU || f3v.tileV != searchTileV) {
                    continue;
                }
                
                // Get the UV position of the base vertex
                VertexCIter baseV = cuttedMeshPtr->vertices.begin() + f3v.baseVertexIndex;
                const auto& baseUV = uvPositions[baseV];
                
                // Apply tile offset to get position in tiled space
                double tileUVx = baseUV.first + f3v.tileU;
                double tileUVy = baseUV.second + f3v.tileV;
                
                // Compute distance
                double dx = pointUV.first - tileUVx;
                double dy = pointUV.second - tileUVy;
                double dist = std::sqrt(dx * dx + dy * dy);
                
                if (dist < minDist) {
                    minDist = dist;
                    nearest = f3v;
                }
            }
        }
    }
    
    return nearest;
}

const F3Vertex& F3TileMesh::getVertex(int f3VertexIndex) const {
    return vertices[f3VertexIndex];
}

std::vector<int> F3TileMesh::getNeighborIndices(int f3VertexIndex) const {
    auto it = vertexToNeighbors.find(f3VertexIndex);
    if (it != vertexToNeighbors.end()) {
        return it->second;
    }
    return {};
}

// Helper function to get neighbors with weights (not a class method to avoid header dependency)
// Takes cuttedMesh as parameter to avoid accessing private members
// Uses 2D UV distances as edge weights instead of 3D distances
std::vector<std::pair<F3Vertex, double>> getNeighborsF3WithMesh(
    const F3TileMesh& f3,
    const F3Vertex& vertex,
    const Mesh& cuttedMesh,
    const VertexData<std::pair<double, double>>& uvPositions) {
    
    std::vector<std::pair<F3Vertex, double>> neighbors;
    
    // Get neighbor F3 vertex indices
    std::vector<int> neighborIndices = f3.getNeighborIndices(vertex.f3VertexIndex);
    if (neighborIndices.empty()) {
        return neighbors;
    }
    
    // Get the UV position and referenceIndex of the current vertex
    VertexCIter baseV = cuttedMesh.vertices.begin() + vertex.baseVertexIndex;
    const auto& baseUV = uvPositions[baseV];
    double vertexUVx = baseUV.first + vertex.tileU;
    double vertexUVy = baseUV.second + vertex.tileV;
    int vertexRefIdx = baseV->referenceIndex;
    
    // For each neighbor F3 vertex, compute the 2D UV edge weight
    for (int neighborF3Idx : neighborIndices) {
        const F3Vertex& neighbor = f3.getVertex(neighborF3Idx);
        
        // Get the UV position and referenceIndex of the neighbor vertex
        VertexCIter neighborBaseV = cuttedMesh.vertices.begin() + neighbor.baseVertexIndex;
        const auto& neighborBaseUV = uvPositions[neighborBaseV];
        int neighborRefIdx = neighborBaseV->referenceIndex;
        
        // Check if these vertices represent the same logical point (same referenceIndex)
        // If so, they should have weight 0 (they're the same point in periodic space)
        double weight;
        if (vertexRefIdx != -1 && neighborRefIdx != -1 && vertexRefIdx == neighborRefIdx) {
            // Both are boundary vertices with the same referenceIndex - same logical point
            weight = 0.0;
        } else {
            // Compute Euclidean distance in 2D UV space
            double neighborUVx = neighborBaseUV.first + neighbor.tileU;
            double neighborUVy = neighborBaseUV.second + neighbor.tileV;
            
            double dx = neighborUVx - vertexUVx;
            double dy = neighborUVy - vertexUVy;
            weight = std::sqrt(dx * dx + dy * dy);
        }
        
        neighbors.emplace_back(neighbor, weight);
    }
    
    return neighbors;
}

int F3TileMesh::getBaseVertexIndexOnT(const F3Vertex& f3Vertex, const Mesh& cuttedMesh) const {
    // The base vertex in the cutted mesh has an index
    // We need to map this to the original torus T vertex index
    // Strategy: If referenceIndex is set, it points to the original vertex in the subdivided mesh.
    // For the cut mesh, vertices with referenceIndex != -1 are duplicates on seams.
    // We always use referenceIndex if available, otherwise use the index itself.
    // Note: This assumes the geometry-central mesh T corresponds to the subdivided mesh (before cutting).
    
    VertexCIter baseV = cuttedMesh.vertices.begin() + f3Vertex.baseVertexIndex;
    if (baseV->referenceIndex != -1) {
        // This is a duplicate vertex on a seam - use the original vertex index
        return baseV->referenceIndex;
    } else {
        // This is a unique vertex - use its index directly
        // Note: This assumes the cut mesh indices for non-seam vertices match the subdivided mesh
        return baseV->index;
    }
}

std::vector<F3Vertex> dijkstraF3(
    const F3TileMesh& f3,
    const F3Vertex& start,
    const F3Vertex& end,
    const Mesh& cuttedMesh,
    const VertexData<std::pair<double, double>>& uvPositions) {
    
    // Dijkstra's algorithm using 2D UV distances as edge weights
    std::unordered_map<int, double> dist;  // F3 vertex index -> distance
    std::unordered_map<int, int> prev;      // F3 vertex index -> previous F3 vertex index
    std::priority_queue<std::pair<double, int>, 
                        std::vector<std::pair<double, int>>,
                        std::greater<std::pair<double, int>>> pq;  // (distance, vertex index)
    
    // Initialize distances
    for (size_t i = 0; i < f3.numVertices(); ++i) {
        dist[static_cast<int>(i)] = std::numeric_limits<double>::max();
    }
    
    dist[start.f3VertexIndex] = 0.0;
    pq.push({0.0, start.f3VertexIndex});
    
    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();
        
        if (d > dist[u]) {
            continue;  // Skip if we've found a better path
        }
        
        if (u == end.f3VertexIndex) {
            break;  // Reached destination
        }
        
        // Get neighbors with 2D UV weights
        const F3Vertex& uVertex = f3.getVertex(u);
        auto neighbors = getNeighborsF3WithMesh(f3, uVertex, cuttedMesh, uvPositions);
        
        for (const auto& [vVertex, weight] : neighbors) {
            int v = vVertex.f3VertexIndex;
            double alt = dist[u] + weight;
            
            if (alt < dist[v]) {
                dist[v] = alt;
                prev[v] = u;
                pq.push({alt, v});
            }
        }
    }
    
    // Reconstruct path
    std::vector<F3Vertex> path;
    if (dist[end.f3VertexIndex] < std::numeric_limits<double>::max()) {
        // Path exists
        int current = end.f3VertexIndex;
        while (current != start.f3VertexIndex) {
            path.push_back(f3.getVertex(current));
            auto it = prev.find(current);
            if (it == prev.end()) {
                break;  // No path found
            }
            current = it->second;
        }
        path.push_back(start);
        std::reverse(path.begin(), path.end());
    }
    
    return path;
}

std::vector<GCS::Halfedge> verticesToHalfedgesOnT(
    GCS::ManifoldSurfaceMesh& meshT,
    const std::vector<int>& vertexPath) {
    
    std::vector<GCS::Halfedge> halfedgePath;
    
    for (size_t i = 0; i + 1 < vertexPath.size(); ++i) {
        int uIdx = vertexPath[i];
        int vIdx = vertexPath[i + 1];
        
        if (uIdx < 0 || uIdx >= static_cast<int>(meshT.nVertices()) ||
            vIdx < 0 || vIdx >= static_cast<int>(meshT.nVertices())) {
            continue;
        }
        
        GCS::Vertex u = meshT.vertex(uIdx);
        GCS::Vertex v = meshT.vertex(vIdx);
        
        // Find the outgoing halfedge from u whose tip is v
        GCS::Halfedge foundHe;
        bool found = false;
        
        for (GCS::Halfedge he : u.outgoingHalfedges()) {
            if (he.tipVertex() == v) {
                foundHe = he;
                found = true;
                break;
            }
        }
        
        if (found) {
            halfedgePath.push_back(foundHe);
        } else {
            // Vertices are not adjacent - this shouldn't happen if vertexPath is valid
            // Skip this edge
            continue;
        }
    }
    
    return halfedgePath;
}

std::vector<GC::Vector3> computeGeodesicPath(
    const PlanarLocator& planarLocator,
    std::unique_ptr<GCS::ManifoldSurfaceMesh>& meshT,
    std::unique_ptr<GCS::VertexPositionGeometry>& geometryT,
    const std::pair<double, double>& a,
    const std::pair<double, double>& b) {
    
    // Step 1: Build F3
    F3TileMesh f3(planarLocator);
    
    // Step 2: Find nearest vertices in F3
    F3Vertex startVertex = f3.nearestVertex(a);
    F3Vertex endVertex = f3.nearestVertex(b);
    
    // Step 3: Run Dijkstra on F3
    const Mesh* cuttedMesh = planarLocator.get_cutted_mesh_ptr().get();
    if (!cuttedMesh) {
        throw std::runtime_error("PlanarLocator has no cutted mesh");
    }
    const auto& uvPositions = *planarLocator.get_uv_positions_ptr();
    std::vector<F3Vertex> f3Path = dijkstraF3(f3, startVertex, endVertex, *cuttedMesh, uvPositions);
    
    if (f3Path.empty()) {
        // Provide more diagnostic information
        std::string errorMsg = "No path found in F3. ";
        errorMsg += "Start vertex: F3 idx=" + std::to_string(startVertex.f3VertexIndex) + 
                    ", tile=(" + std::to_string(startVertex.tileU) + "," + std::to_string(startVertex.tileV) + 
                    "), baseIdx=" + std::to_string(startVertex.baseVertexIndex) + ". ";
        errorMsg += "End vertex: F3 idx=" + std::to_string(endVertex.f3VertexIndex) + 
                    ", tile=(" + std::to_string(endVertex.tileU) + "," + std::to_string(endVertex.tileV) + 
                    "), baseIdx=" + std::to_string(endVertex.baseVertexIndex) + ". ";
        errorMsg += "F3 has " + std::to_string(f3.numVertices()) + " vertices. ";
        
        // Check if vertices have neighbors
        std::vector<int> startNeighbors = f3.getNeighborIndices(startVertex.f3VertexIndex);
        std::vector<int> endNeighbors = f3.getNeighborIndices(endVertex.f3VertexIndex);
        errorMsg += "Start has " + std::to_string(startNeighbors.size()) + " neighbors. ";
        errorMsg += "End has " + std::to_string(endNeighbors.size()) + " neighbors.";
        
        throw std::runtime_error(errorMsg);
    }
    
    // Step 4: Map F3 vertices to base vertices on T
    std::vector<int> baseVertexPath;
    
    for (const auto& f3v : f3Path) {
        int baseVertexIdxOnT = f3.getBaseVertexIndexOnT(f3v, *cuttedMesh);
        baseVertexPath.push_back(baseVertexIdxOnT);
    }
    
    // Step 5: Convert to halfedge path
    std::vector<GCS::Halfedge> halfedgePath = verticesToHalfedgesOnT(*meshT, baseVertexPath);
    
    if (halfedgePath.empty()) {
        throw std::runtime_error("Failed to convert vertex path to halfedge path");
    }
    
    // Step 6: Build FlipEdgeNetwork and shorten
    std::vector<std::vector<GCS::Halfedge>> paths = {halfedgePath};
    GCS::VertexData<bool> pinnedVertices(*meshT, false);
    
    // Pin start and end vertices
    if (!baseVertexPath.empty()) {
        pinnedVertices[meshT->vertex(baseVertexPath[0])] = true;
        pinnedVertices[meshT->vertex(baseVertexPath.back())] = true;
    }
    
    // Create FlipEdgeNetwork
    std::unique_ptr<GCS::FlipEdgeNetwork> edgeNetwork = 
        std::make_unique<GCS::FlipEdgeNetwork>(*meshT, *geometryT, paths, pinnedVertices);
    
    // Shorten the path
    edgeNetwork->iterativeShorten();
    
    // Step 7: Extract the geodesic polyline
    edgeNetwork->posGeom = geometryT.get();
    std::vector<std::vector<GC::Vector3>> polylinePaths = edgeNetwork->getPathPolyline3D();
    
    if (polylinePaths.empty() || polylinePaths[0].empty()) {
        throw std::runtime_error("Failed to extract geodesic polyline");
    }
    
    return polylinePaths[0];
}

Eigen::MatrixXd computeGeodesicPathFromPlanarLocator(
    const PlanarLocator& planarLocator,
    const std::pair<double, double>& a,
    const std::pair<double, double>& b) {
    
    // Get the subdivided mesh from PlanarLocator (which contains the 3D positions)
    const Mesh* subdividedMesh = planarLocator.get_subdivided_mesh_ptr().get();
    if (!subdividedMesh) {
        throw std::runtime_error("PlanarLocator has no subdivided mesh");
    }
    
    // Use the subdivided mesh vertices and faces instead of the input
    Eigen::MatrixXd subdividedVertices = subdividedMesh->getVertices();
    Eigen::MatrixXi subdividedFaces = subdividedMesh->getFaces();
    
    // Convert Eigen matrices to geometry-central format
    GC::DenseMatrix<double> vertsGC(subdividedVertices.rows(), 3);
    GC::DenseMatrix<int64_t> facesGC(subdividedFaces.rows(), 3);
    
    for (int i = 0; i < subdividedVertices.rows(); ++i) {
        for (int j = 0; j < 3; ++j) {
            vertsGC(i, j) = subdividedVertices(i, j);
        }
    }
    
    for (int i = 0; i < subdividedFaces.rows(); ++i) {
        for (int j = 0; j < 3; ++j) {
            facesGC(i, j) = static_cast<int64_t>(subdividedFaces(i, j));
        }
    }
    
    // Create geometry-central mesh
    std::unique_ptr<GCS::ManifoldSurfaceMesh> meshT = 
        std::make_unique<GCS::ManifoldSurfaceMesh>(facesGC);
    std::unique_ptr<GCS::VertexPositionGeometry> geometryT = 
        std::make_unique<GCS::VertexPositionGeometry>(*meshT);
    
    for (size_t i = 0; i < meshT->nVertices(); ++i) {
        for (int j = 0; j < 3; ++j) {
            geometryT->inputVertexPositions[i][j] = vertsGC(i, j);
        }
    }
    
    // Compute geodesic path
    std::vector<GC::Vector3> path = computeGeodesicPath(
        planarLocator, meshT, geometryT, a, b);
    
    // Convert to Eigen matrix
    Eigen::MatrixXd pathMatrix(path.size(), 3);
    for (size_t i = 0; i < path.size(); ++i) {
        pathMatrix(i, 0) = path[i].x;
        pathMatrix(i, 1) = path[i].y;
        pathMatrix(i, 2) = path[i].z;
    }
    
    return pathMatrix;
}

Eigen::MatrixXd computeDijkstraF3Path3D(
    const PlanarLocator& planarLocator,
    const std::pair<double, double>& a,
    const std::pair<double, double>& b) {
    
    // Get the subdivided mesh from PlanarLocator (which contains the 3D positions)
    const Mesh* subdividedMesh = planarLocator.get_subdivided_mesh_ptr().get();
    if (!subdividedMesh) {
        throw std::runtime_error("PlanarLocator has no subdivided mesh");
    }
    
    // Use the subdivided mesh vertices and faces
    Eigen::MatrixXd subdividedVertices = subdividedMesh->getVertices();
    Eigen::MatrixXi subdividedFaces = subdividedMesh->getFaces();
    
    // Convert Eigen matrices to geometry-central format
    GC::DenseMatrix<double> vertsGC(subdividedVertices.rows(), 3);
    GC::DenseMatrix<int64_t> facesGC(subdividedFaces.rows(), 3);
    
    for (int i = 0; i < subdividedVertices.rows(); ++i) {
        for (int j = 0; j < 3; ++j) {
            vertsGC(i, j) = subdividedVertices(i, j);
        }
    }
    
    for (int i = 0; i < subdividedFaces.rows(); ++i) {
        for (int j = 0; j < 3; ++j) {
            facesGC(i, j) = static_cast<int64_t>(subdividedFaces(i, j));
        }
    }
    
    // Create geometry-central mesh
    std::unique_ptr<GCS::ManifoldSurfaceMesh> meshT = 
        std::make_unique<GCS::ManifoldSurfaceMesh>(facesGC);
    std::unique_ptr<GCS::VertexPositionGeometry> geometryT = 
        std::make_unique<GCS::VertexPositionGeometry>(*meshT);
    
    for (size_t i = 0; i < meshT->nVertices(); ++i) {
        for (int j = 0; j < 3; ++j) {
            geometryT->inputVertexPositions[i][j] = vertsGC(i, j);
        }
    }
    
    // Step 1: Build F3
    F3TileMesh f3(planarLocator);
    
    // Step 2: Find nearest vertices in F3
    F3Vertex startVertex = f3.nearestVertex(a);
    F3Vertex endVertex = f3.nearestVertex(b);
    
    // Step 3: Run Dijkstra on F3
    const Mesh* cuttedMesh = planarLocator.get_cutted_mesh_ptr().get();
    if (!cuttedMesh) {
        throw std::runtime_error("PlanarLocator has no cutted mesh");
    }
    const auto& uvPositions = *planarLocator.get_uv_positions_ptr();
    std::vector<F3Vertex> f3Path = dijkstraF3(f3, startVertex, endVertex, *cuttedMesh, uvPositions);
    
    if (f3Path.empty()) {
        std::string errorMsg = "No path found in F3. ";
        errorMsg += "Start vertex: F3 idx=" + std::to_string(startVertex.f3VertexIndex) + 
                    ", tile=(" + std::to_string(startVertex.tileU) + "," + std::to_string(startVertex.tileV) + 
                    "), baseIdx=" + std::to_string(startVertex.baseVertexIndex) + ". ";
        errorMsg += "End vertex: F3 idx=" + std::to_string(endVertex.f3VertexIndex) + 
                    ", tile=(" + std::to_string(endVertex.tileU) + "," + std::to_string(endVertex.tileV) + 
                    "), baseIdx=" + std::to_string(endVertex.baseVertexIndex) + ". ";
        errorMsg += "F3 has " + std::to_string(f3.numVertices()) + " vertices. ";
        
        std::vector<int> startNeighbors = f3.getNeighborIndices(startVertex.f3VertexIndex);
        std::vector<int> endNeighbors = f3.getNeighborIndices(endVertex.f3VertexIndex);
        errorMsg += "Start has " + std::to_string(startNeighbors.size()) + " neighbors. ";
        errorMsg += "End has " + std::to_string(endNeighbors.size()) + " neighbors.";
        
        throw std::runtime_error(errorMsg);
    }
    
    // Step 4: Map F3 vertices to base vertices on T
    std::vector<int> baseVertexPath;
    
    for (const auto& f3v : f3Path) {
        int baseVertexIdxOnT = f3.getBaseVertexIndexOnT(f3v, *cuttedMesh);
        baseVertexPath.push_back(baseVertexIdxOnT);
    }
    
    // Step 5: Convert to halfedge path
    std::vector<GCS::Halfedge> halfedgePath = verticesToHalfedgesOnT(*meshT, baseVertexPath);
    
    if (halfedgePath.empty()) {
        throw std::runtime_error("Failed to convert vertex path to halfedge path");
    }
    
    // Step 6: Convert halfedge path to 3D coordinates
    // Extract vertices along the halfedge path
    std::vector<GC::Vector3> path3D;
    
    if (!halfedgePath.empty()) {
        // Add the tail vertex of the first halfedge (start point)
        GCS::Vertex startV = halfedgePath[0].tailVertex();
        path3D.push_back(geometryT->vertexPositions[startV]);
        
        // Add the tip vertex of each halfedge (end points)
        for (const auto& he : halfedgePath) {
            GCS::Vertex tipV = he.tipVertex();
            path3D.push_back(geometryT->vertexPositions[tipV]);
        }
    }
    
    // Convert to Eigen matrix
    Eigen::MatrixXd pathMatrix(path3D.size(), 3);
    for (size_t i = 0; i < path3D.size(); ++i) {
        pathMatrix(i, 0) = path3D[i].x;
        pathMatrix(i, 1) = path3D[i].y;
        pathMatrix(i, 2) = path3D[i].z;
    }
    
    return pathMatrix;
}

Eigen::MatrixXd computeDijkstraF3Path2D(
    const PlanarLocator& planarLocator,
    const std::pair<double, double>& a,
    const std::pair<double, double>& b) {
    
    // Get the subdivided mesh from PlanarLocator (which contains the 3D positions)
    const Mesh* subdividedMesh = planarLocator.get_subdivided_mesh_ptr().get();
    if (!subdividedMesh) {
        throw std::runtime_error("PlanarLocator has no subdivided mesh");
    }
    
    // Use the subdivided mesh vertices and faces
    Eigen::MatrixXd subdividedVertices = subdividedMesh->getVertices();
    Eigen::MatrixXi subdividedFaces = subdividedMesh->getFaces();
    
    // Convert Eigen matrices to geometry-central format
    GC::DenseMatrix<double> vertsGC(subdividedVertices.rows(), 3);
    GC::DenseMatrix<int64_t> facesGC(subdividedFaces.rows(), 3);
    
    for (int i = 0; i < subdividedVertices.rows(); ++i) {
        for (int j = 0; j < 3; ++j) {
            vertsGC(i, j) = subdividedVertices(i, j);
        }
    }
    
    for (int i = 0; i < subdividedFaces.rows(); ++i) {
        for (int j = 0; j < 3; ++j) {
            facesGC(i, j) = static_cast<int64_t>(subdividedFaces(i, j));
        }
    }
    
    // Create geometry-central mesh
    std::unique_ptr<GCS::ManifoldSurfaceMesh> meshT = 
        std::make_unique<GCS::ManifoldSurfaceMesh>(facesGC);
    std::unique_ptr<GCS::VertexPositionGeometry> geometryT = 
        std::make_unique<GCS::VertexPositionGeometry>(*meshT);
    
    for (size_t i = 0; i < meshT->nVertices(); ++i) {
        for (int j = 0; j < 3; ++j) {
            geometryT->inputVertexPositions[i][j] = vertsGC(i, j);
        }
    }
    
    // Step 1: Build F3
    F3TileMesh f3(planarLocator);
    
    // Step 2: Find nearest vertices in F3
    F3Vertex startVertex = f3.nearestVertex(a);
    F3Vertex endVertex = f3.nearestVertex(b);
    
    // Step 3: Run Dijkstra on F3
    const Mesh* cuttedMesh = planarLocator.get_cutted_mesh_ptr().get();
    if (!cuttedMesh) {
        throw std::runtime_error("PlanarLocator has no cutted mesh");
    }
    const auto& uvPositions = *planarLocator.get_uv_positions_ptr();
    std::vector<F3Vertex> f3Path = dijkstraF3(f3, startVertex, endVertex, *cuttedMesh, uvPositions);
    
    if (f3Path.empty()) {
        std::string errorMsg = "No path found in F3. ";
        errorMsg += "Start vertex: F3 idx=" + std::to_string(startVertex.f3VertexIndex) + 
                    ", tile=(" + std::to_string(startVertex.tileU) + "," + std::to_string(startVertex.tileV) + 
                    "), baseIdx=" + std::to_string(startVertex.baseVertexIndex) + ". ";
        errorMsg += "End vertex: F3 idx=" + std::to_string(endVertex.f3VertexIndex) + 
                    ", tile=(" + std::to_string(endVertex.tileU) + "," + std::to_string(endVertex.tileV) + 
                    "), baseIdx=" + std::to_string(endVertex.baseVertexIndex) + ". ";
        errorMsg += "F3 has " + std::to_string(f3.numVertices()) + " vertices. ";
        
        std::vector<int> startNeighbors = f3.getNeighborIndices(startVertex.f3VertexIndex);
        std::vector<int> endNeighbors = f3.getNeighborIndices(endVertex.f3VertexIndex);
        errorMsg += "Start has " + std::to_string(startNeighbors.size()) + " neighbors. ";
        errorMsg += "End has " + std::to_string(endNeighbors.size()) + " neighbors.";
        
        throw std::runtime_error(errorMsg);
    }
    
    // Step 4: Convert F3 vertices to 2D UV coordinates
    Eigen::MatrixXd path2D(f3Path.size(), 2);
    
    for (size_t i = 0; i < f3Path.size(); ++i) {
        const F3Vertex& f3v = f3Path[i];
        // Get the UV position of the base vertex
        VertexCIter baseV = cuttedMesh->vertices.begin() + f3v.baseVertexIndex;
        const auto& baseUV = uvPositions[baseV];
        
        // Apply tile offset to get the UV coordinate in the tiled space
        path2D(i, 0) = baseUV.first + f3v.tileU;
        path2D(i, 1) = baseUV.second + f3v.tileV;
    }
    
    return path2D;
}

// Helper function to get neighbors with 2D UV weights
std::vector<std::pair<F3Vertex, double>> getNeighborsF3WithUV(
    const F3TileMesh& f3,
    const F3Vertex& vertex,
    const Mesh& cuttedMesh,
    const VertexData<std::pair<double, double>>& uvPositions) {
    
    std::vector<std::pair<F3Vertex, double>> neighbors;
    
    // Get neighbor F3 vertex indices
    std::vector<int> neighborIndices = f3.getNeighborIndices(vertex.f3VertexIndex);
    if (neighborIndices.empty()) {
        return neighbors;
    }
    
    // Get the UV position of the current vertex
    VertexCIter baseV = cuttedMesh.vertices.begin() + vertex.baseVertexIndex;
    const auto& baseUV = uvPositions[baseV];
    double vertexUVx = baseUV.first + vertex.tileU;
    double vertexUVy = baseUV.second + vertex.tileV;
    
    // For each neighbor F3 vertex, compute the 2D UV edge weight
    for (int neighborF3Idx : neighborIndices) {
        const F3Vertex& neighbor = f3.getVertex(neighborF3Idx);
        
        // Get the UV position of the neighbor vertex
        VertexCIter neighborBaseV = cuttedMesh.vertices.begin() + neighbor.baseVertexIndex;
        const auto& neighborBaseUV = uvPositions[neighborBaseV];
        double neighborUVx = neighborBaseUV.first + neighbor.tileU;
        double neighborUVy = neighborBaseUV.second + neighbor.tileV;
        
        // Compute Euclidean distance in 2D UV space
        double dx = neighborUVx - vertexUVx;
        double dy = neighborUVy - vertexUVy;
        double weight = std::sqrt(dx * dx + dy * dy);
        
        neighbors.emplace_back(neighbor, weight);
    }
    
    return neighbors;
}

} // namespace geomp

