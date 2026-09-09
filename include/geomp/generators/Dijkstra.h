#pragma once
#include "geomp/mesh/MeshData.h"
#include "geomp/mesh/Types.h"
#include <queue>
#include <vector>
#include <utility>
#include <unordered_map>

namespace geomp {
class Dijkstra {
  public:
    // return the shortest path (indices of vertices) and the distance between the start and end vertices
    static std::pair<double, std::vector<int>> findShortestPath(const Mesh& mesh, const VertexCIter& start, const VertexCIter& end);
};

}  // namespace geomp