#include "geomp/generators/Dijkstra.h"
#include <cstdint>
#include <limits>

namespace geomp {

std::pair<double, std::vector<int>> Dijkstra::findShortestPath(const Mesh& mesh, const VertexCIter& start, const VertexCIter& end)
{
    using VertexDistance = std::pair<double, VertexCIter>;
    using Comparator = std::greater<VertexDistance>;
    std::priority_queue<VertexDistance, std::vector<VertexDistance>, Comparator> pq;

    // Initialize distances
    VertexData<double> dist(mesh, std::numeric_limits<double>::infinity());
    VertexData<VertexCIter> prev(mesh, VertexCIter());
    VertexData<uint8_t> visited(mesh, 0);

    // Add the start vertex to the priority queue
    dist[start] = 0.0;
    pq.push({0.0, start});

    while (!pq.empty()) {
        // Get the vertex with the smallest distance
        auto [d, v] = pq.top();
        pq.pop();

        // Skip if the vertex has already been visited
        if (visited[v]) {
            continue;
        }

        // Mark the vertex as visited
        visited[v] = 1;

        // If the vertex is the end vertex, break
        if (v == end) {
            break;
        }

        // Iterate over all neighbors of v
        HalfEdgeCIter he = v->halfEdge();
        do {
            const VertexCIter& w = he->flip()->vertex();
            // Calculate the distance to the neighbor
            double alt = dist[v] + length(he->edge());
            // Update the distance to the neighbor if it is smaller
            if (alt < dist[w]) {
                dist[w] = alt;
                prev[w] = v;
                // Add the neighbor to the priority queue(allow duplicate vertices in the queue)
                pq.push({alt, w});
            }
            he = he->flip()->next();
        } while (he != v->halfEdge());
    }
    // Reconstruct path
    std::vector<int> path;
    VertexCIter v = end;
    while (v != start) {
        path.push_back(v->index);
        v = prev[v];
    }
    path.push_back(start->index);
    return {dist[end], path};
}

} // namespace geomp