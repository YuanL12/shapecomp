#pragma once
#include <vector>
#include "geomp/mesh/Mesh.h"

// ReebHanTun includes
#include <ReebHanTun/SimpleMesh.h>
#include <ReebHanTun/psbmReebGraph.h>
#include <ReebHanTun/RenderVector3.h>
#include <map>
#include <set>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Forward declarations for ReebHanTun types
class _SimpleMesh;
class _SimpleMeshVertex;
class _SimpleMeshEdge;
class _SimpleMeshTriangle;
class psbmReebGraph;
class Vector3;
struct myPairCompare;

namespace geomp {

/**
 * Convert geomp::Mesh to _SimpleMesh
 * Maintains vertex index alignment so we can map edges back using vertex indices
 *
 * @param mesh Input geomp::Mesh
 * @param simpleMesh Output _SimpleMesh (will be cleared and populated)
 * @param orientTriangles Output vector of oriented triangle vertex indices (3 per triangle)
 * @param extraVertices Output set of extra vertex indices (for hole closing, initially empty)
 * @param meshNormal Output vector of triangle normals (one per triangle, must persist for lifetime of simpleMesh)
 */
void convertMeshToSimpleMesh(const Mesh& mesh,
    _SimpleMesh& simpleMesh,
    std::vector<int>& orientTriangles,
    std::set<int>& extraVertices,
    std::vector<Vector3>& meshNormal);

/**
 * Find Reeb graph generators for a mesh
 * Uses ReebHanTun to compute handle and tunnel loops, then returns the vertices on the generators
 *
 * @param mesh Input mesh
 * @return Pair of vectors containing the vertices on the handle and tunnel generators
 */

class ReebHanTunLoops {
  private:
    _SimpleMesh simpleMesh;
    std::vector<int> orientTriangles;
    std::set<int> extraVertices;
    std::vector<Vector3> meshNormal;

  public:
    std::vector<std::set<int>> h_basis_loops;  // handle loops
    std::vector<std::set<int>> v_basis_loops;  // tunnel loops

    ReebHanTunLoops() = delete;
    ReebHanTunLoops(const Mesh& mesh);
    // Constructor with a defined direction vector(it should produce a Morse function when product with the vertices)
    ReebHanTunLoops(const Mesh& mesh, const Vector3& distinctDirection);

    /**
     * Label the edges on the generators by Reeb graph algorithm and label them onCut on the original mesh
     *
     * @param mesh Input/output mesh - edges will be labeled with onCut and generatorIndex
     */
    void labelReebGraphGeneratorsOnMesh(Mesh& mesh);

    /**
     * Find the vertices on the generators by Reeb graph algorithm
     * Uses ReebHanTun to compute handle and tunnel loops, then returns the vertices on the generators
     *
     * @return Pair of vectors containing the vertices on the handle and tunnel generators
     */
    std::pair<std::vector<int>, std::vector<int>> findReebGraphGenerators();
};

}  // namespace geomp
