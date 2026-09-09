#pragma once

#include "geomp/mesh/Types.h"
#include "geomp/mesh/MeshData.h"
namespace geomp {

class MeshSurgery {
  public:
    /**
     * Basic version of cutting the mesh along the edges with onCut = Tree
     *
     * @param mesh the mesh to cut in-place
     */
    static void cutMesh(Mesh& mesh);
    /**
     * Cutting the torus along the edges with onCut = Tree
     *
     * @param mesh the mesh to cut in-place
     * @return edges on the generators in the cutted mesh
     */
    static std::vector<EdgeCIter> cutTorus(Mesh& mesh);
    /**
     * Cutting the torus along the edges with onCut = Tree
     *
     * @param mesh the mesh to cut in-place
     * @param representative_edges_indices indices of the representative edges of the generators
     * @param representative_edges_after_cut representative edges and its cutted brothers in the cutted mesh
     */
    static void cutTorus(Mesh& mesh,
        std::vector<int> representative_edges_indices,
        std::vector<std::pair<int, int>>& representative_edges_pairs_after_cut);
    
    /**
     * Cutting the torus along the edges with onCut = Tree
     *
     * @param mesh the mesh to cut in-place
     * @return halfedges on the generators in the original mesh
     */
    static std::vector<HalfEdgeIter> cutTorus2(Mesh& mesh);

    /**
     * Tile the unit square along one side.
     *
     * @param mesh the mesh to tile in-place
     * @param direction the direction to tile (1: left, 2: right, 3: top, 4: bottom)
     */
    static void tileUnitSquareAlongOneSide(Mesh& mesh, int direction);
};

}  // namespace geomp