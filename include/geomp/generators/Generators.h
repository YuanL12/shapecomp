#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <cstdint>
#include <queue>
#include <unordered_set>

// basic mesh 
#include "geomp/mesh/MeshData.h"
#include "geomp/mesh/Types.h"

// matrix operations
#include "geomp/matrix/Laplacian.h"

// homology generators operations
#include "geomp/generators/Dijkstra.h"
#include "geomp/generators/MeshSurgery.h"

namespace geomp {
/*
Class Linked List
Class used to store vertices on a generator.
  support
  1. next: return the next vertex Node on the generator, if exists, and None otherwise
  2. prev: return the previous vertex Node on the generator, if exists, and None otherwise
  3. ++: same as next
  4. --: same as prev
  5. insert: given a edge,
          1) find it two vertices u, v
          2) find one of u/v that still does not have prev or next assigned, say u->next() is emepty
          3) assign the next vertex node of as v, u->setNext(v)
          Remark: in the second step, when trying to insert edge e=(u,v),
                  except the first insertion, at least one of u/v has appeared.
                  Here is algorithm:
                  Iterat over all edges in generator:
                          a. If neither u and v appears, u->next() = v. Continue.
                          b. Find the appeared vertex (u,v vertices can both appear in the past);
                          c. Check where the spot (prev or next) of appeared vertex is available,
                          c. Connect at the available spot;
*/
class VertexNode {
  public:
    VertexIter vertexPtr;
    VertexNode* prev;
    VertexNode* next;

    VertexNode(VertexIter v) : vertexPtr(v), prev(nullptr), next(nullptr) {}

    inline bool hasPrev() const { return prev != nullptr; };
    inline bool hasNext() const { return next != nullptr; };

    // operators
    inline VertexNode* operator++() { return next; };
    inline VertexNode* operator--() { return prev; };
    inline bool operator!=(const VertexNode* other) const { return this != other; };
    // // dereference
    // inline VertexIter operator*() { return vertexPtr;};
    // // -> operator
    // inline VertexIter operator->() { return vertexPtr;};
};

class CircularLinkedList {
  public:
    // constructor
    CircularLinkedList() : head(nullptr), tail(nullptr) {}

    // constructor by a vector of edges(orientation might be wrong)
    CircularLinkedList(const std::vector<EdgeCIter>& edges);

    // constructor by a vector of ordered halfedges (assume the input a vector connected halfedges)
    CircularLinkedList(const std::vector<HalfEdgeIter>& halfEdges);

    // destructor
    ~CircularLinkedList();

    // reverse the direction of the whole linked list
    void reverse();

    // begin return head
    inline VertexNode* begin() const { return head; };
    inline VertexNode* begin() { return head; };
    inline VertexNode* getHead() const { return head; };
    inline VertexNode* getHead() { return head; };

    // set head
    inline void setHead(VertexNode* v) {
        head = v;
        tail = v->prev;
    };

    // insert
    void insert(EdgeCIter e);

    void orientedInsert(HalfEdgeCIter he);

    // next
    VertexNode* next(VertexNode* v);

    // prev
    VertexNode* prev(VertexNode* v);

    // print
    void print() const;

  private:
    VertexNode* head = nullptr;
    VertexNode* tail = nullptr;

    VertexNode* findNode(VertexCIter vertex);
};

class Generators {
  public:
    // Constructor
    Generators(const Mesh& mesh);

    // computes generators
    void labelGeneratorsOnEdgesRaw(Mesh& mesh) const;

    /**
     * Find the edges on homological generators by Tree-Cotree algorithm and label them onCut
     *
     * @param mesh the mesh to label in-place
     * @return representative edges of the generators (left edges in Tree-Cotree algorithm)
     */
    std::vector<int> labelTreeCotreeGeneratorsEdges(Mesh& mesh) const;

    /**
     * Erickson–Whittlesey algorithm.
     * First find and cut the mesh on the generators by Tree-Cotree algorithm,
     * then find the shortest path (Dijkstra's algorithm) between the two endpoints of the representative edges and
     * label them onCut
     *
     * @param mesh the mesh to label in-place
     */
    void labelShortestGeneratorEdges(Mesh& mesh);

    /**
     * Find the edges on homological generators by Reeb graph algorithm and label them onCut
     *
     * @param mesh the mesh to label in-place
     */
    void labelReebGraphGeneratorEdges(Mesh& mesh);

    /**
     * Find the edges on homological generators by Reeb graph algorithm and label them onCut
     *
     * @param mesh the mesh to label in-place
     * @param distinctDirection the direction vector to produce a Reeb graph
     */
    void labelReebGraphGeneratorEdges(Mesh& mesh, const std::vector<double>& distinctDirection);

    /**
     * Find the generators that intersect only at one vertex by barycentric subdivision.
     * label the edges on the first generator as onCut.
     *
     * @param mesh the mesh to find in-place
     * @return faces on the second generator
     */
    std::vector<FaceCIter> findOnceIntersectGenerators(Mesh& mesh);

    // find closed 1-form, i.e. a functional on each edge, where edge is oriented from u to v (u < v).
    std::vector<Eigen::VectorXd> computeHarmonicOneForms(const Mesh& mesh) const;

  private:
    VertexData<VertexCIter> primalParent;  // primal parent
    FaceData<FaceCIter> dualParent, ngon;  // dual parent, ngon

    // find the closed 1-form corresponding to an edge e (in a homological generator)
    EdgeData<double> findClosedOneForm(const Mesh& mesh, EdgeCIter starting_edge) const;

    // find the shortest path between two nodes in a dual spanning tree
    std::vector<FaceCIter> findShortestPathInDualTree(FaceCIter f1, FaceCIter f2) const;

    // find the shortest path between two nodes in a primal spanning tree, return edges
    std::vector<EdgeIter> findShortestPathInSpanningTree(VertexCIter v1, VertexCIter v2) const;

    // builds primal spanning tree
    void buildPrimalSpanningTree(const Mesh& mesh);

    // checks whether an edge is in the primal spanning tree
    bool inPrimalSpanningTree(EdgeCIter e) const;

    // builds dual spanning tree
    void buildDualSpanningTree(const Mesh& mesh);

    // checks whether an edge is in the dual spanning tree
    bool inDualSpanningTree(EdgeCIter e) const;

    // sets new faces for a face f by passing halfedges h1, h2, h3
    void setNewFaces(FaceIter f, HalfEdgeIter h1, HalfEdgeIter h2, HalfEdgeIter h3);
};

// reset the head of nodesOnGenerators to the root of 2 generators of a torus
void resetHeadForTorus(CircularLinkedList& nodesOnGenerators);

// returns shared edge between u and v
EdgeIter sharedEdge(VertexCIter u, VertexCIter v);

}  // namespace geomp
