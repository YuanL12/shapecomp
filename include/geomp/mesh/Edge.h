#pragma once

#include "geomp/mesh/Types.h"

namespace geomp {

class Edge {
public:
	// constructor
	Edge(Mesh *mesh);

	// copy constructor
	Edge(const Edge& e);

	// returns one of the halfedges associated with this edge
	HalfEdgeIter halfEdge() const;

	// sets halfedge
	void setHalfEdge(HalfEdgeCIter he);

	// sets mesh
	void setMesh(Mesh *mesh);

	// checks if this edge is on the boundary
	bool onBoundary() const;

	// return two endpoints of the edge (u, v) where u < v
	std::pair<VertexCIter, VertexCIter> twoEndpoints() const;

	// boolean flag to indicate if edge is on a generator
	bool onGenerator = false;

	// boolean flag to indicate if edge is on a cut
	bool onCut = false;

	// boolean flag to indicate if cut can pass through edge
	bool isCuttable = true;

	// (for torus), int flag to indicate: 1 if edge is on a generator 1, 
	// 2 if edge is on a generator 2, -1 otherwise
	int generatorIndex = -1;

	// id between 0 and |E|-1
	int index;
	
	// index of cutted edge to the original one 
	int referenceIndex = -1;

	// support cout for debugging
	friend std::ostream& operator<<(std::ostream& os, const Edge& e);

private:
	// index of one of the halfedges associated with this edge
	int halfEdgeIndex;

	// pointer to mesh this edge belongs to
	Mesh *mesh;
};

} // namespace geomp
