#include "geomp/mesh/Edge.h"
#include "geomp/mesh/Mesh.h"

namespace geomp {

Edge::Edge(Mesh *mesh_):
onGenerator(false),
onCut(false),
isCuttable(true),
index(-1),
halfEdgeIndex(-1),
mesh(mesh_)
{

}

Edge::Edge(const Edge& e):
onGenerator(e.onGenerator),
onCut(e.onCut),
isCuttable(e.isCuttable),
generatorIndex(e.generatorIndex),
index(e.index),
referenceIndex(e.referenceIndex),
halfEdgeIndex(e.halfEdgeIndex),
mesh(e.mesh)
{

}

HalfEdgeIter Edge::halfEdge() const
{
	return mesh->halfEdges.begin() + halfEdgeIndex;
}

void Edge::setHalfEdge(HalfEdgeCIter he)
{
	halfEdgeIndex = he->index;
}

void Edge::setMesh(Mesh *mesh_)
{
	mesh = mesh_;
}

bool Edge::onBoundary() const
{
	return halfEdge()->onBoundary || halfEdge()->flip()->onBoundary;
}

std::pair<VertexCIter, VertexCIter> Edge::twoEndpoints() const
{
	VertexCIter u = halfEdge()->vertex();
	VertexCIter v = halfEdge()->flip()->vertex();
	if (u->index < v->index) {
		return std::make_pair(u, v);
	} else {
		return std::make_pair(v, u);
	}
}

std::ostream& operator<<(std::ostream& os, const Edge& e) 
{
	// get two endpoints of the edge
	std::pair<VertexCIter, VertexCIter> endpoints = e.twoEndpoints();
	os << "e_" << e.index << "=(" << endpoints.first->index << "," << endpoints.second->index << ")";
	return os;
}

} // namespace geomp
