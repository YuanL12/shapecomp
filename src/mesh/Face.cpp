#include "geomp/mesh/Face.h"
#include "geomp/mesh/Mesh.h"

namespace geomp {

Face::Face(Mesh *mesh_):
fillsHole(false),
inNorthPoleVicinity(false),
index(-1),
halfEdgeIndex(-1),
mesh(mesh_)
{

}

Face::Face(const Face& f):
fillsHole(f.fillsHole),
inNorthPoleVicinity(f.inNorthPoleVicinity),
index(f.index),
halfEdgeIndex(f.halfEdgeIndex),
mesh(f.mesh)
{

}

Vector Face::barycenter() const
{
	Vector b(0.0, 0.0, 0.0);
	HalfEdgeCIter h = halfEdge();
	int n = 0;
	do {
		b += h->vertex()->position;
		h = h->next();
		n++;
	} while (h != halfEdge());

	return b/n;
}

HalfEdgeIter Face::halfEdge() const
{
	return mesh->halfEdges.begin() + halfEdgeIndex;
}

void Face::setHalfEdge(HalfEdgeCIter he)
{
	halfEdgeIndex = he->index;
}

void Face::setMesh(Mesh *mesh_)
{
	mesh = mesh_;
}

bool Face::isReal() const
{
	return !halfEdge()->onBoundary && !inNorthPoleVicinity;
}

} // namespace geomp
