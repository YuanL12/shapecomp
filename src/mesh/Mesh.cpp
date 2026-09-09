#include "geomp/mesh/Mesh.h"

namespace geomp {

Mesh::Mesh():
radius(0.0),
cm(0.0, 0.0, 0.0),
status(ErrorCode::ok)
{

}

Mesh::Mesh(const Mesh& mesh):
radius(mesh.radius),
cm(mesh.cm),
status(mesh.status)
{
	// std::cout << "Call the copy constructor of Mesh" << std::endl;
	// initialize halfEdges
	halfEdges.reserve(mesh.halfEdges.capacity());
	for (HalfEdgeCIter h = mesh.halfEdges.begin(); h != mesh.halfEdges.end(); h++) {
		HalfEdgeIter hNew = halfEdges.emplace(halfEdges.end(), HalfEdge(*h));
		hNew->setMesh(this);
	}

	// initialize vertices
	vertices.reserve(mesh.vertices.capacity());
	for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
		VertexIter vNew = vertices.emplace(vertices.end(), Vertex(*v));
		vNew->setMesh(this);
	}

	// initialize edges
	edges.reserve(mesh.edges.capacity());
	for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
		EdgeIter eNew = edges.emplace(edges.end(), Edge(*e));
		eNew->setMesh(this);
	}

	// initialize faces
	faces.reserve(mesh.faces.capacity());
	for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
		FaceIter fNew = faces.emplace(faces.end(), Face(*f));
		fNew->setMesh(this);
	}

	// initialize corners
	corners.reserve(mesh.corners.capacity());
	for (CornerCIter c = mesh.corners.begin(); c != mesh.corners.end(); c++) {
		CornerIter cNew = corners.emplace(corners.end(), Corner(*c));
		cNew->setMesh(this);
	}

	// initialize boundaries
	boundaries.reserve(mesh.boundaries.capacity());
	for (BoundaryCIter b = mesh.boundaries.begin(); b != mesh.boundaries.end(); b++) {
		BoundaryIter bNew = boundaries.emplace(boundaries.end(), Face(*b));
		bNew->setMesh(this);
	}
}

Mesh::Mesh(Mesh&& mesh) noexcept
    : radius(std::move(mesh.radius)),
      cm(std::move(mesh.cm)),
      status(std::move(mesh.status)),
      halfEdges(std::move(mesh.halfEdges)),
      vertices(std::move(mesh.vertices)),
      edges(std::move(mesh.edges)),
      faces(std::move(mesh.faces)),
      corners(std::move(mesh.corners)),
      boundaries(std::move(mesh.boundaries)) {
    // std::cout << "Call the move constructor of Mesh" << std::endl;

    // Update the mesh pointers for all elements in the moved vectors
    for (HalfEdgeIter h = halfEdges.begin(); h != halfEdges.end(); ++h) {
        h->setMesh(this);
    }

    for (VertexIter v = vertices.begin(); v != vertices.end(); ++v) {
        v->setMesh(this);
    }

    for (EdgeIter e = edges.begin(); e != edges.end(); ++e) {
        e->setMesh(this);
    }

    for (FaceIter f = faces.begin(); f != faces.end(); ++f) {
        f->setMesh(this);
    }

    for (CornerIter c = corners.begin(); c != corners.end(); ++c) {
        c->setMesh(this);
    }

    for (BoundaryIter b = boundaries.begin(); b != boundaries.end(); ++b) {
        b->setMesh(this);
    }
}

int Mesh::eulerCharacteristic() const
{
	return (int)(vertices.size() - edges.size() + faces.size());
}

double Mesh::diameter() const
{
	double maxLimit = std::numeric_limits<double>::max();
	double minLimit = std::numeric_limits<double>::lowest();
	Vector minBounds(maxLimit, maxLimit, maxLimit);
	Vector maxBounds(minLimit, minLimit, minLimit);

	for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
		const Vector& p = v->position;

		minBounds.x = std::min(p.x, minBounds.x);
		minBounds.y = std::min(p.y, minBounds.y);
		minBounds.z = std::min(p.z, minBounds.z);
		maxBounds.x = std::max(p.x, maxBounds.x);
		maxBounds.y = std::max(p.y, maxBounds.y);
		maxBounds.z = std::max(p.z, maxBounds.z);
	}

	return (maxBounds - minBounds).norm();
}

void computeEigenvectors2x2(double a, double b, double c, Vector& v1, Vector& v2)
{
	double disc = std::sqrt((a - c)*(a - c) + 4.0*b*b)/2.0;
	double lambda1 = (a + c)/2.0 - disc;
	double lambda2 = (a + c)/2.0 + disc;

	v1 = b < 0.0 ? Vector(-b, a - lambda1) : Vector(b, lambda2 - a);
	double v1Norm = v1.norm();
	if (v1Norm > 0) v1 /= v1Norm;
	v2 = Vector(-v1.y, v1.x);
}

void Mesh::projectUvsToPcaAxis()
{
	// compute center of mass
	Vector cm;
	int nUvs = 0;
	for (FaceCIter f = faces.begin(); f != faces.end(); f++) {
		if (f->isReal() && !f->fillsHole) {
			Vector centroid = centroidUV(f);

			cm += centroid;
			nUvs++;
		}
	}
	cm /= nUvs;

	// translate UVs to origin
	for (WedgeIter w = wedges().begin(); w != wedges().end(); w++) {
		if (w->isReal()) {
			w->uv -= cm;
		}
	}

	// build covariance matrix
	double a = 0, b = 0, c = 0;
	for (FaceCIter f = faces.begin(); f != faces.end(); f++) {
		if (f->isReal() && !f->fillsHole) {
			Vector centroid = centroidUV(f);

			a += centroid.x*centroid.x;
			b += centroid.x*centroid.y;
			c += centroid.y*centroid.y;
		}
	}

	// compute eigenvectors
	Vector v1, v2;
	computeEigenvectors2x2(a, b, c, v1, v2);

	// project uvs onto principal axes
	for (WedgeIter w = wedges().begin(); w != wedges().end(); w++) {
		if (w->isReal()) {
			Vector& uv = w->uv;
			uv = Vector(dot(v1, uv), dot(v2, uv));
			uv += cm;
		}
	}
}

void Mesh::orientUvsToMinimizeBoundingBox(int nRotations)
{
	// precompute rotations
	double maxRotation = M_PI;
	double minInf = -std::numeric_limits<double>::infinity();
	double maxInf = std::numeric_limits<double>::infinity();
	std::vector<Vector> boxMin(nRotations, Vector(maxInf, maxInf));
	std::vector<Vector> boxMax(nRotations, Vector(minInf, minInf));

	std::vector<Vector> rotations(nRotations);
	for (int i = 0; i < nRotations; i++) {
		double theta = (maxRotation*i)/nRotations;
		rotations[i] = Vector(std::cos(theta), std::sin(theta));
	}

	// try all rotations
	for (WedgeIter w = wedges().begin(); w != wedges().end(); w++) {
		if (w->isReal()) {
			const Vector& uv = w->uv;

			for (int i = 0; i < nRotations; i++) {
				double cosTheta = rotations[i].x;
				double sinTheta = rotations[i].y;
				Vector uvRotated(cosTheta*uv.x - sinTheta*uv.y,
								 sinTheta*uv.x + cosTheta*uv.y);
				boxMin[i].x = std::min(boxMin[i].x, uvRotated.x);
				boxMin[i].y = std::min(boxMin[i].y, uvRotated.y);
				boxMax[i].x = std::max(boxMax[i].x, uvRotated.x);
				boxMax[i].y = std::max(boxMax[i].y, uvRotated.y);
			}
		}
	}

	// find the best rotation
	Vector bestRotation;
	double minScore = std::numeric_limits<double>::infinity();
	for (int i = 0; i < nRotations; i++) {
		Vector extent = boxMax[i] - boxMin[i];
		if (extent.y < minScore) {
			minScore = extent.y;
			bestRotation = rotations[i];
		}
	}

	// apply best rotation
	double cosTheta = bestRotation.x;
	double sinTheta = bestRotation.y;
	for (WedgeIter w = wedges().begin(); w != wedges().end(); w++) {
		if (w->isReal()) {
			const Vector& uv = w->uv;
			Vector uvRotated(cosTheta*uv.x - sinTheta*uv.y,
							 sinTheta*uv.x + cosTheta*uv.y);
			w->uv = uvRotated;
		}
	}
}

EdgeCIter Mesh::getEdge(int u_id, int v_id) const{
	// search by vertex indices first 
	auto u = vertices[u_id];
	HalfEdgeCIter he = u.halfEdge();
	do {
		if (he->flip()->vertex()->index == v_id) {
			return he->edge();
		}
		he = he->flip()->next();
	} while (he != u.halfEdge());

	// if not found, search by iterating over all edges 
	for (EdgeCIter e = edges.begin(); e != edges.end(); e++) {
		int e_0_id = e->halfEdge()->vertex()->index;
		int e_1_id = e->halfEdge()->flip()->vertex()->index;
		if ((e_0_id == u_id && e_1_id == v_id) || (e_0_id == v_id && e_1_id == u_id)) {
			return e;
		}
	}
	
	std::cerr << "Trying to find edge with endpoints u, v, but not found, will stop now!" << std::endl;
	std::abort();
	return edges.end();
}

EdgeIter Mesh::getEdge(int u_id, int v_id) {
	// search by vertex indices first 
	auto u = vertices[u_id];
	HalfEdgeCIter he = u.halfEdge();
	do {
		if (he->flip()->vertex()->index == v_id) {
			return he->edge();
		}
		he = he->flip()->next();
	} while (he != u.halfEdge());

	// if not found, search by iterating over all edges 
	for (EdgeIter e = edges.begin(); e != edges.end(); e++) {
		int e_0_id = e->halfEdge()->vertex()->index;
		int e_1_id = e->halfEdge()->flip()->vertex()->index;
		if ((e_0_id == u_id && e_1_id == v_id) || (e_0_id == v_id && e_1_id == u_id)) {
			return e;
		}
	}
	
	std::cerr << "Trying to find edge with endpoints u, v, but not found, will stop now!" << std::endl;
	std::abort();
	return edges.end();
}

double Mesh::areaRatio() const
{
	double totalArea = 0.0;
	double totalAreaUV = 0.0;
	for (FaceCIter f = faces.begin(); f != faces.end(); f++) {
		if (f->isReal() && !f->fillsHole) {
			totalArea += area(f);
			totalAreaUV += areaUV(f);
		}
	}

	if (std::isinf(totalAreaUV) || std::isnan(totalAreaUV)) return 1.0;
	return totalAreaUV > 0.0 ? totalArea/totalAreaUV : 1.0;
}

CutPtrSet Mesh::cutBoundary() const
{
	if (boundaries.size() == 0) {
		// if there is no boundary, initialize the iterator with the first edge
		// on the cut
		for (EdgeCIter e = edges.begin(); e != edges.end(); e++) {
			if (e->onCut) return CutPtrSet(e->halfEdge());
		}

		return CutPtrSet();
	}

	return CutPtrSet(boundaries[0].halfEdge());
}

std::vector<Wedge>& Mesh::wedges()
{
	return corners;
}

const std::vector<Wedge>& Mesh::wedges() const
{
	return corners;
}

Eigen::MatrixXd Mesh::getVertices() const
{
	Eigen::MatrixXd V(vertices.size(), 3);
	for (size_t i = 0; i < vertices.size(); ++i) {
		auto idx = vertices[i].index;
		V(idx, 0) = vertices[i].position.x;
		V(idx, 1) = vertices[i].position.y;
		V(idx, 2) = vertices[i].position.z;
	}
	return V;
}

Eigen::MatrixXi Mesh::getFaces() const
{
	Eigen::MatrixXi F(faces.size(), 3);
	for (size_t i = 0; i < faces.size(); ++i) {
		auto idx = faces[i].index;
		F(idx, 0) = faces[i].halfEdge()->vertex()->index;
		F(idx, 1) = faces[i].halfEdge()->next()->vertex()->index;
		F(idx, 2) = faces[i].halfEdge()->prev()->vertex()->index;
	}
	return F;
}

Eigen::MatrixXi Mesh::getEdges() const
{
	Eigen::MatrixXi E(edges.size(), 2);
	for (size_t i = 0; i < edges.size(); ++i) {
		auto idx = edges[i].index;
		E(idx, 0) = edges[i].halfEdge()->vertex()->index;
		E(idx, 1) = edges[i].halfEdge()->flip()->vertex()->index;
	}
	return E;
}

} // namespace geomp
