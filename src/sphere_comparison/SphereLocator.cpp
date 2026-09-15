#include "geomp/sphere_comparison/SphereLocator.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <tuple>

#include "geomp/embedding/sphere/Normalize.h"
#include "geomp/embedding/sphere/SphereTutte.h"
#include "geomp/mesh/GeometryUtils.h"

namespace geomp {

SphereLocator::SphereLocator() : mesh_ptr(nullptr), pole_index(-1) {}

SphereLocator::SphereLocator(Mesh& mesh, std::optional<std::string> Tutte_embedding_type)
	: mesh_ptr(nullptr), pole_index(-1) {
	construct(mesh, Tutte_embedding_type, 0);
}

void SphereLocator::construct(Mesh& mesh, std::optional<std::string> Tutte_embedding_type, int pole_type) {
	if (!mesh.boundaries.empty()) {
		throw std::runtime_error("SphereLocator requires a closed mesh (no boundary loops)");
	}
	if (mesh.eulerCharacteristic() != 2) {
		throw std::runtime_error(
			"SphereLocator requires a genus-zero mesh (Euler characteristic 2), got " +
			std::to_string(mesh.eulerCharacteristic()));
	}

	mesh_ptr = &mesh;
	int Tutte_type = 2;  // coTan default
	if (Tutte_embedding_type.has_value()) {
		Tutte_type = get_Tutte_embedding_type(Tutte_embedding_type.value());
	}

	SphereTutte tutte;
	SphereEmbedding embedding = tutte.embed(mesh, pole_type, Tutte_type, false);
	sphere_vertices = std::move(embedding.vertices);
	pole_index = embedding.poleIndex;
	SphereNormalize::normalize(sphere_vertices, false);
}

int SphereLocator::get_Tutte_embedding_type(const std::string& type) const {
	if (type == "Uniform") {
		return 1;
	}
	if (type == "coTan") {
		return 2;
	}
	if (type == "MeanValue") {
		return 3;
	}
	if (type == "Authalic") {
		return 4;
	}
	throw std::runtime_error("Invalid Tutte embedding type: " + type);
}

Mesh* SphereLocator::get_mesh_ptr() const {
	return mesh_ptr;
}

int SphereLocator::get_pole_index() const {
	return pole_index;
}

Eigen::MatrixXd SphereLocator::get_sphere_vertices() const {
	if (!mesh_ptr) {
		throw std::runtime_error("SphereLocator has no mesh");
	}
	return sphere_vertices;
}

Vector SphereLocator::sphereVertex(int index) const {
	return Vector(sphere_vertices(index, 0), sphere_vertices(index, 1), sphere_vertices(index, 2));
}

bool SphereLocator::barycentricOnFace(
	FaceCIter f,
	const Eigen::Vector3d& p,
	Eigen::Vector3d& bary,
	Eigen::Vector3d& raw) const {
	if (!f->isReal()) {
		return false;
	}
	const int ia = f->halfEdge()->vertex()->index;
	const int ib = f->halfEdge()->next()->vertex()->index;
	const int ic = f->halfEdge()->prev()->vertex()->index;
	const Vector A = sphereVertex(ia);
	const Vector B = sphereVertex(ib);
	const Vector C = sphereVertex(ic);

	const Vector u = cross(B - A, C - A);
	const Vector pvec(p[0], p[1], p[2]);
	const double n = dot(u, pvec);
	if (std::abs(n) < 1e-15) {
		return false;
	}
	const double alpha = dot(u, A) / n;
	if (alpha <= 0.0) {
		return false;
	}
	const Vector p2 = Vector(alpha * p[0], alpha * p[1], alpha * p[2]);
	double la, lb, lc;
	std::tie(la, lb, lc) = barycentric_coordinates(A, B, C, p2);
	raw = Eigen::Vector3d(la, lb, lc);

	bary = raw;
	bary[0] = std::max(0.0, bary[0]);
	bary[1] = std::max(0.0, bary[1]);
	bary[2] = std::max(0.0, bary[2]);
	const double s = bary.sum();
	if (s <= std::numeric_limits<double>::epsilon()) {
		int argmax = 0;
		if (raw[1] > raw[argmax]) {
			argmax = 1;
		}
		if (raw[2] > raw[argmax]) {
			argmax = 2;
		}
		bary = Eigen::Vector3d::Zero();
		bary[argmax] = 1.0;
		return true;
	}
	bary /= s;
	return true;
}

std::pair<int, Eigen::Vector3d> SphereLocator::findBaryCoordsBrute(const Eigen::Vector3d& p) const {
	int best_face = -1;
	Eigen::Vector3d best_bary = Eigen::Vector3d(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
	double best_score = -std::numeric_limits<double>::infinity();

	for (FaceCIter f = mesh_ptr->faces.begin(); f != mesh_ptr->faces.end(); ++f) {
		Eigen::Vector3d bary, raw;
		if (!barycentricOnFace(f, p, bary, raw)) {
			continue;
		}
		const double score = std::min(raw[0], std::min(raw[1], raw[2]));
		if (score > best_score) {
			best_score = score;
			best_face = f->index;
			best_bary = bary;
		}
	}

	if (best_face < 0) {
		// Degenerate fallback: nearest face centroid.
		double best_dist = std::numeric_limits<double>::infinity();
		for (FaceCIter f = mesh_ptr->faces.begin(); f != mesh_ptr->faces.end(); ++f) {
			if (!f->isReal()) {
				continue;
			}
			const int ia = f->halfEdge()->vertex()->index;
			const int ib = f->halfEdge()->next()->vertex()->index;
			const int ic = f->halfEdge()->prev()->vertex()->index;
			Eigen::Vector3d c =
				(sphere_vertices.row(ia) + sphere_vertices.row(ib) + sphere_vertices.row(ic)).transpose() / 3.0;
			const double dist = (c - p).squaredNorm();
			if (dist < best_dist) {
				best_dist = dist;
				best_face = f->index;
				best_bary = Eigen::Vector3d(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
			}
		}
	}
	if (best_face < 0) {
		throw std::runtime_error("SphereLocator: unable to locate a spherical triangle");
	}
	return {best_face, best_bary};
}

std::pair<int, Eigen::Vector3d> SphereLocator::findBaryCoords(const Eigen::Vector3d& point) const {
	if (!mesh_ptr) {
		throw std::runtime_error("SphereLocator has no mesh");
	}
	if (sphere_vertices.rows() == 0) {
		throw std::runtime_error("SphereLocator has no sphere vertices");
	}

	Eigen::Vector3d p = point;
	const double nrm = p.norm();
	if (nrm <= 1e-15) {
		throw std::runtime_error("SphereLocator: query point is at the origin");
	}
	p /= nrm;

	FaceCIter f = mesh_ptr->faces.begin();
	while (f != mesh_ptr->faces.end() && !f->isReal()) {
		++f;
	}
	if (f == mesh_ptr->faces.end()) {
		throw std::runtime_error("SphereLocator: mesh has no real faces");
	}

	const int maxWalk = static_cast<int>(mesh_ptr->faces.size()) * 4 + 8;
	int count = 0;
	while (count < maxWalk) {
		if (!f->isReal()) {
			break;
		}
		Eigen::Vector3d bary, raw;
		const bool ok = barycentricOnFace(f, p, bary, raw);
		if (ok && raw[0] >= -1e-8 && raw[1] >= -1e-8 && raw[2] >= -1e-8) {
			return {f->index, bary};
		}

		if (!ok) {
			break;
		}

		HalfEdgeCIter AB = f->halfEdge();
		HalfEdgeCIter BC = AB->next();
		HalfEdgeCIter CA = AB->prev();
		if (raw[0] <= raw[1] && raw[0] <= raw[2]) {
			f = BC->flip()->face();
		} else if (raw[1] <= raw[2]) {
			f = CA->flip()->face();
		} else {
			f = AB->flip()->face();
		}
		++count;
	}

	return findBaryCoordsBrute(p);
}

std::pair<int, Eigen::Vector3d> SphereLocator::findBaryCoords(const std::vector<double>& point) const {
	if (point.size() != 3) {
		throw std::runtime_error("SphereLocator.findBaryCoords expects a length-3 point");
	}
	return findBaryCoords(Eigen::Vector3d(point[0], point[1], point[2]));
}

std::vector<std::pair<int, Eigen::Vector3d>> SphereLocator::findBaryCoordsParallel(
	const Eigen::MatrixXd& points) const {
	if (!mesh_ptr) {
		throw std::runtime_error("SphereLocator has no mesh");
	}
	if (points.cols() != 3) {
		throw std::runtime_error("SphereLocator.findBaryCoordsParallel expects an (N, 3) array");
	}

	const size_t n = static_cast<size_t>(points.rows());
	std::vector<std::pair<int, Eigen::Vector3d>> results(n);

	tf::Executor executor;
	tf::Taskflow taskflow;
	taskflow.for_each_index(static_cast<size_t>(0), n, static_cast<size_t>(1), [&](size_t i) {
		Eigen::Vector3d p(points(static_cast<int>(i), 0), points(static_cast<int>(i), 1), points(static_cast<int>(i), 2));
		results[i] = findBaryCoords(p);
	});
	executor.run(taskflow).get();
	return results;
}

}  // namespace geomp
