#include "geomp/embedding/sphere/SphereTutte.h"

#include <iostream>
#include <stdexcept>
#include <string>

#include <Eigen/SparseLU>

namespace geomp {

namespace {
constexpr double kPi = 3.14159265358979323846;
}

VertexIter SphereTutte::pickPole(Mesh& mesh, int type) const {
	if (type == 0) {
		return mesh.vertices.begin();
	}

	int val_ref = mesh.vertices.begin()->degree();
	VertexIter pole = mesh.vertices.begin();
	for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
		const int valence = v->degree();
		if (type == 2) {
			if (valence < val_ref) {
				val_ref = valence;
				pole = v;
			}
		} else if (type == 1) {
			if (valence > val_ref) {
				val_ref = valence;
				pole = v;
			}
		}
	}
	return pole;
}

bool SphereTutte::isPole(const Vertex& v) const {
	return v.index == poleIndex;
}

bool SphereTutte::isInterior(const Vertex& v) const {
	return !isPole(v) && !v.onBoundaryTutte;
}

void SphereTutte::punctureMesh(Mesh& mesh, VertexIter pole) {
	for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
		v->inNorthPoleVicinity = false;
		v->onBoundaryTutte = false;
		v->indexN = -1;
	}
	for (FaceIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f) {
		f->inNorthPoleVicinity = false;
	}

	pole->inNorthPoleVicinity = true;
	int idx = 0;
	HalfEdgeIter he = pole->halfEdge();
	do {
		he->face()->inNorthPoleVicinity = true;
		VertexIter ring = he->flip()->vertex();
		ring->indexN = idx;
		ring->onBoundaryTutte = true;
		ring->inNorthPoleVicinity = true;
		++idx;
		he = he->flip()->next();
	} while (he != pole->halfEdge());

	Nb = idx;

	idx = 0;
	for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
		if (isInterior(*v)) {
			v->indexN = idx;
			++idx;
		}
	}
}

void SphereTutte::restorePunctureFlags(Mesh& mesh) {
	for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
		v->inNorthPoleVicinity = false;
		v->onBoundaryTutte = false;
	}
	for (FaceIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f) {
		f->inNorthPoleVicinity = false;
	}
}

void SphereTutte::initSystem(Mesh& mesh) {
	if (H_is_set) {
		throw std::runtime_error("H matrix is already set when trying to initialize it");
	}

	std::vector<Triplet> Mat_coefficients;
	const double zero = 0.0;

	for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
		if (isInterior(*v)) {
			const int idx = v->indexN;
			Mat_coefficients.push_back(Triplet(idx, idx, zero));
			Bx[idx] = 0.0;
			By[idx] = 0.0;
		}
	}

	for (EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); ++e) {
		HalfEdgeIter he = e->halfEdge();
		HalfEdgeIter he2 = he->flip();
		VertexIter v_i = he->vertex();
		VertexIter v_j = he2->vertex();
		if (isInterior(*v_i) && isInterior(*v_j)) {
			Mat_coefficients.push_back(Triplet(v_i->indexN, v_j->indexN, zero));
			Mat_coefficients.push_back(Triplet(v_j->indexN, v_i->indexN, zero));
		}
	}

	H.setFromTriplets(Mat_coefficients.begin(), Mat_coefficients.end());
	H_is_set = true;
}

void SphereTutte::resetSystem(Mesh& mesh) {
	for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
		if (isInterior(*v)) {
			const int idx = v->indexN;
			H.coeffRef(idx, idx) = 0;
			Bx[idx] = 0.0;
			By[idx] = 0.0;
		}
	}

	for (EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); ++e) {
		HalfEdgeIter he = e->halfEdge();
		HalfEdgeIter he2 = he->flip();
		VertexIter v_i = he->vertex();
		VertexIter v_j = he2->vertex();
		if (isInterior(*v_i) && isInterior(*v_j)) {
			H.coeffRef(v_i->indexN, v_j->indexN) = 0;
			H.coeffRef(v_j->indexN, v_i->indexN) = 0;
		}
	}
}

void SphereTutte::uniformWeight(Mesh& mesh) {
	std::vector<Triplet> H_coefficients;
	for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
		if (!isInterior(*v)) {
			continue;
		}
		const int idx_i = v->indexN;
		Bx[idx_i] = 0;
		By[idx_i] = 0;
		int deg_i = 0;
		const double weight = 1.0;
		HalfEdgeIter he = v->halfEdge();
		do {
			VertexIter v_j = he->flip()->vertex();
			if (isInterior(*v_j)) {
				H_coefficients.push_back(Triplet(idx_i, v_j->indexN, weight));
			} else if (v_j->onBoundaryTutte) {
				Bx[idx_i] -= weight * Xb[v_j->indexN];
				By[idx_i] -= weight * Yb[v_j->indexN];
			}
			he = he->flip()->next();
			++deg_i;
		} while (he != v->halfEdge());
		H_coefficients.push_back(Triplet(idx_i, idx_i, -deg_i));
	}

	H.setFromTriplets(H_coefficients.begin(), H_coefficients.end());
	H_is_set = true;
}

void SphereTutte::cotanWeight(Mesh& mesh) {
	if (H_is_set) {
		resetSystem(mesh);
	} else {
		initSystem(mesh);
	}

	const double floor_w = 0.01;
	for (EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); ++e) {
		HalfEdgeIter he = e->halfEdge();
		HalfEdgeIter he2 = he->flip();
		VertexIter v_i = he->vertex();
		VertexIter v_j = he2->vertex();
		if (isPole(*v_i) || isPole(*v_j)) {
			continue;
		}

		double weight = he->cotan() + he2->cotan();
		if (weight < 0) {
			weight = floor_w;
		}

		const int idx = v_i->indexN;
		const int jdx = v_j->indexN;
		if (isInterior(*v_i) && isInterior(*v_j)) {
			H.coeffRef(idx, jdx) = weight;
			H.coeffRef(jdx, idx) = weight;
			H.coeffRef(idx, idx) -= weight;
			H.coeffRef(jdx, jdx) -= weight;
		} else if (v_i->onBoundaryTutte && isInterior(*v_j)) {
			H.coeffRef(jdx, jdx) -= weight;
			Bx[jdx] -= weight * Xb[idx];
			By[jdx] -= weight * Yb[idx];
		} else if (isInterior(*v_i) && v_j->onBoundaryTutte) {
			H.coeffRef(idx, idx) -= weight;
			Bx[idx] -= weight * Xb[jdx];
			By[idx] -= weight * Yb[jdx];
		}
	}
}

void SphereTutte::meanValueWeight(Mesh& mesh) {
	if (H_is_set) {
		resetSystem(mesh);
	} else {
		initSystem(mesh);
	}

	const double floor_w = 0.01;
	for (EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); ++e) {
		HalfEdgeIter he = e->halfEdge();
		HalfEdgeIter he2 = he->flip();
		VertexIter v_i = he->vertex();
		VertexIter v_j = he2->vertex();
		if (isPole(*v_i) || isPole(*v_j)) {
			continue;
		}

		const Vector a = he->vertex()->position;
		const Vector b = he->next()->vertex()->position;
		const Vector c = he->prev()->vertex()->position;
		const Vector d = he2->prev()->vertex()->position;

		const double l_ij = (a - b).norm();
		const double l_jk = (b - c).norm();
		const double l_ik = (a - c).norm();
		const double l_jl = (b - d).norm();
		const double l_il = (a - d).norm();
		if (l_ij <= 0.0) {
			continue;
		}

		double s = 0.5 * (l_jk + l_ij + l_ik);
		double N = (s - l_ij) * (s - l_ik);
		double D = s * (s - l_jk);
		const double tan_alpha_1 = (N > 0.0 && D > 0.0) ? std::sqrt(N / D) : 0.0;

		s = 0.5 * (l_jl + l_ij + l_il);
		N = (s - l_ij) * (s - l_il);
		D = s * (s - l_jl);
		const double tan_alpha_2 = (N > 0.0 && D > 0.0) ? std::sqrt(N / D) : 0.0;

		double weight = (tan_alpha_1 + tan_alpha_2) / l_ij;
		if (weight < 0) {
			weight = floor_w;
		}

		const int idx = v_i->indexN;
		const int jdx = v_j->indexN;
		if (isInterior(*v_i) && isInterior(*v_j)) {
			H.coeffRef(idx, jdx) = weight;
			H.coeffRef(jdx, idx) = weight;
			H.coeffRef(idx, idx) -= weight;
			H.coeffRef(jdx, jdx) -= weight;
		} else if (v_i->onBoundaryTutte && isInterior(*v_j)) {
			H.coeffRef(jdx, jdx) -= weight;
			Bx[jdx] -= weight * Xb[idx];
			By[jdx] -= weight * Yb[idx];
		} else if (isInterior(*v_i) && v_j->onBoundaryTutte) {
			H.coeffRef(idx, idx) -= weight;
			Bx[idx] -= weight * Xb[jdx];
			By[idx] -= weight * Yb[jdx];
		}
	}
}

void SphereTutte::authalicWeight(Mesh& mesh) {
	if (H_is_set) {
		resetSystem(mesh);
	} else {
		initSystem(mesh);
	}

	const double floor_w = 0.01;
	for (EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); ++e) {
		HalfEdgeIter he = e->halfEdge();
		HalfEdgeIter he2 = he->flip();
		VertexIter v_i = he->vertex();
		VertexIter v_j = he2->vertex();
		if (isPole(*v_i) || isPole(*v_j)) {
			continue;
		}

		const Vector a = he->vertex()->position;
		const Vector b = he->next()->vertex()->position;
		const Vector c = he->prev()->vertex()->position;
		const Vector d = he2->prev()->vertex()->position;

		const double l_ij = (a - b).norm();
		const double l_jk = (b - c).norm();
		const double l_ik = (a - c).norm();
		const double l_jl = (b - d).norm();
		const double l_il = (a - d).norm();
		if (l_ij <= 0.0) {
			continue;
		}

		double N = l_ij * l_ij + l_jk * l_jk - l_ik * l_ik;
		double D = 2 * l_ij * l_jk;
		double cos_gamma_1 = (D == 0.0) ? 0.0 : N / D;
		cos_gamma_1 = std::max(-1.0, std::min(1.0, cos_gamma_1));
		double sin_gamma_1 = std::sqrt(std::max(0.0, 1.0 - cos_gamma_1 * cos_gamma_1));
		const double cot_gamma_1 = (sin_gamma_1 > 0.0) ? cos_gamma_1 / sin_gamma_1 : 0.0;

		N = l_ij * l_ij + l_jl * l_jl - l_il * l_il;
		D = 2 * l_ij * l_jl;
		double cos_gamma_2 = (D == 0.0) ? 0.0 : N / D;
		cos_gamma_2 = std::max(-1.0, std::min(1.0, cos_gamma_2));
		double sin_gamma_2 = std::sqrt(std::max(0.0, 1.0 - cos_gamma_2 * cos_gamma_2));
		const double cot_gamma_2 = (sin_gamma_2 > 0.0) ? cos_gamma_2 / sin_gamma_2 : 0.0;

		double weight = (cot_gamma_1 + cot_gamma_2) / (l_ij * l_ij);
		if (weight < 0) {
			weight = floor_w;
		}

		const int idx = v_i->indexN;
		const int jdx = v_j->indexN;
		if (isInterior(*v_i) && isInterior(*v_j)) {
			H.coeffRef(idx, jdx) = weight;
			H.coeffRef(jdx, idx) = weight;
			H.coeffRef(idx, idx) -= weight;
			H.coeffRef(jdx, jdx) -= weight;
		} else if (v_i->onBoundaryTutte && isInterior(*v_j)) {
			H.coeffRef(jdx, jdx) -= weight;
			Bx[jdx] -= weight * Xb[idx];
			By[jdx] -= weight * Yb[idx];
		} else if (isInterior(*v_i) && v_j->onBoundaryTutte) {
			H.coeffRef(idx, idx) -= weight;
			Bx[idx] -= weight * Xb[jdx];
			By[idx] -= weight * Yb[jdx];
		}
	}
}

void SphereTutte::setBoundary(double radius) {
	Xb.resize(Nb);
	Yb.resize(Nb);
	const double angle = 2.0 * kPi / static_cast<double>(Nb);
	for (int i = 0; i < Nb; ++i) {
		Xb[i] = radius * std::cos(i * angle);
		Yb[i] = radius * std::sin(i * angle);
	}
}

SphereEmbedding SphereTutte::embed(Mesh& mesh, int pole_type, int Tutte_type, bool verbose) {
	H_is_set = false;
	poleIndex = -1;
	Nb = 0;

	VertexIter pole = pickPole(mesh, pole_type);
	poleIndex = pole->index;
	punctureMesh(mesh, pole);

	const double radius = 100.0;
	setBoundary(radius);

	const int n_vertices = static_cast<int>(mesh.vertices.size());
	const int n_active = n_vertices - 1 - Nb;
	if (n_active < 0) {
		restorePunctureFlags(mesh);
		throw std::runtime_error("Sphere Tutte: negative number of interior vertices");
	}

	H.resize(n_active, n_active);
	Bx.resize(n_active);
	By.resize(n_active);
	X.resize(n_active);
	Y.resize(n_active);
	H.setZero();
	if (n_active > 0) {
		Bx.setZero();
		By.setZero();
	}

	std::string type;
	if (Tutte_type == 1) {
		uniformWeight(mesh);
		type = "Uniform";
	} else if (Tutte_type == 2) {
		cotanWeight(mesh);
		type = "coTan";
	} else if (Tutte_type == 3) {
		meanValueWeight(mesh);
		type = "MeanValue";
	} else if (Tutte_type == 4) {
		authalicWeight(mesh);
		type = "Authalic";
	} else {
		restorePunctureFlags(mesh);
		throw std::runtime_error("Invalid Tutte embedding type number: " + std::to_string(Tutte_type));
	}

	double err1 = 0.0;
	double err2 = 0.0;
	if (n_active > 0) {
		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
		solver.analyzePattern(H);
		solver.factorize(H);
		if (solver.info() != Eigen::Success) {
			restorePunctureFlags(mesh);
			throw std::runtime_error("Failed LU factorization for sphere Tutte linear system");
		}
		X = solver.solve(Bx);
		err1 = (H * X - Bx).norm();
		if (solver.info() != Eigen::Success) {
			restorePunctureFlags(mesh);
			throw std::runtime_error("Failed solving for X coordinates in sphere Tutte");
		}
		Y = solver.solve(By);
		err2 = (H * Y - By).norm();
		if (solver.info() != Eigen::Success) {
			restorePunctureFlags(mesh);
			throw std::runtime_error("Failed solving for Y coordinates in sphere Tutte");
		}
	}

	if (verbose) {
		std::cout << "Sphere Tutte embedding:\n";
		std::cout << "  pole vertex ID       : " << poleIndex << "\n";
		std::cout << "  boundary points      : " << Nb << "\n";
		std::cout << "  interior points      : " << n_active << "\n";
		std::cout << "  edge weights         : " << type << "\n";
		std::cout << "  error X / Y          : " << err1 << " / " << err2 << "\n";
	}

	if (err1 > 1.e-5 || err2 > 1.e-5) {
		std::cout << "Warning: large residual when solving sphere Tutte embedding\n";
	}

	const double ang1 = kPi / 180.0;
	double factor = std::sin(ang1) / (1.0 - std::cos(ang1));
	factor = factor / radius;

	Eigen::MatrixXd sphere(n_vertices, 3);
	sphere.setZero();
	for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v) {
		const int row = v->index;
		double U = 0.0;
		double V = 0.0;
		if (isPole(*v)) {
			sphere(row, 0) = 0.0;
			sphere(row, 1) = 0.0;
			sphere(row, 2) = -1.0;
			continue;
		}
		if (isInterior(*v)) {
			U = factor * X[v->indexN];
			V = factor * Y[v->indexN];
		} else {
			U = factor * Xb[v->indexN];
			V = factor * Yb[v->indexN];
		}
		const double val = U * U + V * V;
		const double den = 1.0 / (1.0 + val);
		sphere(row, 0) = 2.0 * U * den;
		sphere(row, 1) = 2.0 * V * den;
		sphere(row, 2) = (1.0 - val) * den;
	}

	restorePunctureFlags(mesh);

	return SphereEmbedding{std::move(sphere), poleIndex};
}

}  // namespace geomp
