/* ====SPHERE TUTTE =================================================================================
 *
 * Genus-zero Tutte embedding: puncture a pole, embed the rest in the plane, then inverse
 * stereographic projection onto the unit sphere.
 *
 * Rewritten from include/geomp/embedding/sphere/STutte.h against the current Mesh API.
 * Class name SphereTutte avoids colliding with geomp::Tutte (torus).
 *
 * Author: Patrice Koehl (original), Yuan Luo
 =============================================================================================== */

#pragma once

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <cmath>
#include <utility>
#include <vector>

#include <Eigen/Sparse>

#include "geomp/mesh/Mesh.h"
#include "geomp/mesh/Types.h"

namespace geomp {

typedef Eigen::Triplet<double> Triplet;

struct SphereEmbedding {
	Eigen::MatrixXd vertices;  // (V, 3) unit-sphere coordinates, rows indexed by vertex->index
	int poleIndex = -1;
};

class SphereTutte {
public:
	// pole_type: 0 = first vertex, 1 = highest valence, 2 = lowest valence
	// Tutte_type: 1 Uniform, 2 coTan, 3 MeanValue, 4 Authalic
	SphereEmbedding embed(Mesh& mesh, int pole_type = 0, int Tutte_type = 2, bool verbose = false);

private:
	VertexIter pickPole(Mesh& mesh, int type) const;
	void punctureMesh(Mesh& mesh, VertexIter pole);
	void restorePunctureFlags(Mesh& mesh);

	bool isInterior(const Vertex& v) const;
	bool isPole(const Vertex& v) const;

	void initSystem(Mesh& mesh);
	void resetSystem(Mesh& mesh);
	void uniformWeight(Mesh& mesh);
	void cotanWeight(Mesh& mesh);
	void meanValueWeight(Mesh& mesh);
	void authalicWeight(Mesh& mesh);
	void setBoundary(double radius);

	int poleIndex = -1;
	int Nb = 0;
	Eigen::SparseMatrix<double> H;
	Eigen::VectorXd X, Y, Xb, Yb, Bx, By;
	bool H_is_set = false;
};

}  // namespace geomp
