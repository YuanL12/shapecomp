/* ====TUTTE.H ====================================================================================
 *
 * Author: Yuan Luo, Nov 2024, based on the work of Patrice Koehl
 * Department of Mathematics, University of California, Davis
 *
 * This file implements different Tutte embeddings of a torus in the plane.
 * The boundary points are evenly spaced on the unit squre.
 * 
 *
 * A Tutte embedding in the plane satisfies for each vertex i:
 *              Sum_j (w'_ij xj) - xi = 0				(1)
 *              Sum_j (w'_ij yj) - yi = 0				(2)
 * where the sum extends to all vertices j in the star of i, and w'_ij are normalized weights, i.e.
 *              w'_ij = w_ij / (Sum_k w_ik)				(3)
 * where the original weights w_ij defines the type of Tutte embedding. We have included:
 *
 *		1. Embedding with edge weights w_ij set to 1; in this case the Tutte
 *		   weights w'_ij = 1/N(i), valence of i, are not symmetric
 *		2. Embedding with edge weights w_ij = l_ij ; in this case again
 *		   the Tutte weights w'_ij = l_ij / (Sum_k l_ik) lare not symmetric
 *		3. Embedding with cotan weights as w_ij = cot(beta_k) + cot(beta_l) (co-tangent weights)
 *		   the Tutte weights are symmetric, but may not be positive
 *		   This embedding is harmonic
 *		4. Embedding with tan weights as w_ij = (tan(a_i/2) + tan(a'_i/2))/l_ij (half-tangent weights)
 *		   the Tutte weights are positive 
 *		   This is the mean value embedding
 *		5. Embedding with modified cotangent weights
 *		   This is an authalic (area preserving) embedding
 *
 * Linear System view is:
 * Let V1 be the set of interior vertices and V2 the set of boundary vertices
 * Denote by L[V1, V1] the submatirx of Laplacian matrix L of the interior parts
 * Denote by L[V1, V2] the submatirx of Laplacian matrix L of the interior (as rows) and boundary parts (as columns)
 * where L[i,i] = -deg(v_i) and L[i,j] = 1 if v_i and v_j are connected by an edge.
 * We want to solve L[V1, V1] \phi(V1) + L[V1, V2] \phi(V2) = 0.   
 * Algorithmically, L[V1, V1] \phi(V1) = H \phi(V1)  = - L[V1, V2] \phi(V2) = B, 
 * where B is the constraints of the system, H will be a sparse square matrix.
 =============================================================================================== */


#pragma once
#include <vector>
#include <cmath>
#include <Eigen/Sparse>
#include "geomp/mesh/MeshData.h"
#include "geomp/mesh/Types.h"
#include "geomp/generators/Generators.h"
#include "extern/icecream/icecream.hpp"

namespace geomp {

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
typedef Eigen::Triplet<double> Triplet;

class Tutte{
// Careful: the number of interior vertices is 
// int n_active = n_vertices - 1 - Nb; for a sphere (need to remove one pole)
// int n_active = n_vertices - Nb; for a torus

public:
	// Performs planar Tutte embedding
	VertexData<std::pair<double, double>> planarTutte(Mesh& mesh, CircularLinkedList& nodesOnGenerators, 
							int Tutte_type = 1, bool verbose = false);

private:

	// assign indexN to vertices to differentiate between boundary and interior vertices
	void indexVertices(Mesh& mesh, CircularLinkedList& nodesOnGenerators);

	void setSquareBoundary(CircularLinkedList& nodesOnGenerators);

	// init matrix for Tutte linear system
	void initSystem(Mesh& mesh);

	// Reset matrix for Tutte linear system to 0
	void resetSystem(Mesh& mesh);

 	// Set uniform, non symmetric weights for the Tutte embedding based on valence
	void uniformWeight(Mesh& mesh);

	// Set cotangent, symmetric weights for the Tutte embedding
	void cotanWeight(Mesh& mesh);

	// Set half-tangent weights for the Tutte embedding
	void meanValueWeight(Mesh& mesh);

	// Set modified co-tangent weights for the Tutte embedding
	void authalicWeight(Mesh& mesh);

    protected:

	// number of points on the boundary
	int Nb;

	// The indices system for the Tutte embedding is constructed by calling indexVertices()
	// e.g., X[v->indexN] is the x-coordinate of an interior vertex v
	// Xb[v->indexN] is the x-coordinate of a boundary vertex v
	// Laplacian matrix of the interior parts 
	Eigen::SparseMatrix<double> H;

	// X, Y are the x,y coordinates of the interior vertices to be solved
	Eigen::VectorXd X, Y;

	// X, Y coordiantes on the boundary
	Eigen::VectorXd Xb, Yb;

	// Bx, By are the constraints of the system for solving X and Y respectively
	Eigen::VectorXd Bx, By;

	// flag to indicate if matrix H is set
	bool H_is_set = false;

};


}// end of namespace geomp