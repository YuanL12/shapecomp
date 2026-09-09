/* ====TUTTE.H ====================================================================================
 *
 * Author: Patrice Koehl, June 2018
 * Department of Computer Science
 * University of California, Davis
 *
 * This file implements different Tutte embeddings in the plane.
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
 *
 =============================================================================================== */

#pragma once

/* ===== INCLUDES AND TYPEDEFS =================================================================
*
* Third party libraries included: Eigen
*
=============================================================================================== */

#define _USE_MATH_DEFINES // for M_PI

#pragma once
#include <vector>
#include <cmath>
// #include <Eigen/Core>
#include <Eigen/Sparse>
// #include <Eigen/SparseLU>
#include "Types.h"
#include "Mesh.h"

using Triplet = Eigen::Triplet<double>;


  /* =============================================================================================
     Define class
   =============================================================================================== */

class Tutte{

public:
	// select vertex from genus 0 surface that will be remove temperarily
	VertexIter pickPole(Mesh& mesh, int type);

	// Removes a vertex from the mesh and set the tags.
	void punctureMesh(Mesh& mesh, VertexIter pole);
	
	// Performs planar Tutte embedding
	VertexIter planarTutte(Mesh& mesh, int vtype, int Tutte_type, bool *SUCCESS, bool verbose);

private:
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

	// Set position of boundary vertices in the plane
	void setBoundary(double Radius);

	// Projection to sphere
	void stereo2Sphere(Mesh& mesh);

    protected:

	Eigen::SparseMatrix<double> H;
	Eigen::VectorXd Bx, By, X, Y;

	// number of points on the boundary
	int Nb;
	// X, Y coordiantes on the boundary
	std::vector<double> Xb;
	std::vector<double> Yb;

	bool H_is_set = false;

  };

  /* ===== PICK VERTEX    ========================================================================
   * Input:
   * 	@Mesh -  Mesh in OpenMesh format
   *
   * Select vertex that will be removed from the mesh to enable planar Tutte embedding.
   * At this stage, selection is based on finding the vertex with either the highest valence
   * (type = 1), or the lowest valence (type = 2). If type is set to 0, vertex # 0 is picked
   * arbitrarily
   *
   =============================================================================================== */

  VertexIter Tutte::pickPole(Mesh& mesh, int type)
  {
	if(type==0) return mesh.vertices.begin();

	/* Find vertex with largest, or smallest valence */

	int valence;
	int val_ref = mesh.vertices.begin()->degree();
	VertexIter pole=mesh.vertices.begin();

	for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
	{
		valence = v->degree();

		if(type==2) {
			if(valence < val_ref) {
				val_ref = valence;
				pole = v;
			}
		} else if(type==1) {
			if(valence > val_ref) {
				val_ref = valence;
				pole = v;
			}
		}
	}

	return pole;
  }

  /* ===== PUNCTURE MESH  ========================================================================
   * Input:
   * 	@Mesh -  Mesh
   *    VR    -  Vertex to be removed
   *
   * Remove vertex VR.  All incident
   * vertices and faces  are set to inNorthPoleVicinity
   *
   =============================================================================================== */
  	void Tutte::punctureMesh(Mesh& mesh, VertexIter pole){
		for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++){
			v->NorthPole = false;
			v->inNorthPoleVicinity = false;
		}
		// loop over all faces/vertices(star) surrouding the pole
		// to index all points on the boundary
		int idx = 0;
		pole->NorthPole = true;
		HalfEdgeIter he = pole->he;
		do {
			he->face->inNorthPoleVicinity = true;
			he->flip->vertex->indexN = idx;
			he->flip->vertex->inNorthPoleVicinity = true;
			idx++;
			he = he->flip->next;
		} while (he != pole->he); //go back to the starting edge 

		// store the number of points on the boundary
		Nb = idx;

		// index the interior points 
		idx = 0; 
		for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++){
			bool b1 = v->NorthPole;
			bool b2 = v->inNorthPoleVicinity;
			// neither the pole or on the boundary
			if(!b1 && !b2) { 
				v->indexN = idx;
				idx++;
			}
		}

  	}// end punctureMesh

  /* =============================================================================================
   Init matrix
   =============================================================================================== */

  void Tutte::initSystem(Mesh& mesh)
  {

  /* ==================================================================================================
	Initialize all values for regular constraints to 0
   ==================================================================================================*/

	int idx, jdx;
	std::vector<Triplet> Mat_coefficients;
	double zero = 0;
	bool b1, b2, b3, b4;

	for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
	{
		b1 = v->NorthPole;
		b2 = v->inNorthPoleVicinity;
		if(!b1 && !b2) {
			idx = v->indexN;
			Mat_coefficients.push_back(Triplet(idx, idx, zero));
			Bx[idx] = 0.0;
			By[idx] = 0.0;
		}
	}

	HalfEdgeIter he, he2;
	VertexIter v_i, v_j;
	
	for (EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); e++)
	{
		he = e->he;
		he2 = he->flip;

		v_i = he->vertex;
		v_j = he2->vertex;

		b1 = v_i->NorthPole;
		b2 = v_i->inNorthPoleVicinity;
		b3 = v_j->NorthPole;
		b4 = v_j->inNorthPoleVicinity;

		if(!b1 && !b2 && !b3 && !b4) {
			idx = v_i->indexN;
			jdx = v_j->indexN;
			Mat_coefficients.push_back(Triplet(idx, jdx, zero));
			Mat_coefficients.push_back(Triplet(jdx, idx, zero));
		}
	}

	H.setFromTriplets(Mat_coefficients.begin(), Mat_coefficients.end());
	H_is_set = true;

  }

  /* =============================================================================================
   Reset matrix to zero
   =============================================================================================== */

  void Tutte::resetSystem(Mesh& mesh)
  {

	int idx, jdx;
	bool b1, b2, b3, b4;

	for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
	{
		b1 = v->NorthPole;
		b2 = v->inNorthPoleVicinity;
		if(!b1 && !b2) {
			idx = v->indexN;
			H.coeffRef(idx, idx)     = 0;
			Bx[idx] = 0.0;
			By[idx] = 0.0;
		}
	}

	HalfEdgeIter he, he2;
	VertexIter v_i, v_j;
	
	for (EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); e++)
	{
		he = e->he;
		he2 = he->flip;

		v_i = he->vertex;
		v_j = he2->vertex;

		b1 = v_i->NorthPole;
		b2 = v_i->inNorthPoleVicinity;
		b3 = v_j->NorthPole;
		b4 = v_j->inNorthPoleVicinity;

		if(!b1 && !b2 && !b3 && !b4) {
			idx = v_i->indexN;
			jdx = v_j->indexN;
			H.coeffRef(idx, jdx)     = 0;
			H.coeffRef(idx, jdx)     = 0;
		}
	}
  }

  /* ===== Uniform vertex weights ====================================================================
   * Computed as:
   * w_ij = 1/N(i)
   * where N(i) is the valence of vertex i.
   * Note: not symmetric; barycentric embedding
   *
   * Initial embedding from Tutte: 
	W.T. Tutte (1963). How to draw a graph. Proc. Lond. Math. Soc. 13, 743-767
   ==================================================================================================*/

  void Tutte::uniformWeight(Mesh& mesh)
  {
	std::vector<Triplet> H_coefficients;
	int idx, jdx;
	double weight;
	VertexIter v_j;
	for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
	{
		bool b1 = v->NorthPole;
		bool b2 = v->inNorthPoleVicinity;

		if(!b1 && !b2) {

			idx = v->indexN;
			Bx[idx] = 0;
			By[idx] = 0;
			weight = -v->degree();
			H_coefficients.push_back(Triplet(idx, idx, weight));

			weight = 1.0;
			HalfEdgeIter he = v->he;
			do {
				v_j = he->flip->vertex;
				bool b3 = v_j->NorthPole;
				bool b4 = v_j->inNorthPoleVicinity;
				jdx = v_j->indexN;

				if(!b3 && !b4) { // interior vertex
					H_coefficients.push_back(Triplet(idx, jdx, weight));
				} else if (b4) { // boundary vertex
					Bx[idx] = Bx[idx] - weight*Xb[jdx];
					By[idx] = By[idx] - weight*Yb[jdx];
				}
				he= he->flip->next;
			} while (he != v->he);
				
		}
	}

	H.setFromTriplets(H_coefficients.begin(), H_coefficients.end());
	H_is_set = true;

  }

  /* ===== Cotan-based edge weights =================================================================
    Computed by considering the two triangles incident to (ij): (ijk) and (ijl)
    w_ij = (cotan(beta_k) + cotan(beta_l))
    This is an example of harmonic map; note that w_ij are symmetric, but may be negative
   
    From U. Pinkall and K. Polthier (1993). Computing Discrete Minimal Surfaces and Their Conjugates.
	  Experiment. Math., 2, 15-36  
   ==================================================================================================*/

  void Tutte::cotanWeight(Mesh& mesh)
  {

	if(H_is_set) {
		resetSystem(mesh);
	} else {
		initSystem(mesh);
		H_is_set = true;
	}

	double weight, cot_beta_1, cot_beta_2;
	double floor = 0.01;

        HalfEdgeIter he, he2;
        VertexIter v_i, v_j;
	bool b1, b2, b3, b4;
	int idx, jdx;

	for (EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); e++)
	{
		he = e->he;
		he2 = he->flip;

		v_i = he->vertex;
		v_j = he2->vertex;
		idx = v_i->index;
		jdx = v_j->index;

		b1 = v_i->NorthPole;
		b2 = v_i->inNorthPoleVicinity;
		b3 = v_j->NorthPole;
		b4 = v_j->inNorthPoleVicinity;

		if(!b1 && !b3) {

			cot_beta_1 = he->cotan();
			cot_beta_2 = he2->cotan();
			weight = cot_beta_1 + cot_beta_2;
			if(weight < 0) weight = floor;

                       	idx = v_i->indexN;
                       	jdx = v_j->indexN;
			if( !b2 && !b4) {
				H.coeffRef(idx, jdx) = weight;
				H.coeffRef(jdx, idx) = weight;
				H.coeffRef(idx, idx) -= weight;
				H.coeffRef(jdx, jdx) -= weight;
			} else {
				if(b2 && !b4) {
					H.coeffRef(jdx, jdx) -= weight;
					Bx[jdx] -= weight*Xb[idx];
					By[jdx] -= weight*Yb[idx];
				} else if(!b2 && b4) {
					H.coeffRef(idx, idx) -= weight;
					Bx[idx] -= weight*Xb[jdx];
					By[idx] -= weight*Yb[jdx];
				}
			}

		}
	}
  }

  /* ===== Mean value weights ========================================================================
   * Based on: Floater, M. 1997. Parametrization and smooth approximation of surface triangulations. 
		CAGD 14, 3, 231–250
   ==================================================================================================*/

  void Tutte::meanValueWeight(Mesh& mesh)
  {

	if(H_is_set) {
		resetSystem(mesh);
	} else {
		initSystem(mesh);
		H_is_set = true;
	}

	int idx, jdx;
	double l_ij, l_ik, l_jk;
	double l_il, l_jl;
	double tan_alpha_1, tan_alpha_2;
	double weight, s, N, D;
	double floor = 0.01;

        HalfEdgeIter he, he2;
        VertexIter v_i, v_j;
	bool b1, b2, b3, b4;
	Vector a, b, c, d;

	for (EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); e++)
	{
		he = e->he;
		he2 = he->flip;

		v_i = he->vertex;
		v_j = he2->vertex;

		b1 = v_i->NorthPole;
		b2 = v_i->inNorthPoleVicinity;
		b3 = v_j->NorthPole;
		b4 = v_j->inNorthPoleVicinity;

		if(!b1 && !b3) {

			a = he->vertex->position;
			b = he->next->vertex->position;
			c = he->prev->vertex->position;
			d = he2->prev->vertex->position;

			l_ij = (a-b).norm();
			l_jk = (b-c).norm();
			l_ik = (a-c).norm();
			l_jl = (b-d).norm();
			l_il = (a-d).norm();

			/* Calculate alpha_1 */
			s = 0.5*(l_jk + l_ij + l_ik);
			N = (s-l_ij)*(s-l_ik);
			D = s*(s-l_jk);
			tan_alpha_1 = std::sqrt(N/D);

			/* Calculate alpha_2 */
			s = 0.5*(l_jl + l_ij + l_il);
			N = (s-l_ij)*(s-l_il);
			D = s*(s-l_jl);
			tan_alpha_2 = std::sqrt(N/D);

			/* Weight as calculated */
			weight = (tan_alpha_1 + tan_alpha_2)/l_ij;
			if(weight < 0) weight = floor;

                        idx = v_i->indexN;
                        jdx = v_j->indexN;
			if( !b2 && !b4) {
				H.coeffRef(idx, jdx) = weight;
				H.coeffRef(jdx, idx) = weight;
				H.coeffRef(idx, idx) -= weight;
				H.coeffRef(jdx, jdx) -= weight;
			} else {
				if(b2 && !b4) {
					H.coeffRef(jdx, jdx) -= weight;
					Bx[jdx] -= weight*Xb[idx];
					By[jdx] -= weight*Yb[idx];
				} else if(!b2 && b4) {
					H.coeffRef(idx, idx) -= weight;
					Bx[idx] -= weight*Xb[jdx];
					By[idx] -= weight*Yb[jdx];
				}
			}
		}
	}
  }

  /* ===== Authalic edge weights =================================================================
   * Computed by considering the two triangles incident to (ij): (ijk) and (ijl)
   * w_ij = (cotan(gamma_j) + cotan(gamma'_j))/l_ij^2
   * This is an example of an authalic map (area preserving map); 
   *	note that w_ij are symmetric, but may be negative
   *
   * From  M. Desbrun, M. Meyer, and P. Alliez. (2002). Intrinsic Parametrizations of Surface Meshes
   * 	   Eurographics, 21, 209-218.
   ==================================================================================================*/

  void Tutte::authalicWeight(Mesh& mesh)
  {

	if(H_is_set) {
		resetSystem(mesh);
	} else {
		initSystem(mesh);
		H_is_set = true;
	}

	int idx, jdx;
	double l_ij, l_ik, l_jk;
	double l_il, l_jl;
	double cos_gamma_1, cos_gamma_2;
	double sin_gamma_1, sin_gamma_2;
	double cot_gamma_1, cot_gamma_2;
	double weight, N, D;
	double floor = 0.01;

        HalfEdgeIter he, he2;
        VertexIter v_i, v_j;
	bool b1, b2, b3, b4;
	Vector a, b, c, d;

	for (EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); e++)
	{
		he = e->he;
		he2 = he->flip;

		v_i = he->vertex;
		v_j = he2->vertex;

		b1 = v_i->NorthPole;
		b2 = v_i->inNorthPoleVicinity;
		b3 = v_j->NorthPole;
		b4 = v_j->inNorthPoleVicinity;

		if(!b1 && !b3) {

			a = he->vertex->position;
			b = he->next->vertex->position;
			c = he->prev->vertex->position;
			d = he2->prev->vertex->position;

			l_ij = (a-b).norm();
			l_jk = (b-c).norm();
			l_ik = (a-c).norm();
			l_jl = (b-d).norm();
			l_il = (a-d).norm();

			/* Calculate gamma_1 */
			N = l_ij*l_ij + l_jk*l_jk - l_ik*l_ik;
			D = 2*l_ij*l_jk;
			cos_gamma_1 = N/D;
			sin_gamma_1 = std::sqrt(1-cos_gamma_1*cos_gamma_1);
			cot_gamma_1 = cos_gamma_1/sin_gamma_1;

			/* Calculate gamma_2 */
			N = l_ij*l_ij + l_jl*l_jl - l_il*l_il;
			D = 2*l_ij*l_jl;
			cos_gamma_2 = N/D;
			sin_gamma_2 = std::sqrt(1-cos_gamma_2*cos_gamma_2);
			cot_gamma_2 = cos_gamma_2/sin_gamma_2;

			/* Weight as calculated */
			weight = (cot_gamma_1 + cot_gamma_2)/(l_ij*l_ij);
			if(weight < 0) weight = floor;

                        idx = v_i->indexN;
                        jdx = v_j->indexN;
			if( !b2 && !b4) {
				H.coeffRef(idx, jdx) = weight;
				H.coeffRef(jdx, idx) = weight;
				H.coeffRef(idx, idx) -= weight;
				H.coeffRef(jdx, jdx) -= weight;
			} else {
				if(b2 && !b4) {
					H.coeffRef(jdx, jdx) -= weight;
					Bx[jdx] -= weight*Xb[idx];
					By[jdx] -= weight*Yb[idx];
				} else if(!b2 && b4) {
					H.coeffRef(idx, idx) -= weight;
					Bx[idx] -= weight*Xb[jdx];
					By[idx] -= weight*Yb[jdx];
				}
			}
		}
	}
  }
  /* ===== SET BOUNDARY   ============================================================================
   * Set coordinates of boundary points on the plane
   ==================================================================================================*/

  void Tutte::setBoundary(double Radius){
	double angle = 2.0*M_PI/Nb;
	for(int i = 0; i < Nb; i++){
		Xb.push_back(Radius*std::cos(i*angle));
		Yb.push_back(Radius*std::sin(i*angle));
	}
  }

  	/*
	Project planar mesh onto unit sphere, 
	Using inverse of South pole stereographic projection, the plane is through the equator
	*/
	void Tutte::stereo2Sphere(Mesh& mesh){
		// (u,v) -> (2u/(u^2+v^2+1), , )
		double U, V, val, den;
		for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++){
			bool b1 = v->NorthPole;
			if(b1) {
				v->position2.x = 0.;
				v->position2.y = 0.;
				v->position2.z = -1.;
			} else {
				U = v->position2.x;
				V = v->position2.y;
				val = U*U + V*V;
				den = 1.0/(1.0+val);
				v->position2.x = 2.0*U*den;
				v->position2.y = 2.0*V*den;
				v->position2.z = (1.0 - val)*den;
			}
		}
  	}

  /* ===== Tutte Embedding ============================================================================
   * Solve for vertex coordinates on the plane using Tutte Embedding
   ==================================================================================================*/

  VertexIter Tutte::planarTutte(Mesh& mesh, int vtype, int Tutte_type, bool *SUCCESS, bool verbose = false)
  {

	*SUCCESS = true;

  	// Picking a pole and label its boundary vertices
	VertexIter pole = pickPole(mesh, vtype);
	punctureMesh(mesh, pole);

  	// Set boundary on the plane based on number of points that have been removed
	double Radius = 100.0;
	setBoundary(Radius);

	/*
	====================================================================================
	Define Tutte Embedding System
	====================================================================================
	*/
	// Set matrices for solving system of linear equations 
	int n_vertices = mesh.vertices.size();
	int n_active = n_vertices - 1 - Nb;
	H.resize(n_active, n_active);
	Bx.resize(n_active); // constraints for x coordinates of interior vertices
	By.resize(n_active);
	X.resize(n_active);
	Y.resize(n_active);
	H.setZero();

  	// Based on type of Tutte embedding, define Vertex Weights and Edge Weights
	std::string type;
	if(Tutte_type == 1) {
		uniformWeight(mesh);
		type="Uniform";
	} else if(Tutte_type == 2) {
		cotanWeight(mesh);
		type="coTan";
	} else if(Tutte_type == 3) {
		meanValueWeight(mesh);
		type="MeanValue";
	} else if(Tutte_type == 4) {
		authalicWeight(mesh);
		type="Authalic";
	}

  	// Solve linear systems using Sparse LU from Eigen
	Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
	// Eigen::UmfPackLU<Eigen::SparseMatrix<double>> solver; // needs suitesparse
	solver.analyzePattern(H);
	solver.factorize(H);
	if (solver.info() != Eigen::Success){
		std::cout << "Failed LU factorization for Tutte linear system" << std::endl;
		exit(1);
	}
	// x coordinates on plane
	X = solver.solve(Bx); 
	double err1 = (H*X - Bx).norm();
	if (solver.info() != Eigen::Success){
		std::cout << "Failed Solving for X coordinates" << std::endl;
		exit(1);
	}
	// y coordinates on plane
	Y = solver.solve(By);
	double err2 = (H*Y - By).norm();
	if (solver.info() != Eigen::Success){
		throw std::runtime_error("Tutte Embedding failed because unable to solving for Y coordinates");
	}

	// print results
	if (verbose){
		std::cout << " " << std::endl;
		std::cout << "Tutte embedding: " << std::endl;
		std::cout << "=================" << std::endl;
		std::cout << " " << std::endl;
		std::cout << "Mesh punctured at vertex ID                          : " << pole->index << std::endl;
		std::cout << "Number of boundary points                            : " << Nb << std::endl;
		std::cout << "Edge weight for embedding                            : " << type << std::endl;
		std::cout << "Error on X coordinates                               :" << err1 << std::endl;
		std::cout << "Error on Y coordinates                               :" << err2 << std::endl;
		std::cout << " " << std::endl;	
	}

	if(err1 > 1.e-5 || err2 > 1.e-5) {
		*SUCCESS = false;
		return pole;
	}

  	/* 
	====================================================================================
	Update coordinates of vertices: they should be in the plane now
	Scale so that the initial circle of radius Radius has radius
        sin(PI)/(1-cos(PI))
   	====================================================================================
	*/
	double ang1 = M_PI/180.0;
	double factor = std::sin(ang1)/(1.0-std::cos(ang1));
	factor = factor/Radius;
	int idx;
    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++){
		bool b1 = v->NorthPole;
		bool b2 = v->inNorthPoleVicinity;
		idx = v->indexN;
		if(!b1 && !b2) {
			v->position2.x = factor*X[idx];
			v->position2.y = factor*Y[idx];
			v->position2.z = 0;
		} else if(b2) {
			v->position2.x = factor*Xb[idx];
			v->position2.y = factor*Yb[idx];
			v->position2.z = 0;
		}
	}

    // Project onto sphere
	stereo2Sphere(mesh);

    // restore vertex star
	pole->NorthPole = false;
	HalfEdgeIter he = pole->he;
	do {
		he->face->inNorthPoleVicinity = false;
		he = he->flip->next;
	} while (he != pole->he);

	return pole;
  }
