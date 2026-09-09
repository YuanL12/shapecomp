#include "geomp/embedding/torus/Tutte.h"

namespace geomp{

void Tutte::indexVertices(Mesh& mesh, CircularLinkedList& nodesOnGenerators){
	resetHeadForTorus(nodesOnGenerators);
	// store the index mapping: index -> indexN
	std::vector<int> indexMapping(mesh.vertices.size());
	
	// give index to vertices on the boundary first
	int idx = 0;
	VertexNode* nd = nodesOnGenerators.getHead(); 
	do{
		VertexIter v = nd->vertexPtr;
		v->indexN = idx;
		v->onBoundaryTutte = true;
		indexMapping[v->index] = v->indexN;
		idx++;
		nd = nd->next;
	}while(nd != nodesOnGenerators.getHead());

	// index now is the number of boundary vertices
	Nb = idx;

	// index the interior vertices
	idx = 0;
	for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++){
		if(!v->onBoundaryTutte){
			v->indexN = idx;
			indexMapping[v->index] = v->indexN;
			idx++;
		}
	}
	// return indexMapping;
}

void Tutte::initSystem(Mesh& mesh){
	if(H_is_set) {
		throw std::runtime_error("H matrix is already set when trying to initialize it");
	}

	// Construct H \phi = B
	int idx;
	std::vector<Triplet> Mat_coefficients;
	double zero = 0;

	for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++){
		if(!v->onBoundaryTutte) {
			idx = v->indexN;
			Mat_coefficients.push_back(Triplet(idx, idx, zero));
			Bx[idx] = 0.0;
			By[idx] = 0.0;
		}
	}

	HalfEdgeIter he, he2;
	VertexIter v_i, v_j;
	int idx_i, idx_j;
	for (EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); e++){
		he = e->halfEdge();
		he2 = he->flip();
		v_i = he->vertex();
		v_j = he2->vertex();

		if(!v_i->onBoundaryTutte && !v_j->onBoundaryTutte){
			idx_i = v_i->indexN;
			idx_j = v_j->indexN;
			Mat_coefficients.push_back(Triplet(idx_i, idx_j, zero));
			Mat_coefficients.push_back(Triplet(idx_j, idx_i, zero));
		}
	}

	H.setFromTriplets(Mat_coefficients.begin(), Mat_coefficients.end());
	H_is_set = true;
}

/* =============================================================================================
Reset matrix to zero
=============================================================================================== */

void Tutte::resetSystem(Mesh& mesh){
	int idx;
	for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++){
		if(!v->onBoundaryTutte) {
			idx = v->indexN;
			H.coeffRef(idx, idx) = 0;
			Bx[idx] = 0.0;
			By[idx] = 0.0;
		}
	}

	HalfEdgeIter he, he2;
	VertexIter v_i, v_j;
	int idx_i, idx_j;
	for (EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); e++)
	{
		he = e->halfEdge();
		he2 = he->flip();
		v_i = he->vertex();
		v_j = he2->vertex();

		if(!v_i->onBoundaryTutte && !v_j->onBoundaryTutte){
			idx_i = v_i->indexN;
			idx_j = v_j->indexN;
			H.coeffRef(idx_i, idx_j) = 0;
			H.coeffRef(idx_j, idx_i) = 0;
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

void Tutte::uniformWeight(Mesh& mesh){
	std::vector<Triplet> H_coefficients;
	int idx_i, idx_j;
	double weight;
	VertexIter v_j;
	for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++){
		if(!v->onBoundaryTutte){
			idx_i = v->indexN;
			Bx[idx_i] = 0;
			By[idx_i] = 0;
			// deg_i = -v->degree();
			int deg_i = 0;
			
			// loop over all neighbors of v
			weight = 1.0;
			HalfEdgeIter he = v->halfEdge(); 
			do {
				v_j = he->flip()->vertex();
				idx_j = v_j->indexN;

				if(!v_j->onBoundaryTutte) {
					H_coefficients.push_back(Triplet(idx_i, idx_j, weight));
				} else {
					Bx[idx_i] = Bx[idx_i] - weight * Xb[idx_j];
					By[idx_i] = By[idx_i] - weight * Yb[idx_j];
				}
				he= he->flip()->next();
				deg_i++;
			} while (he != v->halfEdge());
			
			H_coefficients.push_back(Triplet(idx_i, idx_i, -deg_i));
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

void Tutte::cotanWeight(Mesh& mesh){
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

	for (EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); e++){
		he = e->halfEdge();
		he2 = he->flip();

		v_i = he->vertex();
		v_j = he2->vertex();

		cot_beta_1 = he->cotan();
		cot_beta_2 = he2->cotan();
		weight = cot_beta_1 + cot_beta_2;
		if(weight < 0) weight = floor;
		
		idx = v_i->indexN;
		jdx = v_j->indexN;

		if( !v_i->onBoundaryTutte && !v_j->onBoundaryTutte) {
			H.coeffRef(idx, jdx) = weight;
			H.coeffRef(jdx, idx) = weight;
			H.coeffRef(idx, idx) -= weight;
			H.coeffRef(jdx, jdx) -= weight;
		} else {
			if(v_i->onBoundaryTutte && !v_j->onBoundaryTutte) {
				H.coeffRef(jdx, jdx) -= weight;
				Bx[jdx] -= weight*Xb[idx];
				By[jdx] -= weight*Yb[idx];
			} else if(!v_i->onBoundaryTutte && v_j->onBoundaryTutte) {
				H.coeffRef(idx, idx) -= weight;
				Bx[idx] -= weight*Xb[jdx];
				By[idx] -= weight*Yb[jdx];
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
	Vector a, b, c, d;

	for (EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); e++)
	{
		he = e->halfEdge();
		he2 = he->flip();

		v_i = he->vertex();
		v_j = he2->vertex();

		a = he->vertex()->position;
		b = he->next()->vertex()->position;
		c = he->prev()->vertex()->position;
		d = he2->prev()->vertex()->position;

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

		if( !v_i->onBoundaryTutte && !v_j->onBoundaryTutte) {
			H.coeffRef(idx, jdx) = weight;
			H.coeffRef(jdx, idx) = weight;
			H.coeffRef(idx, idx) -= weight;
			H.coeffRef(jdx, jdx) -= weight;
		} else {
			if(v_i->onBoundaryTutte && !v_j->onBoundaryTutte) {
				H.coeffRef(jdx, jdx) -= weight;
				Bx[jdx] -= weight*Xb[idx];
				By[jdx] -= weight*Yb[idx];
			} else if(!v_i->onBoundaryTutte && v_j->onBoundaryTutte) {
				H.coeffRef(idx, idx) -= weight;
				Bx[idx] -= weight*Xb[jdx];
				By[idx] -= weight*Yb[jdx];
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
	Vector a, b, c, d;

	for (EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); e++)
	{
		he = e->halfEdge();
		he2 = he->flip();

		v_i = he->vertex();
		v_j = he2->vertex();

		
		a = he->vertex()->position;
		b = he->next()->vertex()->position;
		c = he->prev()->vertex()->position;
		d = he2->prev()->vertex()->position;

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
		if(!v_i->onBoundaryTutte && !v_j->onBoundaryTutte) {
			H.coeffRef(idx, jdx) = weight;
			H.coeffRef(jdx, idx) = weight;
			H.coeffRef(idx, idx) -= weight;
			H.coeffRef(jdx, jdx) -= weight;
		} else {
			if(v_i->onBoundaryTutte && !v_j->onBoundaryTutte) {
				H.coeffRef(jdx, jdx) -= weight;
				Bx[jdx] -= weight*Xb[idx];
				By[jdx] -= weight*Yb[idx];
			} else if(!v_i->onBoundaryTutte && v_j->onBoundaryTutte) {
				H.coeffRef(idx, idx) -= weight;
				Bx[idx] -= weight*Xb[jdx];
				By[idx] -= weight*Yb[jdx];
			}
		}
	}
}

void Tutte::setSquareBoundary(CircularLinkedList& nodesOnGenerators)
{
	/*
	1. iterating over all nodes in nodesOnGenerators,
	   find the root of nodesOnGenerators by finding the vertex with reference index appears 4 times	
	2. start from the root, iterate over all nodes in nodesOnGenerators
	   stop at meet the root at the third time 
	   count the number of nodes traversed between meeting the root at the second time and the third time
	   the two integers n1, n2 are the number of nodes on the first and second generators
	3. use dx = 1/(n1+1) and dy = 1/(1+n2) as the distance between each pair of nodes on the boundary
	4. start from the root, iterate over all nodes in nodesOnGenerators
	   set the coordinates of the nodes on the boundary:
	   1) (0, 0) (1dx, 0) (2dx, 0) ... (1, 0)
	   2) (1, 0) (1, dy) (1, 2dy) ... (1, 1)
	   3) (1, 1) (1-dx, 1) (1-2dx, 1) ... (0, 1)
	   4) (0, 1) (0, 1-dy) (0, 1-2dy) ... (0, 0)
	*/
	 
	
	// find n1, n2 are the number of nodes on the first and second generators except the root
	int n1 = 0, n2 = 0;
	auto nd = nodesOnGenerators.getHead();
	auto v = nd->vertexPtr;
	int rootIndex = (v->referenceIndex == -1) ? v->index : v->referenceIndex;
	// start with the node next to the root
	nd = nd->next;
	bool haveMeetRootOnce = false;
	while (true) {
		v = nd->vertexPtr;
		int index = (v->referenceIndex == -1) ? v->index : v->referenceIndex;
		if (index == rootIndex) {
			// break meet the root at the second time
			if (haveMeetRootOnce) {
				break;
			}
			haveMeetRootOnce = true;
		}
		if (haveMeetRootOnce) {
			n2++;
		} else {
			n1++;
		}
		nd = nd->next;
	}
	n2--; // the root node is counted in n2

	// set the coordinates Xb, Yb of the nodes on the boundary
	Xb.resize(Nb);
	Yb.resize(Nb);
	double dx = 1.0/(n1+1);
	double dy = 1.0/(n2+1);
	nd = nodesOnGenerators.getHead();
	int index = (v->referenceIndex == -1) ? v->index : v->referenceIndex;
	
	// x-axis first, (0, 0) (1dx, 0) (2dx, 0) ... (n1*dx, 0)
	for (int j = 0; j < n1+1; j++) {
		auto v = nd->vertexPtr;
		Xb[v->indexN] = j*dx;
		Yb[v->indexN] = 0;
		nd = nd->next;
	}

	// vertical second, (1, 0) (1, dy) (1, 2dy) ... (1, n2*dy)
	for (int j = 0; j < n2+1; j++) {
		auto v = nd->vertexPtr;
		Xb[v->indexN] = 1;
		Yb[v->indexN] = j*dy;
		nd = nd->next;
	}

	// horizontal third, (1, 1) (1-dx, 1) (1-2dx, 1) ... (0, 1 - n1*dx)
	for (int j = 0; j < n1+1; j++) {
		auto v = nd->vertexPtr;
		Xb[v->indexN] = 1 - j*dx;
		Yb[v->indexN] = 1;
		nd = nd->next;
	}

	// vertical fourth, (0, 1) (0, 1-dy) (0, 1-2dy) ... (0, 1-n2*dy)
	for (int j = 0; j < n2+1; j++) {
		auto v = nd->vertexPtr;
		Xb[v->indexN] = 0;
		Yb[v->indexN] = 1 - j*dy;
		nd = nd->next;
	}

	if (nd != nodesOnGenerators.getHead()) {
		throw std::runtime_error("Error: after assigning coordinates to boundary nodes, the head of nodesOnGenerators is not reached");
	}
}

void printSparseMatrix(const Eigen::SparseMatrix<double>& matrix) {
    std::cout << "Sparse Matrix (" << matrix.rows() << " x " << matrix.cols() << "):" << std::endl;

    for (int k = 0; k < matrix.outerSize(); ++k) { // Iterate over outer dimension
        for (Eigen::SparseMatrix<double>::InnerIterator it(matrix, k); it; ++it) {
            // Print row, column, and value of non-zero entry
            std::cout << "(" << it.row() << ", " << it.col() << ") = " << it.value() << std::endl;
        }
    }
}


/* ===== Tutte Embedding ============================================================================
* Solve for vertex coordinates on the plane using Tutte Embedding
==================================================================================================*/

VertexData<std::pair<double, double>> Tutte::planarTutte(Mesh& mesh, 
CircularLinkedList& nodesOnGenerators, int Tutte_type, bool verbose)
{
	VertexData<std::pair<double, double>> vertexUV_positions(mesh, std::make_pair(0.0, 0.0));

	// assign index to vertices to differentiate between boundary and interior vertices
	indexVertices(mesh, nodesOnGenerators);
	setSquareBoundary(nodesOnGenerators);

	/*
	====================================================================================
	Define Tutte Embedding System
	====================================================================================
	*/
	// Set matrices for solving system of linear equations 
	int n_vertices = mesh.vertices.size();
	// # of interior points
	int n_active = n_vertices - Nb; 
	H.resize(n_active, n_active);
	Bx.resize(n_active);
	By.resize(n_active);
	// Bx.setZero();
	// By.setZero();
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

	// print Bx and By
	if (verbose){
		std::cout << " " << std::endl;
		std::cout << "Tutte embedding: " << std::endl;
		std::cout << "=================" << std::endl;
		std::cout << " " << std::endl;
		std::cout << "Number of boundary points                            : " << Nb << std::endl;
		std::cout << "Edge weight for embedding                            : " << type << std::endl;
		std::cout << " " << std::endl;
		// TODO: Bx is not set correctly in Python but is fine in C++
		IC(Bx);
		IC(By);
		// print the sparse matrix H
		std::cout << "Sparse matrix H: " << std::endl;
		std::cout << "=================" << std::endl;
		printSparseMatrix(H);
		std::cout << " " << std::endl;
	}
	
	

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
		std::cout << "Number of boundary points                            : " << Nb << std::endl;
		std::cout << "Edge weight for embedding                            : " << type << std::endl;
		std::cout << "Error on X coordinates                               :" << err1 << std::endl;
		std::cout << "Error on Y coordinates                               :" << err2 << std::endl;
		std::cout << " " << std::endl;	
	}

	if(err1 > 1.e-5 || err2 > 1.e-5) {
		std::cout << "Warning: Large error when solving Tutte Embedding" << std::endl;
	}

	// store coordinates of all vertices
	for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++){
		if(!v->onBoundaryTutte){
			vertexUV_positions[v] = std::make_pair(X[v->indexN], Y[v->indexN]);
		} else {
			vertexUV_positions[v] = std::make_pair(Xb[v->indexN], Yb[v->indexN]);
		}
	}
	return vertexUV_positions;
}// end of Tutte::planarTutte

}// end of namespace geomp