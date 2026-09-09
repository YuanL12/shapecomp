#include "geomp/matrix/Laplacian.h"

namespace geomp {

/*
*/
VertexData<std::pair<double, double>> meshHarmonic::findUVpositions(
	const Mesh& originalMesh, 
	const Mesh& cuttedMesh, 
	const std::pair<Eigen::VectorXd, Eigen::VectorXd>& twoHarmonicOneFormVectors)
{
	// assign each vertex a value given by the harmonic 1-form
	VertexData<std::pair<double, double>> vertexUV_positions(cuttedMesh, std::make_pair(0.0, 0.0));
	// run bfs of the cutted mesh
	VertexData<char> visited(cuttedMesh, 0); // 0 for unvisited, 1 for visited
	std::queue<VertexCIter> que;
	que.push(cuttedMesh.vertices.begin());
	visited[cuttedMesh.vertices.begin()] = 1;
	while (!que.empty()) {
		VertexCIter u = que.front(); // base vertex
		que.pop();
		// push all unvisited neighbors of u to the top of stack
		HalfEdgeCIter he = u->halfEdge();
		do {
			VertexCIter v = he->flip()->vertex();
			if (visited[v] == 0) {
				que.push(v);

				// get the edge between u and v in the uncutted original mesh
				int u_original_index = u->referenceIndex == -1 ? u->index : u->referenceIndex;
				int v_original_index = v->referenceIndex == -1 ? v->index : v->referenceIndex;
				EdgeCIter e = originalMesh.getEdge(u_original_index, v_original_index);

				// get the harmonic 1-form value of the edge (e_0, e_1) where e_0 < e_1
				double edgeValueX = twoHarmonicOneFormVectors.first[e->index];
				double edgeValueY = twoHarmonicOneFormVectors.second[e->index];

				// assign the value to the unvisited vertex v
				// pn is +1 if the edge is from u to v, -1 otherwise
				double pn = u_original_index < v_original_index ? 1.0 : -1.0;
				vertexUV_positions[v].first = vertexUV_positions[u].first + pn * edgeValueX;
				vertexUV_positions[v].second = vertexUV_positions[u].second + pn * edgeValueY;
				
				// set v as visited
				visited[v] = 1;
			}
			he = he->flip()->next();
		} while (he != u->halfEdge());
	}

	return vertexUV_positions;


}


VertexData<double> meshHarmonic::computeVertexValue(
	const Mesh& originalMesh, 
	const Mesh& cuttedMesh, 
	const Eigen::VectorXd& originalHarmonicOneFormVector) const
{
	// assign each vertex a value given by the harmonic 1-form
	VertexData<double> vertexValue(cuttedMesh, 0.0);
	// run dfs of the cutted mesh
	VertexData<char> visited(cuttedMesh, 0); // 0 for unvisited, 1 for visited
	std::stack<VertexCIter> stack;
	stack.push(cuttedMesh.vertices.begin());
	visited[cuttedMesh.vertices.begin()] = 1;
	while (!stack.empty()) {
		VertexCIter u = stack.top(); // base vertex
		stack.pop();
		// push all unvisited neighbors of u to the top of stack
		HalfEdgeCIter he = u->halfEdge();
		do {
			VertexCIter v = he->flip()->vertex();
			if (visited[v] == 0) {
				stack.push(v);

				// get the edge between u and v in the uncutted original mesh
				int u_original_index = u->referenceIndex == -1 ? u->index : u->referenceIndex;
				int v_original_index = v->referenceIndex == -1 ? v->index : v->referenceIndex;
				EdgeCIter e = originalMesh.getEdge(u_original_index, v_original_index);

				// get the harmonic 1-form value of the edge (e_0, e_1) where e_0 < e_1
				double edgeValue = originalHarmonicOneFormVector[e->index];

				// assign the value to the unvisited vertex v
				// pn is +1 if the edge is from u to v, -1 otherwise
				double pn = u->index < v->index ? 1.0 : -1.0;
				vertexValue[v] = vertexValue[u] + pn * edgeValue;

				// set v as visited
				visited[v] = 1;
			}
			he = he->next();
		} while (he != u->halfEdge());
	}

	return vertexValue;
}

Eigen::VectorXd meshHarmonic::computeHarmonicOneForm(const Eigen::VectorXd& closedOneFormVector) const
{
	// solve the equation Ax = Ab, where A is the vertex-edge incidence matrix
	// and b is the closed 1-form vector, such that x is not equal to b
	// harmonic 1-form is then b - x
	Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solver;
	solver.compute(vertexEdgeIncidenceMatrix);
	Eigen::VectorXd solution = solver.solve(vertexEdgeIncidenceMatrix * closedOneFormVector);
	Eigen::VectorXd harmonicOneFormVector = closedOneFormVector - solution;

	// sanity check
	// check if the harmonic 1-form is non-zero
	if (harmonicOneFormVector.norm() < 1e-10) {
		std::cerr << "The harmonic 1-form is zero, program will stop now!" << std::endl;
		std::abort();
	}	
	// check if the result is zero
	Eigen::VectorXd result1 = faceEdgeIncidenceMatrix * harmonicOneFormVector;
	Eigen::VectorXd result2 = vertexEdgeIncidenceMatrix * harmonicOneFormVector;
	if (result1.norm() > 1e-10 || result2.norm() > 1e-10) {
		std::cerr << "The `harmonic` 1-form does not satisfy the condition, program will stop now!" << std::endl;
		std::abort();
	}
	return harmonicOneFormVector;
}

/*
Construct the face-edge incidence matrix F,
where each column is a dircted edge with direction (u, v) where u < v. 
For a 1-form vector v, Fv[i] is the sum of directed edge values counter-clockwise around face i.
*/
Eigen::SparseMatrix<double> meshHarmonic::constructFaceEdgeIncidenceMatrix(const Mesh& mesh)
{
	Eigen::SparseMatrix<double> faceEdgeIncidenceMatrix;
	// Resize the matrix to |F| rows and |E| columns
	faceEdgeIncidenceMatrix.resize(mesh.faces.size(), mesh.edges.size());

	// Reserve space for non-zero elements
    // each row will have 3 non-zero elements
	faceEdgeIncidenceMatrix.reserve(mesh.faces.size() * 3);
	
	// Loop over all faces to fill the matrix
	for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
		HalfEdgeCIter he = f->halfEdge();
		do {
			EdgeCIter e = he->edge();
			// get the two endpoints of the edge
			VertexCIter u = he->vertex(); // base vertex of the half-edge
			VertexCIter v = he->flip()->vertex();
			if (u->index > v->index){
				faceEdgeIncidenceMatrix.insert(f->index, e->index) = -1.0;
			}else{
				faceEdgeIncidenceMatrix.insert(f->index, e->index) = 1.0;
			}
			// move to the next half-edge in the face f
			he = he->next();
		} while (he != f->halfEdge());
	}
	return faceEdgeIncidenceMatrix;
}


/*
Construct the vertex-edge incidence matrix M,
where each column is a dircted edge with direction (u, v) where u < v. 
For a 1-form vector v, Mv[i] is the sum of directed edge values of emanating from vertex i.
*/ 
Eigen::SparseMatrix<double> meshHarmonic::constructVertexEdgeIncidenceMatrix(const Mesh& mesh)
{
	Eigen::SparseMatrix<double> vertexEdgeIncidenceMatrix;
	// Resize the matrix to |V| rows and |E| columns
	vertexEdgeIncidenceMatrix.resize(mesh.vertices.size(), mesh.edges.size());

	// Reserve space for non-zero elements
	vertexEdgeIncidenceMatrix.reserve(mesh.edges.size() * 2);
	
	// Loop over all edges to fill the matrix
	for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
		// get the two endpoints of the edge by calling e->twoEndpoints()
		std::pair<VertexCIter, VertexCIter> endpoints = e->twoEndpoints();
		VertexCIter u = endpoints.first;
		VertexCIter v = endpoints.second;

		// smaller index will be assigned +1 and larger index will be assigned -1
		vertexEdgeIncidenceMatrix.insert(u->index, e->index) = 1.0;
		vertexEdgeIncidenceMatrix.insert(v->index, e->index) = -1.0;
	}

	return vertexEdgeIncidenceMatrix;
}

} // namespace geomp