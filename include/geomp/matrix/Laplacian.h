#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <stack>
#include <queue>
#include "geomp/mesh/MeshData.h"
#include "geomp/mesh/Types.h"
#include "extern/icecream/icecream.hpp"
namespace geomp {
/*
Class to compute harmonic 1-form on a mesh.
The two matrices faceEdgeIncidenceMatrix and vertexEdgeIncidenceMatrix
share the same column label and index system:
each column is a directed edge with direction (u, v) where u < v.
*/
class meshHarmonic {
public:
	// face-edge incidence matrix
	Eigen::SparseMatrix<double> faceEdgeIncidenceMatrix;
	// vertex-edge incidence matrix
	Eigen::SparseMatrix<double> vertexEdgeIncidenceMatrix; 

	// add static function to construct the two incidence matrices
	static Eigen::SparseMatrix<double> constructFaceEdgeIncidenceMatrix(const Mesh& mesh);
	static Eigen::SparseMatrix<double> constructVertexEdgeIncidenceMatrix(const Mesh& mesh);

	// constructor
	meshHarmonic(const Mesh& mesh) : faceEdgeIncidenceMatrix(constructFaceEdgeIncidenceMatrix(mesh)),
								  vertexEdgeIncidenceMatrix(constructVertexEdgeIncidenceMatrix(mesh)) {}

	// compute harmonic 1-form from a closed 1-form
	Eigen::VectorXd computeHarmonicOneForm(const Eigen::VectorXd& closedOneFormVector) const;

	// compute a UV value on each vertex of the cutted mesh from two harmonic 1-forms
	static VertexData<std::pair<double, double>> findUVpositions(const Mesh& originalMesh, const Mesh& cuttedMesh, 
		const std::pair<Eigen::VectorXd, Eigen::VectorXd>& twoHarmonicOneFormVectors);
private:
	// create a functional f on vertices from a harmonic 1-form h on direced edges, 
	// e.g., for a directed edge (u, v), the value is f(v) = f(u) + h(e)
	// or f(u) = f(v) - h(e)
	VertexData<double> computeVertexValue(const Mesh& originalMesh, const Mesh& cuttedMesh, 
		const Eigen::VectorXd& originalHarmonicOneFormVector) const;
};



} // namespace geomp
