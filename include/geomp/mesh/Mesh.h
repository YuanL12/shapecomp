#pragma once

#include "geomp/mesh/GeometryUtils.h"
#include "geomp/mesh/CutIterator.h"

namespace geomp {

class Mesh {
public:
	// constructor
	Mesh();

	// copy constructor
	Mesh(const Mesh& mesh);

	// move constructor
	Mesh(Mesh&& mesh) noexcept;

	// returns euler characteristic
	int eulerCharacteristic() const;

	// returns mesh diameter
	double diameter() const;

	// projects UVs to PCA axis
	void projectUvsToPcaAxis();

	// orients UVs to minimize bounding box size
	void orientUvsToMinimizeBoundingBox(int nRotations);

	// computes the ratio of the surface areas of the mesh and flattened mesh
	double areaRatio() const;

	// range based for loop over the cut boundary
	CutPtrSet cutBoundary() const; // valid only for 1 boundary loop, by default, this is the first loop

	// returns reference to wedges (a.k.a. corners)
	std::vector<Wedge>& wedges();
	const std::vector<Wedge>& wedges() const;

	// error (and warning) codes
	enum class ErrorCode {
		ok,
		nonManifoldEdges,
		nonManifoldVertices,
		isolatedVertices
	};

	// get the edge between vertices indices u_id and v_id
	EdgeCIter getEdge(int u_id, int v_id) const;
	EdgeIter getEdge(int u_id, int v_id); 

	// get the vertices, return Eigen::MatrixXd for python
	Eigen::MatrixXd getVertices() const;

	// get the faces, return Eigen::MatrixXi for python
	Eigen::MatrixXi getFaces() const;

	// get the edges, return Eigen::MatrixXi for python
	Eigen::MatrixXi getEdges() const;

	// member variables
	std::vector<Vertex> vertices;
	std::vector<Edge> edges;
	std::vector<Face> faces;
	std::vector<Corner> corners;
	std::vector<HalfEdge> halfEdges;
	std::vector<Face> boundaries;
	double radius;
	Vector cm;
	ErrorCode status;
};

class Model {
public:
	// returns the number of meshes in the model
	int size() const { return (int)meshes.size(); }

	// returns the number of vertices in the model; NOTE: duplicated dofs not counted
	int nVertices() const { return (int)modelToMeshMap.size(); }

	// mesh access
	Mesh& operator[](int index) { return meshes[index]; }
	const Mesh& operator[](int index) const { return meshes[index]; }

	// vertex access
	std::pair<int, int> localVertexIndex(int index) const { return modelToMeshMap[index]; }
	int globalVertexIndex(int mesh, int index) const { return meshToModelMap[mesh][index]; }

	// clears members
	void clear() {
		meshes.clear();
		modelToMeshMap.clear();
		meshToModelMap.clear();
	}

	// members
	std::vector<Mesh> meshes;
	std::vector<std::pair<int, int>> modelToMeshMap;
	std::vector<std::vector<int>> meshToModelMap;
};

} // namespace geomp
