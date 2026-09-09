#pragma once

#include <fstream>
#include <sstream>
#include "geomp/mesh/PolygonSoup.h"
#include "geomp/mesh/MeshData.h"
#include "geomp/mesh/Types.h"
#include <unordered_set>
#include <queue>
#include <optional>

namespace geomp {

class MeshIO {
public:
	// reads polygon soup from file and builds model
	static bool read(
		const std::string& fileName, 
		Model& model, 
		std::string& error, 
		std::optional<bool> normalize = std::nullopt);

	// reads polygon soup from obj file
	static bool readOBJ(const std::string& fileName, PolygonSoup& soup,
						std::vector<std::pair<int, int>>& uncuttableEdges,
						std::string& error);
#ifdef USE_USD
	// reads polygon soup from USD file
	static bool readUSD(const std::string& fileName, PolygonSoup& soup,
						std::vector<std::pair<int, int>>& uncuttableEdges,
						std::string& error);
#endif
	static bool writeMeshOBJ(const std::string& fileName, const Mesh& mesh);
	static bool writeMeshOBJ_uv(const VertexData<std::pair<double, double>>& uv_positions, 
								const std::string& fileName, const Mesh& mesh);

	static void BarycenterDivideMesh(const Mesh& mesh, Mesh& newMesh);
	static Mesh BarycenterDivideMesh(const Mesh& mesh){
		Mesh newMesh;
		BarycenterDivideMesh(mesh, newMesh);
		return newMesh;
	}
	static std::pair<std::vector<EdgeCIter>, std::vector<EdgeCIter>> 
				BarycenterDivideLabelMesh(
				const Mesh& mesh, Mesh& newMesh,
				const std::vector<FaceCIter>& dualCycle);
					
	// builds model
	static bool buildModel(const std::vector<std::pair<int, int>>& uncuttableEdges,
						   PolygonSoup& soup, Model& model, std::string& error, 
						   std::optional<bool> normalize = std::nullopt);

	// writes model positions and UVs to obj file
	static bool writeOBJ(const std::string& fileName, bool writeOnlyUvs,
						 const std::vector<Vector>& positions,
						 const std::vector<Vector>& uvs,
						 const std::vector<int>& vIndices,
						 const std::vector<int>& uvIndices,
						 const std::vector<int>& indicesOffset);
						 
	// Read mesh from Eigen matrices
	static bool readEigenMatrices(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces,
									PolygonSoup& soup, std::vector<std::pair<int, int>>& uncuttableEdges,
									std::string& error);

private:
	static void setNewFace( FaceIter f, 
							HalfEdgeIter h1, 
							HalfEdgeIter h2,
							HalfEdgeIter h3);

	// separates model into components
	static void separateComponents(const PolygonSoup& soup, int nComponents,
								   const std::vector<int>& faceComponent,
								   const std::vector<uint8_t>& isCuttableModelEdge,
								   std::vector<PolygonSoup>& soups,
								   std::vector<std::vector<uint8_t>>& isCuttableSoupEdge,
								   std::vector<std::pair<int, int>>& modelToMeshMap,
								   std::vector<std::vector<int>>& meshToModelMap);

	// preallocates mesh elements
	static void preallocateElements(const PolygonSoup& soup, Mesh& mesh);

	// checks if mesh has isolated vertices
	static bool hasIsolatedVertices(const Mesh& mesh);

	// checks if mesh has non-manifold vertices
	static bool hasNonManifoldVertices(const Mesh& mesh);

	// builds a halfedge mesh
	static bool buildMesh(const PolygonSoup& soup,
						  const std::vector<uint8_t>& isCuttableEdge,
						  Mesh& mesh, std::string& error);

	// centers model around origin through center of mass and rescales to unit radius
	static void normalize_radius(Model& model, bool rescale_to_unit_radius = false);

	// centers model around origin through center of mass and rescales to unit area
	static void normalize_area(Model& model);
};

} // namespace geomp
