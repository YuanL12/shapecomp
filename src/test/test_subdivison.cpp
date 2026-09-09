#include <iostream>
#include <vector>
#include "geomp/mesh/MeshIO.h"
#include "geomp/generators/Generators.h" // Include the header for Generators
#include "geomp/matrix/Laplacian.h"
#include "geomp/mesh/Types.h"
#include "extern/icecream/icecream.hpp"
#include <chrono>
using namespace geomp;


void loadModel(const std::string& inputPath, Model& model,
			   std::vector<uint8_t>& isSurfaceClosed)
{
	std::string error;
	if (MeshIO::read(inputPath, model, error)) {
		int nMeshes = model.size();
		isSurfaceClosed.resize(nMeshes, 0);

		for (int i = 0; i < nMeshes; i++) {
			Mesh& mesh = model[i];
			int nBoundaries = (int)mesh.boundaries.size();
			if (nBoundaries >= 1) {
				// mesh has boundaries, not deal with it for now
				std::cerr << "The input " << inputPath << 
                    " has boundaries, not able to deal with it for now." << std::endl;
                exit(EXIT_FAILURE);
			} else if (nBoundaries == 0) {
				if (mesh.eulerCharacteristic() == 2) {
					// mesh is closed (sphere)
					isSurfaceClosed[i] = 1;
				} else {
					// mesh has handles (non-zero genus)
					// Generators::compute(mesh);
					std::cout << "Model["<<i<<"] has non-zero number of genus." << std::endl;
				}
			}
		}
	} else {
		std::cerr << "Unable to load file: " << inputPath << ". " << error << std::endl;
		exit(EXIT_FAILURE);
	}
}

// find 2 period vectors of torus
Eigen::Matrix2d findTwoPeriodVectors(
	const Mesh& originalMesh,
	const Mesh& cuttedMesh, 
	const VertexData<std::pair<double, double>>& vertices_uv_positions)
{
	// loop over all vertices in the cutted mesh
	// find a vertex with reference index appears >= 3 times
	// it will be used at the 'base' point of the fundamental domain
	std::vector<int> vertexCount;
	vertexCount.resize(originalMesh.vertices.size(), 0);
	int basePointIndex = -1;
	for (VertexCIter v = cuttedMesh.vertices.begin(); v != cuttedMesh.vertices.end(); v++) {
		int index = v->referenceIndex == -1 ? v->index : v->referenceIndex;
		vertexCount[index]++;
		if (vertexCount[index] >= 3) {
			// base point found
			basePointIndex = index;
			break;
		}
	}
	if (basePointIndex == -1) {
		std::cerr << "No vertex appears at least 3 times to searve as the base point of fundamental domain, will stop now!" << std::endl;
		throw std::runtime_error("No vertex appears at least 3 times to searve as the base point of fundamental domain, will stop now!");
	}

	// find the 3 vertices with the same reference index
	std::vector<Vector> base_vertex_uv_positions;
	for (VertexCIter v = cuttedMesh.vertices.begin(); v != cuttedMesh.vertices.end(); v++) {
		int index = v->referenceIndex == -1 ? v->index : v->referenceIndex;
		if (index == basePointIndex) {
			base_vertex_uv_positions.emplace_back(vertices_uv_positions[v]);
		}
	}

	if (base_vertex_uv_positions.size() <= 2) {
		// std::cerr << "The base point does not have 3 vertices, will stop now!" << std::endl;
		// std::abort();
		throw std::runtime_error("Base_vertex_uv_positions should have >=3 vertices, will stop now!");
	}

	Eigen::Matrix2d periodMatrix;
	Vector periodVector1 = base_vertex_uv_positions[1] - base_vertex_uv_positions[0];
    Vector periodVector2 = base_vertex_uv_positions[2] - base_vertex_uv_positions[0];

    // Convert custom Vector to Eigen::Vector2d
    Eigen::Vector2d eigenPeriodVector1 = periodVector1.to2DEigenVector();
    Eigen::Vector2d eigenPeriodVector2 = periodVector2.to2DEigenVector();
    periodMatrix.col(0) = eigenPeriodVector1;
    periodMatrix.col(1) = eigenPeriodVector2;
    return periodMatrix;
}

Eigen::Matrix2d computeLinearTransformationMatrix(
	const Eigen::Vector2d& u1, 
	const Eigen::Vector2d& u2, 
	const Eigen::Vector2d& v1, 
	const Eigen::Vector2d& v2) 
{
    // // Define the original vectors u1 and u2
    // Eigen::Vector2d u1(1.0, 2.0);
    // Eigen::Vector2d u2(3.0, 4.0);

    // // Define the transformed vectors v1 and v2
    // Eigen::Vector2d v1(2.0, 3.0);
    // Eigen::Vector2d v2(4.0, 5.0);

    // Create matrices U and V
    Eigen::Matrix2d U;
    U << u1, u2;

    Eigen::Matrix2d V;
    V << v1, v2;

    // Compute the transformation matrix A
    Eigen::Matrix2d A = V * U.inverse();

    // Print the transformation matrix A
    std::cout << "Transformation matrix A:\n" << A << std::endl;
	return A;
}

void write_uv_obj(
	const std::string& inputPath, 
	const VertexData<std::pair<double, double>>& vertices_uv_positions, 
	const Mesh& cuttedMesh)
{
	// Determine the output path by adding _uv before .obj
    std::string outputPath = inputPath;
    size_t pos = outputPath.find(".obj");
    if (pos != std::string::npos) {
        outputPath.insert(pos, "_uv");
    } else {
        outputPath += "_uv.obj"; // Handle case where inputPath does not end with .obj
    }
	std::cout << "Output Path: "<< outputPath << std::endl;
	std::string error;
	if (MeshIO::writeMeshOBJ_uv(vertices_uv_positions, outputPath, cuttedMesh)) {
		std::cout << "UV mesh is written to " << outputPath << std::endl;
	} else {
		std::cerr << "Unable to write file: " << outputPath << ". " << error << std::endl;
		exit(EXIT_FAILURE);
	}
}

int main(int argc, char* argv[]) {
    std::cout << "Number of arguments: " << argc << std::endl;
    for (int i = 0; i < argc; ++i) {
        std::cout << "Argument " << i << ": " << argv[i] << std::endl;
    }
	std::string inputPath = "../simple_torus.obj"; // path to the input model
	// Override default input path if provided as a command-line argument
    if (argc > 1) {
        inputPath = argv[1];
    }
	std::cout << "Input path: " << inputPath << std::endl;

    // load model
	geomp::Model model;
	std::vector<uint8_t> isSurfaceClosed;
	loadModel(inputPath, model, isSurfaceClosed);
	std::cout << "success load model." << std::endl;
	
	Mesh subdivided_mesh;
	MeshIO::BarycenterDivideMesh(model[0], subdivided_mesh);
    std::cout << "success barycenter divide mesh." << std::endl;

	//Check the sizes of the subdivided mesh
	std::cout << "Original mesh: " << std::endl;
	int nV = model[0].vertices.size();
	int nE = model[0].edges.size();
	int nF = model[0].faces.size();
	std::cout << "Vertices: " << nV << std::endl;
	std::cout << "Edges: " << nE << std::endl;
	std::cout << "Faces: " << nF << std::endl;

	std::cout << "Subdivided mesh: " << std::endl;
	std::cout << "Vertices: " << subdivided_mesh.vertices.size() << std::endl;
	int expected_nV = nV + nE + nF;
	std::cout << "Expected vertices: " << expected_nV;
	std::cout << ", T/F: " << (expected_nV == subdivided_mesh.vertices.size()) << std::endl;
	std::cout << "Edges: " << subdivided_mesh.edges.size() << std::endl;
	int expected_nE = 2*nE + 6*nF;
	std::cout << "Expected edges: " << expected_nE;
	std::cout << ", T/F: " << (expected_nE == subdivided_mesh.edges.size()) << std::endl;
	std::cout << "Faces: " << subdivided_mesh.faces.size() << std::endl;
	int expected_nF = 6*nF;
	std::cout << "Expected faces: " << expected_nF;
	std::cout << ", T/F: " << (expected_nF == subdivided_mesh.faces.size()) << std::endl;

    // write the mesh out
    // std::string outputPath = "../simple_torus_subdivided.obj";
	// MeshIO::writeMeshOBJ(outputPath, subdivided_mesh);
    
    return 0;
}