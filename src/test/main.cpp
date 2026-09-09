#include <iostream>
#include <vector>
#include "geomp/mesh/MeshIO.h"
#include "geomp/generators/Generators.h"
#include "geomp/matrix/Laplacian.h"
#include "geomp/mesh/Types.h"
#include "geomp/torus_comparison/Planar.h"
#include "extern/icecream/icecream.hpp"
#include "geomp/embedding/torus/Tutte.h"
#include <chrono>
#include "geomp/generators/MeshSurgery.h"
using namespace geomp;

// print mesh size info
void printMeshSizeInfo(const Mesh& mesh) {
	std::cout << "Vertices: " << mesh.vertices.size() << std::endl;
	std::cout << "Edges: " << mesh.edges.size() << std::endl;
	std::cout << "Faces: " << mesh.faces.size() << std::endl;
}


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
		std::abort();
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
	const std::string& outputPath, 
	const VertexData<std::pair<double, double>>& vertices_uv_positions, 
	const Mesh& cuttedMesh)
{
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
	// std::string inputPath = "test/input/torus/torus.obj";
	// std::string inputPath = "../test/input/torus/torus.obj"; // path to the input model
	std::string inputPath = "../test/input/torus/torus_subdivision_igl.obj"; // path to the input model
	// Override default input path if provided as a command-line argument
    if (argc > 1) {
        inputPath = argv[1];
    }
	std::cout << "Input path: " << inputPath << std::endl;

    {	
		// load model
		geomp::Model model;
		std::vector<uint8_t> isSurfaceClosed;
		loadModel(inputPath, model, isSurfaceClosed);
		std::cout << "success load model." << std::endl;
		
		printMeshSizeInfo(model[0]);

		// build primal and dual spanning trees
		Generators gen(model[0]);
		std::vector<FaceCIter> dualCycle = gen.findOnceIntersectGenerators(model[0]);

		// barycenter divide the mesh and label the edges
		Mesh subdivided_mesh; // no copy here because BarycenterDivideLabelMesh will make a copy first
		std::pair<std::vector<EdgeCIter>, std::vector<EdgeCIter>> edgeGenerators = MeshIO::BarycenterDivideLabelMesh(model[0], subdivided_mesh, dualCycle);
		
		std::cout << "success barycenter divide mesh. Mesh size info:" << std::endl;
		printMeshSizeInfo(subdivided_mesh);

		// find the vertices on boundary 
		Mesh cutted_mesh(subdivided_mesh); // copy constructor
		CircularLinkedList linked_vertices_bd = MeshSurgery::cutTorus(cutted_mesh);
		std::cout << "success cut torus. Mesh size info:" << std::endl;
		printMeshSizeInfo(cutted_mesh);
		
		// std::cout << "Linked vertices on boundary" << std::endl;
		// linked_vertices_bd.print();

		// Tutte Embedding
		Tutte tutte;
		VertexData<std::pair<double, double>> uv_postiosns = tutte.planarTutte(cutted_mesh, linked_vertices_bd);

		// Planar Locator
		PlanarLocator planarLocator(cutted_mesh, uv_postiosns);

		// PlanarLocator planarLocator = buildPlanarLocator(model[0]);
		std::pair<int, Eigen::Vector3d> res = planarLocator.findBaryCoords(Vector(0.5, 0.5));
		IC(res.first);
		IC(res.second);
	}

	{
		std::cout << "test build planar locator directly from mesh" << std::endl;
		// load model
		geomp::Model model;
		std::vector<uint8_t> isSurfaceClosed;
		loadModel(inputPath, model, isSurfaceClosed);
		std::cout << "success load model." << std::endl;
		
		printMeshSizeInfo(model[0]);

		// build planar locator directly
		PlanarLocator planarLocator(model[0]);
		std::cout << "success build planar locator." << std::endl;

		// get the cutted mesh
		auto cutted_mesh_ptr = planarLocator.get_cutted_mesh_ptr();
		printMeshSizeInfo(*cutted_mesh_ptr);

		// get the subdivided mesh
		auto subdivided_mesh_ptr = planarLocator.get_subdivided_mesh_ptr();
		printMeshSizeInfo(*subdivided_mesh_ptr);
	}
	
	// // write the uv to an obj file
	// std::string outputPath = inputPath;
    // size_t pos = outputPath.find(".obj");
    // if (pos != std::string::npos) {
    //     outputPath.insert(pos, "_square_uv");
    // } else {
    //     outputPath += "_square_uv.obj"; // Handle case where inputPath does not end with .obj
    // }
	// write_uv_obj(outputPath, uv_postiosns, subdivided_mesh);

    return 0;
}