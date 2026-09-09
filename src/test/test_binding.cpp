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

int main(int argc, char* argv[]) {
    std::cout << "Number of arguments: " << argc << std::endl;
    for (int i = 0; i < argc; ++i) {
        std::cout << "Argument " << i << ": " << argv[i] << std::endl;
    }
	std::string inputPath = "../test/input/torus/simple_torus.obj"; // path to the input model
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

	PlanarLocator planarLocator(model[0]);
	std::pair<int, Eigen::Vector3d> res = planarLocator.findBaryCoords(Vector(0.5, 0.5));
	IC(res.first);
	IC(res.second);
	
	auto UV_pos = planarLocator.get_uv_positions();
	IC(UV_pos.size());
    return 0;
}
