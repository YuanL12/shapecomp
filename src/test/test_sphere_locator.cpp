#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "geomp/mesh/MeshIO.h"
#include "geomp/sphere_comparison/SphereLocator.h"

using namespace geomp;

namespace {

constexpr const char* kDefaultMeshPath = "test/input/sphere/icosahedron.obj";

Mesh& loadSingleMesh(const std::string& inputPath, Model& model) {
	std::string error;
	if (!MeshIO::read(inputPath, model, error)) {
		throw std::runtime_error("Unable to load mesh: " + inputPath + ". " + error);
	}
	if (model.size() != 1) {
		throw std::runtime_error("Expected exactly one connected component in " + inputPath);
	}
	return model[0];
}

}  // namespace

int main(int argc, char** argv) {
	try {
		std::string meshPath = kDefaultMeshPath;
		if (argc > 1) {
			meshPath = argv[1];
		}

		std::cout << "[test_sphere_locator] Loading mesh: " << meshPath << std::endl;
		Model model;
		Mesh& mesh = loadSingleMesh(meshPath, model);
		std::cout << "  vertices=" << mesh.vertices.size()
		          << " edges=" << mesh.edges.size()
		          << " faces=" << mesh.faces.size()
		          << " chi=" << mesh.eulerCharacteristic() << std::endl;

		if (mesh.eulerCharacteristic() != 2) {
			throw std::runtime_error("Expected Euler characteristic 2");
		}

		SphereLocator locator(mesh, std::string("coTan"));
		Eigen::MatrixXd sphere = locator.get_sphere_vertices();
		if (sphere.rows() != static_cast<int>(mesh.vertices.size()) || sphere.cols() != 3) {
			throw std::runtime_error("Sphere vertex matrix has the wrong shape");
		}
		for (int i = 0; i < sphere.rows(); ++i) {
			const double n = sphere.row(i).norm();
			if (!std::isfinite(n) || std::abs(n - 1.0) > 1e-6) {
				throw std::runtime_error("Sphere vertex is not on the unit sphere");
			}
		}

		const int pole = locator.get_pole_index();
		if (pole < 0 || pole >= sphere.rows()) {
			throw std::runtime_error("Invalid pole index");
		}

		Eigen::Vector3d query = sphere.row(0).transpose();
		auto [face, bary] = locator.findBaryCoords(query);
		if (face < 0) {
			throw std::runtime_error("Barycentric lookup failed");
		}
		const double bary_sum = bary.sum();
		if (std::abs(bary_sum - 1.0) > 1e-6) {
			throw std::runtime_error("Barycentric coordinates do not sum to 1");
		}

		std::cout << "  pole=" << pole << " sample face=" << face
		          << " bary=[" << bary[0] << ", " << bary[1] << ", " << bary[2] << "]" << std::endl;
		std::cout << "[test_sphere_locator] ok" << std::endl;
		return 0;
	} catch (const std::exception& e) {
		std::cerr << "[test_sphere_locator] " << e.what() << std::endl;
		return 1;
	}
}
