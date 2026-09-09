#pragma once
#include <iostream>
#include <vector>
#include "geomp/mesh/MeshIO.h"
#include "geomp/mesh/Types.h"
#include <limits>
#include <Eigen/Dense>

using namespace geomp;

void loadModel(
	const std::string& inputPath, 
	Model& model,
	std::vector<uint8_t>& isSurfaceClosed,
	std::optional<bool> normalize = std::nullopt
)
{
	std::string error;
	if (MeshIO::read(inputPath, model, error, normalize)) {
		int nMeshes = model.size();
		isSurfaceClosed.resize(nMeshes, 0);

		for (int i = 0; i < nMeshes; i++) {
			Mesh& mesh = model[i];
			int nBoundaries = (int)mesh.boundaries.size();
			if (nBoundaries >= 1) {
				// mesh has boundaries, not deal with it for now
				throw std::runtime_error("The input " + inputPath + " has boundaries, not able to deal with it for now.");
			} else if (nBoundaries == 0) {
				if (mesh.eulerCharacteristic() == 2) {
					// mesh is closed (sphere)
					isSurfaceClosed[i] = 1;
				} 
				// mesh has handles (non-zero genus)
				// Generators::compute(mesh);
			}
		}
	} else {
		throw std::runtime_error("Unable to load file: " + inputPath + ". " + error);
	}
}

void loadModel(
	const Eigen::MatrixXd& vertices, 
	const Eigen::MatrixXi& faces, 
	Model& model,
	std::optional<bool> normalize = std::nullopt)
{
	// read polygon soup from Eigen::MatrixXd and Eigen::MatrixXi
	PolygonSoup soup;
	std::vector<std::pair<int, int>> uncuttableEdges;
	std::string error;
	MeshIO::readEigenMatrices(vertices, faces, soup, uncuttableEdges, error);
	if (error.size() > 0) {
		std::cerr << "Error reading polygon soup from Eigen::MatrixXd and Eigen::MatrixXi: " << error << std::endl;
		throw std::runtime_error(error);
	}

	// build model
		if (!MeshIO::buildModel(uncuttableEdges, soup, model, error, normalize)) {
		std::cerr << "Error building model from polygon soup: " << error << std::endl;
		throw std::runtime_error(error);
	}

	// check if there are more than one connected components
	if (model.size() != 1) {
        throw std::runtime_error("There are " + std::to_string(model.size()) + " connected components in the input path.");
    }
}