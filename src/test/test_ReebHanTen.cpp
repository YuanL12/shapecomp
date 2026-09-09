#include <iostream>
#include <stdexcept>
#include <string>

#include "geomp/generators/ReebLoop.h"
#include "geomp/generators/Generators.h"
#include "geomp/mesh/MeshIO.h"

using namespace geomp;

namespace {

constexpr const char* kDefaultMeshPath = "test/input/torus/torus.obj";
// TODO: fix the problem when using the simple torus (it will not work even if I perturb the vertices or subdivide the mesh)

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
    std::string meshPath = kDefaultMeshPath;
    if (argc > 1) {
        meshPath = argv[1];
    }

    std::cout << "[test_ReebHanTen] Loading mesh: " << meshPath << std::endl;
    Model model;
    Mesh& mesh = loadSingleMesh(meshPath, model);
    
    ReebHanTunLoops reebHanTunLoops(mesh);
    auto [handleVertices, tunnelVertices] = reebHanTunLoops.findReebGraphGenerators();
    std::cout << "Handle vertices: (size = " << handleVertices.size() << ")" << std::endl;
    for (int vertex : handleVertices) {
        std::cout << vertex << ", ";
    }
    std::cout << std::endl;
    std::cout << "Tunnel vertices: (size = " << tunnelVertices.size() << ")" << std::endl;
    for (int vertex : tunnelVertices) {
        std::cout << vertex << ", ";
    }
    std::cout << std::endl;
    return 0;
}
