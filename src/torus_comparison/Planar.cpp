#include "geomp/torus_comparison/Planar.h"

namespace geomp {

PlanarLocator::PlanarLocator(Mesh& mesh, std::optional<std::string> Tutte_embedding_type) {
    mesh_ptr = &mesh;

    // build primal and dual spanning trees
    Generators gen(mesh);
    std::vector<FaceCIter> dualCycle = gen.findOnceIntersectGenerators(mesh);

    // barycenter divide the mesh and label the edges
    subdivided_mesh_ptr = std::make_shared<Mesh>();
    std::pair<std::vector<EdgeCIter>, std::vector<EdgeCIter>> edgeGenerators =
        MeshIO::BarycenterDivideLabelMesh(mesh, *subdivided_mesh_ptr, dualCycle);
    // find the edges on the boundary
    cutted_mesh_ptr = std::make_shared<Mesh>(*subdivided_mesh_ptr);
    std::vector<EdgeCIter> generators_edges_on_planar_mesh = MeshSurgery::cutTorus(*cutted_mesh_ptr);
    // convert edges to ordered linked list of vertices on the boundary
    CircularLinkedList linked_vertices_on_bd = CircularLinkedList(generators_edges_on_planar_mesh);

    // find the two generators by iterating over the linked_vertices_bd
    resetHeadForTorus(linked_vertices_on_bd);
    getGeneratorsFromBoundary(linked_vertices_on_bd);

    // compute Tutte embedding
    Tutte tutte;
    bool verbose = false;
    int Tutte_embedding_type_number = 1;
    if (Tutte_embedding_type.has_value()) {
        Tutte_embedding_type_number = get_Tutte_embedding_type(Tutte_embedding_type.value());
    }
    const auto& uv_postiosns =
        tutte.planarTutte(*cutted_mesh_ptr, linked_vertices_on_bd, Tutte_embedding_type_number, verbose);
    uv_positions_ptr = std::make_shared<UV_TYPE>(uv_postiosns);
}

void PlanarLocator::constructPlanarLocatorShortest(Mesh& mesh, std::optional<std::string> Tutte_embedding_type) {
    // build primal and dual spanning trees
    Generators gen(mesh);
    // label edges on generators on input mesh
    gen.labelShortestGeneratorEdges(mesh);

    // copy the input mesh to cut
    this->cutted_mesh_ptr = std::make_shared<Mesh>(mesh);
    std::vector<EdgeCIter> generators_edges_on_planar_mesh = MeshSurgery::cutTorus(*this->cutted_mesh_ptr);

    // convert edges to ordered linked list of vertices on the boundary
    CircularLinkedList linked_vertices_on_bd = CircularLinkedList(generators_edges_on_planar_mesh);
    resetHeadForTorus(linked_vertices_on_bd);
    getGeneratorsFromBoundary(linked_vertices_on_bd);

    // compute Tutte embedding
    Tutte tutte;
    bool verbose = false;
    int Tutte_embedding_type_number = 1;
    if (Tutte_embedding_type.has_value()) {
        Tutte_embedding_type_number = get_Tutte_embedding_type(Tutte_embedding_type.value());
    }
    const auto& uv_postiosns =
        tutte.planarTutte(*this->cutted_mesh_ptr, linked_vertices_on_bd, Tutte_embedding_type_number, verbose);
    this->uv_positions_ptr = std::make_shared<UV_TYPE>(uv_postiosns);
}

void PlanarLocator::constructPlanarLocatorReebGraph(Mesh& mesh,
    std::optional<std::string> Tutte_embedding_type,
    bool swap_generators) {
    // build primal and dual spanning trees
    Generators gen(mesh);
    // label edges on generators on input mesh
    gen.labelReebGraphGeneratorEdges(mesh);

    // copy the input mesh to cut
    this->cutted_mesh_ptr = std::make_shared<Mesh>(mesh);
    std::vector<EdgeCIter> generators_edges_on_planar_mesh = MeshSurgery::cutTorus(*this->cutted_mesh_ptr);

    // convert edges to ordered linked list of vertices on the boundary
    CircularLinkedList linked_vertices_on_bd = CircularLinkedList(generators_edges_on_planar_mesh);
    resetHeadForTorus(linked_vertices_on_bd);
    getGeneratorsFromBoundary(linked_vertices_on_bd, swap_generators);

    // compute Tutte embedding
    Tutte tutte;
    bool verbose = false;
    int Tutte_embedding_type_number = 1;
    if (Tutte_embedding_type.has_value()) {
        Tutte_embedding_type_number = get_Tutte_embedding_type(Tutte_embedding_type.value());
    }
    const auto& uv_postiosns =
        tutte.planarTutte(*this->cutted_mesh_ptr, linked_vertices_on_bd, Tutte_embedding_type_number, verbose);
    this->uv_positions_ptr = std::make_shared<UV_TYPE>(uv_postiosns);
}

void PlanarLocator::constructPlanarLocatorReebGraphOriented(Mesh& mesh,
    std::optional<std::string> Tutte_embedding_type, 
    std::optional<std::vector<double>> distinctDirection,
    bool swap_generators) {
    // build primal and dual spanning trees
    Generators gen(mesh);
    // label edges on generators on input mesh
    if (distinctDirection.has_value()) {
        gen.labelReebGraphGeneratorEdges(mesh, distinctDirection.value());
    } else {
        gen.labelReebGraphGeneratorEdges(mesh);
    }

    // copy the input mesh to cut
    this->cutted_mesh_ptr = std::make_shared<Mesh>(mesh);
    std::vector<HalfEdgeIter> generators_halfEdges_on_planar_mesh = MeshSurgery::cutTorus2(*this->cutted_mesh_ptr);

    // convert edges to ordered linked list of vertices on the boundary
    CircularLinkedList linked_vertices_on_bd = CircularLinkedList(generators_halfEdges_on_planar_mesh);
    resetHeadForTorus(linked_vertices_on_bd);
    getGeneratorsFromBoundary(linked_vertices_on_bd, swap_generators);

    // compute Tutte embedding
    Tutte tutte;
    bool verbose = false;
    int Tutte_embedding_type_number = 1;
    if (Tutte_embedding_type.has_value()) {
        Tutte_embedding_type_number = get_Tutte_embedding_type(Tutte_embedding_type.value());
    }
    const auto& uv_postiosns =
        tutte.planarTutte(*this->cutted_mesh_ptr, linked_vertices_on_bd, Tutte_embedding_type_number, verbose);
    this->uv_positions_ptr = std::make_shared<UV_TYPE>(uv_postiosns);
}

int PlanarLocator::get_Tutte_embedding_type(std::string type) const {
    int type_number = 1;
    if (type == "Uniform") {
        type_number = 1;
    } else if (type == "coTan") {
        type_number = 2;
    } else if (type == "MeanValue") {
        type_number = 3;
    } else if (type == "Authalic") {
        type_number = 4;
    } else {
        throw std::runtime_error("Invalid Tutte embedding type: " + type);
    }
    return type_number;
}

std::pair<int, Eigen::Vector3d> PlanarLocator::findBaryCoords(const geomp::Vector& P_uv) const {
    // pointer safe check
    if (!cutted_mesh_ptr) {
        throw std::runtime_error("Runtime Error: No mesh contained in PlanarLocator!");
    }

    const auto& uv_positions = *uv_positions_ptr;

    // start with a initial face
    FaceCIter f = cutted_mesh_ptr->faces.begin();
    double area_ABC = area(f);
    double lambda_a, lambda_b, lambda_c;

    /*
    compute the signed area of each triangle from ABC, PBC, PCA, and PAB
    the barycentric coordinates of P = lamda_a * A + lamda_b * B + lamda_c * C
    where lamda_a, lamda_b, and lamda_c are the fractions PBC/ABC, PCA/ABC, and PAB/ABC respectively.
    if lamda_a, lamda_b, and lamda_c are all positive, then P is inside the triangle ABC
    if one of the lamda_a, lamda_b, and lamda_c is negative, say lamda_a,
    then move f to the face adjacent to f opposite to A (the one adjacent to edge BC)
    repeat the process until P is inside the triangle, i.e. all lambda_i are non-negative
    if P is outside the mesh, return -1, raise an error
    */
    int count = 0;
    do {
        // get 3 vertices of the face f
        VertexCIter A = f->halfEdge()->vertex();
        geomp::Vector A_uv = uv_positions[A];
        VertexCIter B = f->halfEdge()->next()->vertex();
        geomp::Vector B_uv = uv_positions[B];
        VertexCIter C = f->halfEdge()->prev()->vertex();
        geomp::Vector C_uv = uv_positions[C];
        std::tie(lambda_a, lambda_b, lambda_c) = barycentric_coordinates(A_uv, B_uv, C_uv, P_uv);

        // get the half edge of the face f
        HalfEdgeCIter AB = f->halfEdge();
        HalfEdgeCIter BC = f->halfEdge()->next();
        HalfEdgeCIter CA = f->halfEdge()->prev();

        // check if P is inside the triangle
        if (lambda_a < 0) {
            f = BC->flip()->face();
        } else if (lambda_b < 0) {
            f = CA->flip()->face();
        } else if (lambda_c < 0) {
            f = AB->flip()->face();
        } else {
            return std::make_pair(f->index, Eigen::Vector3d(lambda_a, lambda_b, lambda_c));
        }
        // stop if iteration is too long
        count++;
        if (count > 100000) {
            std::cout << "Have seen over 100000 faces to find the bary center coordinates, should stop" << std::endl;
            break;
        }
    } while (true);

    throw std::runtime_error("At PlanarLocator::findBaryCoords, Cannot locate the face containing the point P");
    return std::make_pair(-1, Eigen::Vector3d(0.0, 0.0, 0.0));
}

std::pair<int, Eigen::Vector3d> PlanarLocator::findBaryCoords(const std::pair<double, double>& P_uv) const {
    return findBaryCoords(geomp::Vector(P_uv.first, P_uv.second));
};

std::pair<int, Eigen::Vector3d> PlanarLocator::findBaryCoords(const Eigen::Vector2d& point) const {
    return findBaryCoords(geomp::Vector(point(0), point(1)));
};

// (Parallel) find the barycentric coordinates of triangles containing given points
std::vector<std::pair<int, Eigen::Vector3d>> PlanarLocator::findBaryCoordsParallel(
    const Eigen::MatrixXd& points) const {
    // pointer safe check
    if (!cutted_mesh_ptr) {
        throw std::runtime_error("Runtime Error: No mesh contained in PlanarLocator!");
    }

    // find the barycentric coordinates of the triangles containing the points
    // Taskflow loop count
    size_t iteration_count = points.rows();

    std::vector<std::pair<int, Eigen::Vector3d>> results(iteration_count);

    tf::Executor executor;
    tf::Taskflow taskflow;

    taskflow.for_each_index((size_t)0, iteration_count, (size_t)1, [&](size_t i) {
        auto point = points.row(i);
        // convert to geomp::Vector
        geomp::Vector point_vec(point(0), point(1));
        results[i] = findBaryCoords(point_vec);
    });

    executor.run(taskflow).get();  // Execute all tasks
    return results;
}

std::vector<std::pair<int, Eigen::Vector3d>> PlanarLocator::findBaryCoordsParallel(
    const std::vector<std::pair<double, double>>& points) const {
    // pointer safe check
    if (!cutted_mesh_ptr) {
        throw std::runtime_error("Runtime Error: No mesh contained in PlanarLocator!");
    }

    // find the barycentric coordinates of the triangles containing the points
    // Taskflow loop count
    size_t iteration_count = points.size();
    std::vector<std::pair<int, Eigen::Vector3d>> results(iteration_count);

    tf::Executor executor;
    tf::Taskflow taskflow;

    taskflow.for_each_index((size_t)0, iteration_count, (size_t)1, [&](size_t i) {
        std::pair<double, double> point = points[i];
        results[i] = findBaryCoords(point);
    });

    executor.run(taskflow).get();  // Execute all tasks
    return results;
}

std::unordered_map<int, std::vector<int>> PlanarLocator::getIdentificationMap() const {
    // pointer safe check
    if (!cutted_mesh_ptr) {
        throw std::runtime_error(
            "Unable to get getIdentificationMap, because the mesh contained in PlanarLocator is null.");
    }

    std::unordered_map<int, std::vector<int>> identification_map;
    for (VertexCIter v = cutted_mesh_ptr->vertices.begin(); v != cutted_mesh_ptr->vertices.end(); v++) {
        if (v->referenceIndex != -1) {  // e.g. (index = 20, referenceIndex = 4)
            identification_map[v->referenceIndex].push_back(v->index);
        }
    }
    return identification_map;
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> PlanarLocator::getGeneratorPaths3D() const {
    // pointer safe check
    if (!cutted_mesh_ptr) {
        throw std::runtime_error("No cutted mesh in PlanarLocator!");
    }

    const Mesh& mesh = *cutted_mesh_ptr;

    // Helper function to convert edge indices to 3D path
    auto verticesToPath = [&mesh](const std::vector<int>& vertexIndices) -> Eigen::MatrixXd {
        // Convert to Eigen matrix
        Eigen::MatrixXd path(vertexIndices.size(), 3);
        for (size_t i = 0; i < vertexIndices.size(); ++i) {
            path(i, 0) = mesh.vertices[vertexIndices[i]].position.x;
            path(i, 1) = mesh.vertices[vertexIndices[i]].position.y;
            path(i, 2) = mesh.vertices[vertexIndices[i]].position.z;
        }

        return path;
    };

    // Convert both generators to paths
    Eigen::MatrixXd path1 = verticesToPath(generator1);
    Eigen::MatrixXd path2 = verticesToPath(generator2);

    return std::make_pair(path1, path2);
}

void PlanarLocator::getGeneratorsFromBoundary(CircularLinkedList& nodesOnGenerators, bool swap_generators) {
    /*
    1. iterating over all nodes in nodesOnGenerators,
       find the root of nodesOnGenerators by finding the vertex with reference index appears 4 times
    2. start from the root, iterate over all nodes in nodesOnGenerators
       stop at meet the root at the third time
       count the number of nodes traversed between meeting the root at the second time and the third time
       the two integers n1, n2 are the number of nodes on the first and second generators
    3. use dx = 1/(n1+1) and dy = 1/(1+n2) as the distance between each pair of nodes on the boundary
    4. start from the root, iterate over all nodes in nodesOnGenerators
       set the coordinates of the nodes on the boundary:
       1) (0, 0) (1dx, 0) (2dx, 0) ... (1, 0)
       2) (1, 0) (1, dy) (1, 2dy) ... (1, 1)
       3) (1, 1) (1-dx, 1) (1-2dx, 1) ... (0, 1)
       4) (0, 1) (0, 1-dy) (0, 1-2dy) ... (0, 0)
    */

    // find n1, n2 are the number of nodes on the first and second generators except the root
    int n1 = 0, n2 = 0;
    auto nd = nodesOnGenerators.getHead();
    auto v = nd->vertexPtr;
    int rootIndex = (v->referenceIndex == -1) ? v->index : v->referenceIndex;
    // start with the node next to the root
    nd = nd->next;
    bool haveMeetRootOnce = false;
    while (true) {
        v = nd->vertexPtr;
        int index = (v->referenceIndex == -1) ? v->index : v->referenceIndex;
        if (index == rootIndex) {
            // break meet the root at the second time
            if (haveMeetRootOnce) {
                break;
            }
            haveMeetRootOnce = true;
        }
        if (haveMeetRootOnce) {
            n2++;
        } else {
            n1++;
        }
        nd = nd->next;
    }
    n2--;  // the root node is counted in n2

    nd = nodesOnGenerators.getHead();
    generator1.clear();
    generator2.clear();
    // x-axis first, (0, 0) (1dx, 0) (2dx, 0) ... (n1*dx, 0)
    for (int j = 0; j < n1 + 1; j++) {
        auto v = nd->vertexPtr;
        int v_index = (v->referenceIndex == -1) ? v->index : v->referenceIndex;
        generator1.push_back(v_index);
        nd = nd->next;
    }

    auto end_node_of_generator1 = nd;

    // vertical second, (1, 0) (1, dy) (1, 2dy) ... (1, n2*dy)
    for (int j = 0; j < n2 + 1; j++) {
        auto v = nd->vertexPtr;
        int v_index = (v->referenceIndex == -1) ? v->index : v->referenceIndex;
        generator2.push_back(v_index);
        nd = nd->next;
    }

    // determine generator1 and generator2 by the label of the starting edge
    // swap them if needed
    auto e = sharedEdge(nodesOnGenerators.getHead()->vertexPtr, nodesOnGenerators.getHead()->next->vertexPtr);
    bool should_swap = (e->generatorIndex == 2) != swap_generators;
    if (should_swap) {
        std::swap(generator1, generator2);
        // reset the head of nodesOnGenerators to the end node of generator1
        nodesOnGenerators.setHead(end_node_of_generator1);
    }
}

}  // namespace geomp
