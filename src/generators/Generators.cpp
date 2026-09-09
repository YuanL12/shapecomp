// #include <queue>
// #include <unordered_set>
// #include <cstdint>
// #include <Eigen/Sparse>
// #include <Eigen/SparseLU>
// #include <Eigen/Dense>
#include "geomp/generators/Generators.h"
#include "geomp/generators/ReebLoop.h"
// #include "geomp/matrix/laplacian.h"

namespace geomp {

CircularLinkedList::CircularLinkedList(const std::vector<EdgeCIter>& edges) {
    for (EdgeCIter e : edges) {
        insert(e);
    }
}

CircularLinkedList::CircularLinkedList(const std::vector<HalfEdgeIter>& halfEdges) {
    if (halfEdges.empty() || halfEdges.size() < 2) {
        throw std::runtime_error("Cannot build CircularLinkedList from an empty sequence or a single half-edge.");
    }

    HalfEdgeIter he = halfEdges[0]; // u->v
    tail = new VertexNode(he->vertex()); // u
    head = new VertexNode(he->flip()->vertex()); // v
    tail->next = head;
    head->prev = tail;
    VertexNode* lastNode = head; 
    for (int i = 1; i < halfEdges.size(); i++) {
        he = halfEdges[i]; // u->v
        VertexIter u = he->vertex();
        VertexIter v = he->flip()->vertex();

        // sanity check the next starts at the current end of the walk.
        if (lastNode->vertexPtr != u) {
            throw std::runtime_error("The two halfedges are not connected");
        }

        if (i == halfEdges.size() - 1) {
            // now for the last halfedge, no need to create a new node, do 
            // sanity check the last v is tail and u is head
            if (v != tail->vertexPtr) {
                throw std::runtime_error("The last halfedge is not connected to the tail");
            }

            // connect the last node to the tail
            lastNode->next = tail;
            tail->prev = lastNode;
        } else {
            // create an new node for v
            VertexNode* newNode = new VertexNode(v);

            // connect the new node to the last node
            newNode->prev = lastNode;
            lastNode->next = newNode;

            // update the last node to the new node
            lastNode = newNode;
        }
    }
}

void CircularLinkedList::print() const {
    if (head == nullptr) {
        std::cout << "Empty CircularLinkedList" << std::endl;
        return;
    }
    int size = 0;
    VertexNode* current = head;
    do {
        int index =
            (current->vertexPtr->referenceIndex == -1) ? current->vertexPtr->index : current->vertexPtr->referenceIndex;
        std::cout << index << "->";
        // int index = current->vertexPtr->index;
        // int refIndex = current->vertexPtr->referenceIndex;
        // std::cout << index << "/" << refIndex << " ";
        current = current->next;
        size++;
    } while (current != head);
    std::cout << ", size = " << size << std::endl;
}

CircularLinkedList::~CircularLinkedList() {
    if (head == nullptr) return;

    // Break the circular reference
    if (tail != nullptr) {
        tail->next = nullptr;
    }

    // Properly delete all nodes
    VertexNode* current = head;
    while (current != nullptr) {
        VertexNode* nextNode = current->next;
        delete current;
        current = nextNode;
    }
}

VertexNode* CircularLinkedList::next(VertexNode* v) { return v->next; }

VertexNode* CircularLinkedList::prev(VertexNode* v) { return v->prev; }

VertexNode* CircularLinkedList::findNode(VertexCIter vertex) {
    VertexNode* current = head;
    while (current != nullptr) {
        if (current->vertexPtr == vertex) {
            return current;
        }
        current = current->next;
    }
    return nullptr;
}

void CircularLinkedList::insert(EdgeCIter e) {
    VertexIter u = e->halfEdge()->vertex();
    VertexIter v = e->halfEdge()->flip()->vertex();

    VertexNode* uNode = findNode(u);
    VertexNode* vNode = findNode(v);

    if (uNode == nullptr && vNode == nullptr) {
        // Neither u nor v appears, create new nodes
        uNode = new VertexNode(u);
        vNode = new VertexNode(v);
        uNode->next = vNode;
        vNode->prev = uNode;
        if (head == nullptr) {
            head = uNode;
            tail = vNode;
        } else {
            tail->next = uNode;
            uNode->prev = tail;
            tail = vNode;
        }
    } else if (uNode != nullptr && vNode == nullptr) {
        // Only u appears
        vNode = new VertexNode(v);
        if (!uNode->hasNext()) {
            uNode->next = vNode;
            vNode->prev = uNode;
            if (uNode == tail) {
                tail = vNode;
            }
        } else if (!uNode->hasPrev()) {
            uNode->prev = vNode;
            vNode->next = uNode;
            if (uNode == head) {
                head = vNode;
            }
        }
    } else if (uNode == nullptr && vNode != nullptr) {
        // Only v appears
        uNode = new VertexNode(u);
        if (!vNode->hasNext()) {
            vNode->next = uNode;
            uNode->prev = vNode;
            if (vNode == tail) {
                tail = uNode;
            }
        } else if (!vNode->hasPrev()) {
            vNode->prev = uNode;
            uNode->next = vNode;
            if (vNode == head) {
                head = uNode;
            }
        }
    } else {
        // Both u and v appear, connect them if possible
        if (!uNode->hasNext() && !vNode->hasPrev()) {
            uNode->next = vNode;
            vNode->prev = uNode;
        } else if (!uNode->hasPrev() && !vNode->hasNext()) {
            uNode->prev = vNode;
            vNode->next = uNode;
        }
    }
}

void CircularLinkedList::orientedInsert(HalfEdgeCIter he) {
    throw std::runtime_error("Not implemented of orientedInsert");
}

// Constructor
Generators::Generators(const Mesh& mesh) : primalParent(mesh), dualParent(mesh), ngon(mesh) {
    // build primal and dual spanning trees
    buildPrimalSpanningTree(mesh);
    buildDualSpanningTree(mesh);
}

void Generators::buildPrimalSpanningTree(const Mesh& mesh) {
    // mark each vertex as its own parent
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        primalParent[v] = v;
    }

    // build primal spanning tree
    VertexCIter root = mesh.vertices.begin();
    std::queue<VertexCIter> q;
    q.push(root);

    while (!q.empty()) {
        VertexCIter u = q.front();
        q.pop();

        HalfEdgeCIter he = u->halfEdge();
        do {
            HalfEdgeCIter flip = he->flip();
            EdgeCIter e = he->edge();

            if (e->isCuttable) {
                VertexCIter v = flip->vertex();

                if (primalParent[v] == v && v != root) {
                    primalParent[v] = u;
                    q.push(v);
                }
            }

            he = flip->next();
        } while (he != u->halfEdge());
    }
}

bool Generators::inPrimalSpanningTree(EdgeCIter e) const {
    HalfEdgeCIter he = e->halfEdge();
    VertexCIter u = he->vertex();
    VertexCIter v = he->flip()->vertex();

    return primalParent[u] == v || primalParent[v] == u;
}

void Generators::buildDualSpanningTree(const Mesh& mesh) {
    int uncuttableEdges = 0;
    FaceData<std::vector<HalfEdgeCIter>> ngonHalfEdges(mesh);

    // mark each ngon face as its own parent and collect halfedges
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (uncuttableEdges > 0) {
            uncuttableEdges--;
            continue;
        }
        // get the first cuttable edge
        HalfEdgeCIter he = f->halfEdge();
        while (!he->edge()->isCuttable) he = he->next();

        // collect halfedges and ngon face
        HalfEdgeCIter fhe = he;
        std::unordered_set<int> seenUncuttableEdges;
        do {
            ngon[he->face()] = f;
            ngonHalfEdges[f].emplace_back(he);

            he = he->next();
            while (!he->edge()->isCuttable) {
                seenUncuttableEdges.emplace(he->edge()->index);
                he = he->flip()->next();
            }

        } while (he != fhe);

        uncuttableEdges = (int)seenUncuttableEdges.size();
        if (f->fillsHole) uncuttableEdges--;
        dualParent[f] = f;
    }

    // build dual spanning tree by BFS
    FaceCIter root = mesh.faces.begin();
    std::queue<FaceCIter> q;
    q.push(root);

    while (!q.empty()) {
        FaceCIter f = q.front();
        q.pop();

        const std::vector<HalfEdgeCIter>& halfEdges = ngonHalfEdges[f];
        for (int i = 0; i < (int)halfEdges.size(); i++) {
            HalfEdgeCIter he = halfEdges[i];
            EdgeCIter e = he->edge();

            if (!inPrimalSpanningTree(e)) {
                FaceCIter g = ngon[he->flip()->face()];

                if (dualParent[g] == g && g != root) {
                    dualParent[g] = f;
                    q.push(g);
                }
            }
        }
    }
}
bool Generators::inDualSpanningTree(EdgeCIter e) const {
    HalfEdgeCIter he = e->halfEdge();
    FaceCIter f = ngon[he->face()];
    FaceCIter g = ngon[he->flip()->face()];

    return dualParent[f] == g || dualParent[g] == f;
}

EdgeIter sharedEdgeByFaces(FaceCIter f, FaceCIter g) {
    HalfEdgeCIter he = f->halfEdge();
    do {
        if (he->flip()->face() == g) {
            return he->edge();
        }
        he = he->next();
    } while (he != f->halfEdge());

    std::cerr << "Trying to find edge with endpoints f, g, but not found, will stop now!" << std::endl;
    std::abort();
    return he->edge();
}

EdgeIter sharedEdge(VertexCIter u, VertexCIter v) {
    HalfEdgeCIter he = u->halfEdge();

    do {
        if (he->flip()->vertex() == v) {
            return he->edge();
        }

        he = he->flip()->next();
    } while (he != u->halfEdge());

    std::cerr << "Trying to find edge with endpoints u, v, but not found, will stop now!" << std::endl;
    std::abort();
    return he->edge();
}

void Generators::labelShortestGeneratorEdges(Mesh& mesh) {
    // find the edges on the generators by Tree-Cotree algorithm and label them onCut on the original mesh
    std::vector<int> representative_edges_indices = labelTreeCotreeGeneratorsEdges(mesh);
    if (representative_edges_indices.size() != 2) {
        throw std::runtime_error("The number of representative edges is not 2");
    }
    
	// copy the mesh to cut it (onCut Edges are correctly copied too)
    Mesh planarMesh(mesh);
    std::vector<std::pair<int, int>> representative_edges_pairs_after_cut;
    MeshSurgery::cutTorus(planarMesh, representative_edges_indices, representative_edges_pairs_after_cut);

    // reset all edge onCut to false, we will cut the input mesh based on new edges 
    for (EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        e->onCut = false;
    }

    // for each representative edge (u^+,v^+) and its copy (u^-,v^-) after cut
    // find the shortest path btw u^+ and u^- on the planar mesh by Dijkstra's algorithm
    // and v^+ and v^- in the planar mesh
    // pick the shorter one and label them on the original mesh
    int generatorIndex = 1;
    for (std::pair<int, int> edges_pair : representative_edges_pairs_after_cut) {
        EdgeIter e1 = planarMesh.edges.begin() + edges_pair.first;
        EdgeIter e2 = planarMesh.edges.begin() + edges_pair.second;
        VertexIter u_plus = e1->halfEdge()->vertex();
        int u_plus_index_before_cut = u_plus->referenceIndex != -1 ? u_plus->referenceIndex : u_plus->index;
        VertexIter v_plus = e1->halfEdge()->flip()->vertex();
        int v_plus_index_before_cut = v_plus->referenceIndex != -1 ? v_plus->referenceIndex : v_plus->index;
        VertexIter u_minus = e2->halfEdge()->vertex();
        int u_minus_index_before_cut = u_minus->referenceIndex != -1 ? u_minus->referenceIndex : u_minus->index;
        VertexIter v_minus = e2->halfEdge()->flip()->vertex();
        int v_minus_index_before_cut = v_minus->referenceIndex != -1 ? v_minus->referenceIndex : v_minus->index;
        // the u and v might be switched, so we swap them if needed
        if (u_plus_index_before_cut != u_minus_index_before_cut) {
            std::swap(u_minus, v_minus);
        }

        std::pair<double, std::vector<int>> shortestPath = Dijkstra::findShortestPath(planarMesh, u_plus, u_minus);
        std::pair<double, std::vector<int>> shortestPath2 = Dijkstra::findShortestPath(planarMesh, v_plus, v_minus);

        // label the edges on the shorter path to true onCut on the original mesh
        if (shortestPath.first > shortestPath2.first) {
            std::swap(shortestPath, shortestPath2);
        }

        // label the edges on the shorter path to true onCut on the original mesh
        int path_size = shortestPath.second.size();
        for (int i = 0; i < path_size - 1; i++) {
            int v_i_index = shortestPath.second[i];
            int v_i_plus_one_index = shortestPath.second[i + 1];
            VertexIter v_i = planarMesh.vertices.begin() + v_i_index;
            VertexIter v_i_plus_one = planarMesh.vertices.begin() + v_i_plus_one_index;
            // find the edge on the cutted planar mesh
            EdgeIter e_on_planar_mesh = sharedEdge(v_i, v_i_plus_one);
            // get the index of the edge on the original uncutted mesh
            int edge_index_before_cut =
                e_on_planar_mesh->referenceIndex != -1 ? e_on_planar_mesh->referenceIndex : e_on_planar_mesh->index;
            // label its onCut
			EdgeIter e_on_original_mesh = mesh.edges.begin() + edge_index_before_cut;
            e_on_original_mesh->onCut = true;
            e_on_original_mesh->generatorIndex = generatorIndex;
        }
		generatorIndex++;
    }
}

void Generators::labelReebGraphGeneratorEdges(Mesh& mesh) {
    // find the edges on the generators by Reeb graph algorithm and label them onCut on the original mesh
    ReebHanTunLoops reebHanTunLoops(mesh);
    reebHanTunLoops.labelReebGraphGeneratorsOnMesh(mesh);
}

void Generators::labelReebGraphGeneratorEdges(Mesh& mesh, const std::vector<double>& distinctDirection) {
    // find the edges on the generators by Reeb graph algorithm and label them onCut on the original mesh
    Vector3 direction(distinctDirection[0], distinctDirection[1], distinctDirection[2]);
    ReebHanTunLoops reebHanTunLoops(mesh, direction);
    reebHanTunLoops.labelReebGraphGeneratorsOnMesh(mesh);
}

void printPath(const std::vector<FaceCIter>& path) {
    for (FaceCIter f : path) {
        std::cout << f->index << " ";
    }
    std::cout << std::endl;
}

// find the shortest path between two nodes in a dual spanning tree, return [f1, ..., f2]
std::vector<FaceCIter> Generators::findShortestPathInDualTree(FaceCIter f1, FaceCIter f2) const {
    // Find paths from f1 and f2 to the root
    std::vector<FaceCIter> path1, path2;
    while (dualParent[f1] != f1) {
        path1.emplace_back(f1);
        f1 = dualParent[f1];
    }
    path1.emplace_back(f1);  // Add the root

    while (dualParent[f2] != f2) {
        path2.emplace_back(f2);
        f2 = dualParent[f2];
    }
    path2.emplace_back(f2);  // Add the root

    // Reverse the paths to go from root to f1 and f2
    std::reverse(path1.begin(), path1.end());
    std::reverse(path2.begin(), path2.end());

    // Find the lowest common ancestor (LCA)
    int lca_index = 0;
    while (lca_index < path1.size() && lca_index < path2.size() && path1[lca_index] == path2[lca_index]) {
        lca_index++;
    }
    lca_index--;  // The LCA is the last common node

    // Construct the shortest path
    std::vector<FaceCIter> shortestPath;

    // Insert path1 from f1 to LCA
    shortestPath.insert(shortestPath.begin(), path1.rbegin(), path1.rend() - lca_index);

    // Insert path2 from LCA (exclusive) to f2
    shortestPath.insert(shortestPath.end(), path2.begin() + lca_index + 1, path2.end());

    return shortestPath;
}

// find the shortest path between two nodes in a primal spanning tree, return edges
std::vector<EdgeIter> Generators::findShortestPathInSpanningTree(VertexCIter v1, VertexCIter v2) const {
    // Find paths from v1 and v2 to the root, expressed as edges
    std::vector<EdgeIter> path1, path2;
    VertexCIter u = v1;
    while (primalParent[u] != u) {
        VertexCIter v = primalParent[u];
        path1.emplace_back(sharedEdge(u, v));
        u = v;
    }

    u = v2;
    while (primalParent[u] != u) {
        VertexCIter v = primalParent[u];
        path2.emplace_back(sharedEdge(u, v));
        u = v;
    }

    // Reverse the paths to go from root to v1 and v2
    std::reverse(path1.begin(), path1.end());
    std::reverse(path2.begin(), path2.end());

    // Find the lowest common ancestor (LCA)
    int lca_index = 0;
    while (lca_index < path1.size() && lca_index < path2.size() && path1[lca_index] == path2[lca_index]) {
        lca_index++;
    }
    lca_index--;  // The LCA is the last common node

    // Construct the shortest path
    std::vector<EdgeIter> shortestPath;
    if (lca_index == -1) {
        // if no common ancestor, concatenate path1(reversed) and path2
        shortestPath.insert(shortestPath.begin(), path1.rbegin(), path1.rend());
        shortestPath.insert(shortestPath.end(), path2.begin(), path2.end());
        return shortestPath;
    }

    // Insert path1 from v1 to LCA(exclusive)
    shortestPath.insert(shortestPath.begin(), path1.rbegin(), path1.rend() - lca_index - 1);

    // Insert path2 from LCA (exclusive) to v2
    shortestPath.insert(shortestPath.end(), path2.begin() + lca_index + 1, path2.end());

    return shortestPath;
}

EdgeData<double> Generators::findClosedOneForm(const Mesh& mesh, EdgeCIter starting_edge) const {
    // store +1/-1 assign to halfedges in a closed 1-form
    HalfEdgeData<double> halfEdgeSign(mesh, 0.0);

    HalfEdgeCIter he = starting_edge->halfEdge();
    // query two nodes(faces) of the dual spanning tree,
    FaceCIter f1 = he->face();
    FaceCIter f2 = he->flip()->face();

    // find the shortest path from face f1 to f2 (faces)
    std::vector<FaceCIter> shortestPath = findShortestPathInDualTree(f1, f2);

    // give +1/-1 to halfedges in the shortest path
    halfEdgeSign[he] = 1.0;
    halfEdgeSign[he->flip()] = -1.0;
    for (int i = 0; i < shortestPath.size() - 1; i++) {
        FaceCIter f = shortestPath[i];
        FaceCIter g = shortestPath[i + 1];
        // find the halfedge that is shared by faces f and g
        HalfEdgeCIter start_he_of_f = f->halfEdge();
        HalfEdgeCIter he_of_f = f->halfEdge();
        do {
            if (he_of_f->flip()->face() == g) {
                break;
            }
            he_of_f = he_of_f->next();
        } while (he_of_f != start_he_of_f);

        // if not found, stop the program and return error
        if (he_of_f->flip()->face() != g) {
            std::cerr << "Trying to find edge shared by faces f and g at index " << i
                      << ", but not found, program will stop now!" << std::endl;
            std::abort();
        }

        // assign -1 to the halfedge in f and +1 to the one in g
        halfEdgeSign[he_of_f] = -1.0;
        halfEdgeSign[he_of_f->flip()] = 1.0;
    }  // done with assigning +1/-1 to halfedges in the shortest path

    // convert the halfedge sign to a closed 1-form on edges
    EdgeData<double> closedOneForm(mesh, 0.0);
    for (EdgeCIter e = mesh.edges.cbegin(); e != mesh.edges.cend(); e++) {
        HalfEdgeCIter he = e->halfEdge();
        if (halfEdgeSign[he] == 0.0) {
            // if halfedge sign is not assigned, continue to the next edge
            continue;
        }
        VertexCIter u = he->vertex();  // base vertex of the half-edge
        VertexCIter v = he->flip()->vertex();
        if (u->index < v->index) {
            // check before assignment
            if (closedOneForm[e] != 0.0 && closedOneForm[e] != halfEdgeSign[he]) {
                std::cerr << "The closed 1-form is not properly assigned, program will stop now!" << std::endl;
                std::abort();
            }
            closedOneForm[e] = halfEdgeSign[he];
        } else {
            if (closedOneForm[e] != 0.0 && closedOneForm[e] != -halfEdgeSign[he]) {
                std::cerr << "The closed 1-form is not properly assigned, program will stop now!" << std::endl;
                std::abort();
            }
            closedOneForm[e] = -halfEdgeSign[he];
        }
    }
    return closedOneForm;
}

// get all closed 1-forms, each closed 1-form corresponds to a homological generator
// from the Tree-Cotree algorithm
std::vector<Eigen::VectorXd> Generators::computeHarmonicOneForms(const Mesh& mesh) const {
    // initialize the meshHarmonic class
    meshHarmonic meshHarmonic(mesh);

    // loop over edges
    std::vector<Eigen::VectorXd> harmonicOneForms;
    harmonicOneForms.reserve(mesh.edges.size());
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        if (!inPrimalSpanningTree(e) && !inDualSpanningTree(e) && e->isCuttable) {
            // find the closed 1-form corresponding to the homological generator with edge e
            EdgeData<double> closedOneForm = findClosedOneForm(mesh, e);

            // convert the closed 1-form to a dense vector
            Eigen::VectorXd closedOneFormVector(mesh.edges.size());
            for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
                closedOneFormVector[e->index] = closedOneForm[e];
            }

            // compute the harmonic 1-form by the meshHarmonic class
            Eigen::VectorXd harmonicOneFormVector = meshHarmonic.computeHarmonicOneForm(closedOneFormVector);
            std::cout << "Done with computing a harmonic 1-form" << std::endl;
            harmonicOneForms.emplace_back(harmonicOneFormVector);
        }
    }
    return harmonicOneForms;
}

void Generators::labelGeneratorsOnEdgesRaw(Mesh& mesh) const {
    // label generators
    int nGenerators = 0;
    for (EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        if (!inPrimalSpanningTree(e) && !inDualSpanningTree(e) && e->isCuttable) {
            nGenerators++;
            HalfEdgeCIter he = e->halfEdge();

            // track vertices back to the root
            std::vector<EdgeIter> temp1;
            VertexCIter u = he->vertex();
            while (primalParent[u] != u) {
                VertexCIter v = primalParent[u];
                temp1.emplace_back(sharedEdge(u, v));
                u = v;
            }

            std::vector<EdgeIter> temp2;
            u = he->flip()->vertex();
            while (primalParent[u] != u) {
                VertexCIter v = primalParent[u];
                temp2.emplace_back(sharedEdge(u, v));
                u = v;
            }

            // temporarily label edges as cuts
            e->onCut = true;
            for (int i = 0; i < temp1.size(); i++) temp1[i]->onCut = true;
            for (int i = 0; i < temp2.size(); i++) temp2[i]->onCut = true;
        }
    }
    std::cout << "nGenerators = " << nGenerators << std::endl;
}

std::vector<FaceCIter> Generators::findOnceIntersectGenerators(Mesh& mesh) {
    std::vector<FaceCIter> dualCycle;
    std::vector<EdgeIter> cycle;
    // label the first generator
    for (EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        if (!inPrimalSpanningTree(e) && !inDualSpanningTree(e) && e->isCuttable) {
            e->onCut = true;
            e->generatorIndex = 1;
            // get two end points of edge e
            HalfEdgeCIter he = e->halfEdge();
            VertexCIter u = he->vertex();
            VertexCIter v = he->flip()->vertex();

            // find the cyle of edges in the primal spanning tree by nearest common ancestor
            cycle = findShortestPathInSpanningTree(u, v);

            for (EdgeIter e : cycle) {
                e->onCut = true;
            }

            // find the dual cycle of faces starting from the two faces of the edge e
            FaceCIter f1 = he->face();
            FaceCIter f2 = he->flip()->face();
            dualCycle = findShortestPathInDualTree(f1, f2);

            // the next generator will be constructed by subdivision, no need to continue
            break;
        }
    }
    return dualCycle;
}

std::vector<int> Generators::labelTreeCotreeGeneratorsEdges(Mesh& mesh) const {
    std::vector<EdgeIter> cycle;
    int nGenerators = 1;
    // representative edges of generators, used to find the shortest generators after cutting
    std::vector<int> representative_edges_indices;

    for (EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        if (!inPrimalSpanningTree(e) && !inDualSpanningTree(e) && e->isCuttable) {
            e->onCut = true;
            e->generatorIndex = nGenerators;
            // get two end points of edge e
            HalfEdgeCIter he = e->halfEdge();
            VertexCIter u = he->vertex();
            VertexCIter v = he->flip()->vertex();

            // find the cycle of edges in the primal spanning tree by nearest common ancestor
            cycle = findShortestPathInSpanningTree(u, v);
            // label edges as onCut
            for (EdgeIter e : cycle) e->onCut = true;
            nGenerators++;
            representative_edges_indices.emplace_back(e->index);
        }
    }
    return representative_edges_indices;
}

void Generators::setNewFaces(FaceIter f, HalfEdgeIter h1, HalfEdgeIter h2, HalfEdgeIter h3) {
    h1->setNext(h2);
    h2->setNext(h3);
    h3->setNext(h1);
    h1->setPrev(h3);
    h2->setPrev(h1);
    h3->setPrev(h2);
    h1->setFace(f);
    h2->setFace(f);
    h3->setFace(f);
    f->setHalfEdge(h1);
}

// reset the head of nodesOnGenerators to the root of 2 generators of a torus
void resetHeadForTorus(CircularLinkedList& nodesOnGenerators) {
    // find the root of nodesOnGenerators
    std::unordered_map<int, int> vertexCount;
    auto nd = nodesOnGenerators.getHead();
    do {
        auto v = nd->vertexPtr;
        int index = (v->referenceIndex == -1) ? v->index : v->referenceIndex;
        vertexCount[index]++;
        nd = nd->next;
    } while (nd != nodesOnGenerators.getHead());
    // only the root appears 4 times
    int rootIndex = -1;
    for (auto& pair : vertexCount) {
        if (pair.second == 4) {
            rootIndex = pair.first;
            // std::cout << "find root with index: " << rootIndex << std::endl;
        }
        if (pair.second > 4) {
            throw std::runtime_error("Vertex v_" + std::to_string(pair.first) + " appears more than " + std::to_string(pair.second) + " times on the boundary, program will stop now!");
        }
    }
    if (rootIndex == -1) {
        throw std::runtime_error("Could not find the root of the torus, program will stop now!");
    }
    // set as the head of nodesOnGenerators
    nd = nodesOnGenerators.getHead();
    while (true) {
        auto v = nd->vertexPtr;
        int index = (v->referenceIndex == -1) ? v->index : v->referenceIndex;
        if (index == rootIndex) {
            nodesOnGenerators.setHead(nd);
            break;
        }
        nd = nd->next;
    }
}

void CircularLinkedList::reverse() {
    if (head == nullptr) return;
    // start from the head
    VertexNode* current = head;
    VertexNode* prev = head->prev;
    do {
        // fetch the next node
        VertexNode* next = current->next;
        current->next = prev;
        current->prev = next;
        prev = current;
        current = next;
    } while (current != head);
}

}  // namespace geomp
