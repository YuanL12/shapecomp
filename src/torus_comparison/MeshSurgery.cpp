#include "geomp/generators/MeshSurgery.h"
namespace geomp {

void MeshSurgery::cutTorus(Mesh& mesh,
    std::vector<int> representative_edges_indices,
    std::vector<std::pair<int, int>>& representative_edges_pairs_after_cut) {
    int nVertices = (int)mesh.vertices.size();
    int nEdges = (int)mesh.edges.size();
    int nHalfEdges1 = (int)mesh.halfEdges.size();
    int nV = 0;
    int nE = 0;
    int nHe = 0;

    // insert vertices, edges and halfedges
    VertexData<uint8_t> seenVertex(mesh, 0);
    EdgeData<uint8_t> seenEdge(mesh, 0);
    const CutPtrSet& cutBoundary = mesh.cutBoundary();
    for (const WedgeIter& w : cutBoundary) {
        VertexIter v = w->vertex();
        HalfEdgeIter he = w->halfEdge()->next();
        EdgeIter e = he->edge();
        int generatorIndex = e->generatorIndex;
        int e_index = e->index;
        // insert vertex
        if (!seenVertex[v]) {
            seenVertex[v] = 1;
        } else {
            int index = v->index;
            const Vector& position = v->position;
            v = mesh.vertices.emplace(mesh.vertices.end(), Vertex(&mesh));
            v->position = position;
            v->index = nVertices + nV++;
            v->referenceIndex = index;
        }

        //  insert edge
        if (!seenEdge[e]) {
            seenEdge[e] = 1;
        } else {
            // find the copy of the edge in the representativeEdges
            bool found = false;
            int representative_index = -1;
            for (int e2_index : representative_edges_indices) {
                if (e->index == e2_index) {
                    found = true;
                    representative_index = e2_index;
                    break;
                }
            }
            // insert the edge to the mesh
            e = mesh.edges.emplace(mesh.edges.end(), Edge(&mesh));
            e->index = nEdges + nE++;
            e->referenceIndex = e_index;
            e->generatorIndex = generatorIndex;
            if (found) {
                representative_edges_pairs_after_cut.emplace_back(representative_index, e->index);
            }
        }

        // insert halfedge
        HalfEdgeIter newHe = mesh.halfEdges.emplace(mesh.halfEdges.end(), HalfEdge(&mesh));
        newHe->index = nHalfEdges1 + nHe++;
        newHe->setVertex(v);
        newHe->setEdge(e);
        newHe->setFlip(he);
    }

    // insert boundary face
    FaceIter newF = mesh.boundaries.emplace(mesh.boundaries.end(), Face(&mesh));
    newF->setHalfEdge(mesh.halfEdges.begin() + nHalfEdges1);
    newF->index = 0;

    // update connectivity
    // set the prev and next of halfedges on the boundary
    int nHalfEdges2 = (int)mesh.halfEdges.size();
    for (int j = nHalfEdges1; j < nHalfEdges2; j++) {
        int i = j == nHalfEdges1 ? nHalfEdges2 - 1 : j - 1;
        int k = j == nHalfEdges2 - 1 ? nHalfEdges1 : j + 1;

        HalfEdgeIter newPrev = mesh.halfEdges.begin() + i;
        HalfEdgeIter newHe = mesh.halfEdges.begin() + j;
        HalfEdgeIter newNext = mesh.halfEdges.begin() + k;
        HalfEdgeIter flip = newHe->flip();

        // update vertex connectivity
        VertexIter v = newHe->vertex();
        v->setHalfEdge(newHe);
        HalfEdgeIter h = flip->next();
        do {
            h->setVertex(v);

            if (h->edge()->onCut) break;
            h = h->flip()->next();
        } while (true);

        // update edge connectivity
        EdgeIter e = newHe->edge();
        e->setHalfEdge(flip);
        e->onGenerator = true;

        // update face connectivity
        newHe->setFace(newF);

        // update halfedge connectivity
        flip->setFlip(newHe);
        newHe->setPrev(newPrev);
        newHe->setNext(newNext);
        newHe->onBoundary = true;
    }

    // assign halfedge edges and remove cut label from all edges
    for (EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        e->halfEdge()->setEdge(e);
        e->onCut = false;
    }
}

std::vector<EdgeCIter> MeshSurgery::cutTorus(Mesh& mesh) {
    int nVertices = (int)mesh.vertices.size();
    int nEdges = (int)mesh.edges.size();
    int nHalfEdges1 = (int)mesh.halfEdges.size();
    int nV = 0;
    int nE = 0;
    int nHe = 0;

    // insert vertices, edges and halfedges
    VertexData<uint8_t> seenVertex(mesh, 0);
    EdgeData<uint8_t> seenEdge(mesh, 0);
    const CutPtrSet& cutBoundary = mesh.cutBoundary();
    for (const WedgeIter& w : cutBoundary) {
        VertexIter v = w->vertex();
        HalfEdgeIter he = w->halfEdge()->next();
        EdgeIter e = he->edge();
        int generatorIndex = e->generatorIndex;

        // insert vertex
        if (!seenVertex[v]) {
            seenVertex[v] = 1;
        } else {
            int index = v->index;
            const Vector& position = v->position;
            v = mesh.vertices.emplace(mesh.vertices.end(), Vertex(&mesh));
            v->position = position;
            v->index = nVertices + nV++;
            v->referenceIndex = index;
        }

        //  insert edge
        if (!seenEdge[e]) {
            seenEdge[e] = 1;
        } else {
            e = mesh.edges.emplace(mesh.edges.end(), Edge(&mesh));
            e->index = nEdges + nE++;
            e->generatorIndex = generatorIndex;
        }

        // insert halfedge
        HalfEdgeIter newHe = mesh.halfEdges.emplace(mesh.halfEdges.end(), HalfEdge(&mesh));
        newHe->index = nHalfEdges1 + nHe++;
        newHe->setVertex(v);
        newHe->setEdge(e);
        newHe->setFlip(he);
    }

    // insert boundary face
    FaceIter newF = mesh.boundaries.emplace(mesh.boundaries.end(), Face(&mesh));
    newF->setHalfEdge(mesh.halfEdges.begin() + nHalfEdges1);
    newF->index = 0;

    std::vector<EdgeCIter> edgeOnGenerators;
    // update connectivity
    // set the prev and next of halfedges on the boundary
    int nHalfEdges2 = (int)mesh.halfEdges.size();
    for (int j = nHalfEdges1; j < nHalfEdges2; j++) {
        int i = j == nHalfEdges1 ? nHalfEdges2 - 1 : j - 1;
        int k = j == nHalfEdges2 - 1 ? nHalfEdges1 : j + 1;

        HalfEdgeIter newPrev = mesh.halfEdges.begin() + i;
        HalfEdgeIter newHe = mesh.halfEdges.begin() + j;
        HalfEdgeIter newNext = mesh.halfEdges.begin() + k;
        HalfEdgeIter flip = newHe->flip();

        // update vertex connectivity
        VertexIter v = newHe->vertex();
        v->setHalfEdge(newHe);
        HalfEdgeIter h = flip->next();
        do {
            h->setVertex(v);

            if (h->edge()->onCut) break;
            h = h->flip()->next();
        } while (true);

        // update edge connectivity
        EdgeIter e = newHe->edge();
        e->setHalfEdge(flip);
        e->onGenerator = true;
        edgeOnGenerators.emplace_back(e);

        // update face connectivity
        newHe->setFace(newF);

        // update halfedge connectivity
        flip->setFlip(newHe);
        newHe->setPrev(newPrev);
        newHe->setNext(newNext);
        newHe->onBoundary = true;
    }

    // assign halfedge edges and remove cut label from all edges
    for (EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        e->halfEdge()->setEdge(e);
        e->onCut = false;
    }

    // return the edges on the two generators/boundaries after cutting
    return edgeOnGenerators;
}


std::vector<HalfEdgeIter> MeshSurgery::cutTorus2(Mesh& mesh) {
    int nVertices = (int)mesh.vertices.size();
    int nEdges = (int)mesh.edges.size();
    int nHalfEdges1 = (int)mesh.halfEdges.size();
    int nV = 0;
    int nE = 0;
    int nHe = 0;

    // insert vertices, edges and halfedges
    VertexData<uint8_t> seenVertex(mesh, 0);
    EdgeData<uint8_t> seenEdge(mesh, 0);
    const CutPtrSet& cutBoundary = mesh.cutBoundary();
    for (const WedgeIter& w : cutBoundary) {
        VertexIter v = w->vertex();
        HalfEdgeIter he = w->halfEdge()->next();
        EdgeIter e = he->edge();
        int generatorIndex = e->generatorIndex;

        // insert vertex
        if (!seenVertex[v]) {
            seenVertex[v] = 1;
        } else {
            int index = v->index;
            const Vector& position = v->position;
            v = mesh.vertices.emplace(mesh.vertices.end(), Vertex(&mesh));
            v->position = position;
            v->index = nVertices + nV++;
            v->referenceIndex = index;
        }

        //  insert edge
        if (!seenEdge[e]) {
            seenEdge[e] = 1;
        } else {
            e = mesh.edges.emplace(mesh.edges.end(), Edge(&mesh));
            e->index = nEdges + nE++;
            e->generatorIndex = generatorIndex;
        }

        // insert halfedge
        HalfEdgeIter newHe = mesh.halfEdges.emplace(mesh.halfEdges.end(), HalfEdge(&mesh));
        newHe->index = nHalfEdges1 + nHe++;
        newHe->setVertex(v);
        newHe->setEdge(e);
        newHe->setFlip(he);
    }

    // insert boundary face
    FaceIter newF = mesh.boundaries.emplace(mesh.boundaries.end(), Face(&mesh));
    newF->setHalfEdge(mesh.halfEdges.begin() + nHalfEdges1);
    newF->index = 0;

    std::vector<HalfEdgeIter> halfEdgesOnGenerators;
    // update connectivity
    // set the prev and next of halfedges on the boundary
    int nHalfEdges2 = (int)mesh.halfEdges.size();
    for (int j = nHalfEdges1; j < nHalfEdges2; j++) {
        int i = j == nHalfEdges1 ? nHalfEdges2 - 1 : j - 1;
        int k = j == nHalfEdges2 - 1 ? nHalfEdges1 : j + 1;

        HalfEdgeIter newPrev = mesh.halfEdges.begin() + i;
        HalfEdgeIter newHe = mesh.halfEdges.begin() + j;
        HalfEdgeIter newNext = mesh.halfEdges.begin() + k;
        HalfEdgeIter flip = newHe->flip();

        // update vertex connectivity
        VertexIter v = newHe->vertex();
        v->setHalfEdge(newHe);
        HalfEdgeIter h = flip->next();
        do {
            h->setVertex(v);

            if (h->edge()->onCut) break;
            h = h->flip()->next();
        } while (true);

        // update edge connectivity
        EdgeIter e = newHe->edge();
        e->setHalfEdge(flip);
        e->onGenerator = true;
        halfEdgesOnGenerators.emplace_back(flip);

        // update face connectivity
        newHe->setFace(newF);

        // update halfedge connectivity
        flip->setFlip(newHe);
        newHe->setPrev(newPrev);
        newHe->setNext(newNext);
        newHe->onBoundary = true;
    }

    // assign halfedge edges and remove cut label from all edges
    for (EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        e->halfEdge()->setEdge(e);
        e->onCut = false;
    }

    // reverse those halfedges inside of the boundary
    std::reverse(halfEdgesOnGenerators.begin(), halfEdgesOnGenerators.end());
    // return the edges on the two generators/boundaries after cutting
    return halfEdgesOnGenerators;
}

void MeshSurgery::cutMesh(Mesh& mesh) {
    int nVertices = (int)mesh.vertices.size();
    int nEdges = (int)mesh.edges.size();
    int nHalfEdges1 = (int)mesh.halfEdges.size();
    int nV = 0;
    int nE = 0;
    int nHe = 0;

    // insert vertices, edges and halfedges
    VertexData<uint8_t> seenVertex(mesh, 0);
    EdgeData<uint8_t> seenEdge(mesh, 0);
    const CutPtrSet& cutBoundary = mesh.cutBoundary();
    for (const WedgeIter& w : cutBoundary) {
        VertexIter v = w->vertex();
        HalfEdgeIter he = w->halfEdge()->next();
        EdgeIter e = he->edge();

        // insert vertex
        if (!seenVertex[v]) {
            seenVertex[v] = 1;
        } else {
            int index = v->index;
            const Vector& position = v->position;
            v = mesh.vertices.emplace(mesh.vertices.end(), Vertex(&mesh));
            v->position = position;
            v->index = nVertices + nV++;
            v->referenceIndex = index;
        }

        //  insert edge
        if (!seenEdge[e]) {
            seenEdge[e] = 1;
        } else {
            // VertexIter e1 = he->vertex();
            // VertexIter e0 = he->flip()->vertex();
            // // flip to print, because wedge halfedge is the flip of he
            // std::cout << "e_"<<e->index << "=(" << e0->index << ", " << e1->index << "), ";
            e = mesh.edges.emplace(mesh.edges.end(), Edge(&mesh));
            e->index = nEdges + nE++;
        }

        // insert halfedge
        HalfEdgeIter newHe = mesh.halfEdges.emplace(mesh.halfEdges.end(), HalfEdge(&mesh));
        newHe->index = nHalfEdges1 + nHe++;
        newHe->setVertex(v);
        newHe->setEdge(e);
        newHe->setFlip(he);
    }

    // insert boundary face
    FaceIter newF = mesh.boundaries.emplace(mesh.boundaries.end(), Face(&mesh));
    newF->setHalfEdge(mesh.halfEdges.begin() + nHalfEdges1);
    newF->index = 0;

    // update connectivity
    // set the prev and next of halfedges on the boundary
    int nHalfEdges2 = (int)mesh.halfEdges.size();
    for (int j = nHalfEdges1; j < nHalfEdges2; j++) {
        int i = j == nHalfEdges1 ? nHalfEdges2 - 1 : j - 1;
        int k = j == nHalfEdges2 - 1 ? nHalfEdges1 : j + 1;

        HalfEdgeIter newPrev = mesh.halfEdges.begin() + i;
        HalfEdgeIter newHe = mesh.halfEdges.begin() + j;
        HalfEdgeIter newNext = mesh.halfEdges.begin() + k;
        HalfEdgeIter flip = newHe->flip();

        // update vertex connectivity
        VertexIter v = newHe->vertex();
        v->setHalfEdge(newHe);
        HalfEdgeIter h = flip->next();
        do {
            h->setVertex(v);

            if (h->edge()->onCut) break;
            h = h->flip()->next();
        } while (true);

        // update edge connectivity
        EdgeIter e = newHe->edge();
        e->setHalfEdge(flip);
        e->onGenerator = true;

        // update face connectivity
        newHe->setFace(newF);

        // update halfedge connectivity
        flip->setFlip(newHe);
        newHe->setPrev(newPrev);
        newHe->setNext(newNext);
        newHe->onBoundary = true;
    }

    // assign halfedge edges and remove cut label from all edges
    for (EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        e->halfEdge()->setEdge(e);
        e->onCut = false;
    }
}

}  // namespace geomp