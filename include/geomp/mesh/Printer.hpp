// printer functions to print halfedge
#pragma once
#include <iostream>

inline void printHalfEdge(HalfEdgeIter he) {
    std::cout << "HalfEdge: " << he->index << " ";
    std::cout << "Vertex: " << he->vertex()->index << " ";
    std::cout << "Edge: " << he->edge()->index << " ";
    std::cout << "Face: " << he->face()->index << " ";
    std::cout << "Flip: " << he->flip()->index << " ";
    std::cout << "Next: " << he->next()->index << " ";
    std::cout << "Prev: " << he->prev()->index << std::endl;
}

inline void printFaces(const Mesh& mesh){
	std::cout << "Print updated faces" << std::endl;
	for (FaceIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
		// get 3 vertcies of the face
		HalfEdgeIter he = f->halfEdge();
		HalfEdgeIter he2 = he->next();
		HalfEdgeIter he3 = he2->next();
		VertexIter v0 = he->vertex();
		VertexIter v1 = he->next()->vertex();
		VertexIter v2 = he->next()->next()->vertex();
		std::cout << "f_" << f->index << "=("<< v0->index << ", " << v1->index << ", " << v2->index << ")";
		// the face has an edge onGenerator, print its two halfedges
		std::cout << ", half egdes on Generators are ";
		if (he->edge()->onGenerator) {
			std::cout << "(" << he->vertex()->index << "->" << he->flip()->vertex()->index << "), ";
		}
		if (he2->edge()->onGenerator) {
			std::cout << "(" << he2->vertex()->index << "->" << he2->flip()->vertex()->index << "), ";
		}
		if (he3->edge()->onGenerator) {
			std::cout << "(" << he3->vertex()->index << "->" << he3->flip()->vertex()->index << "), ";
		}
		std::cout << "\n";
	}
	std::cout << std::endl;
}