#include "geomp/mesh/MeshIO.h"


#ifdef USE_USD
	#include <pxr/usd/usd/prim.h>
	#include <pxr/usd/usd/stage.h>
	#include <pxr/usd/usd/primRange.h>
	#include <pxr/usd/usdGeom/mesh.h>
	#include <pxr/usd/usdGeom/primvar.h>
	#include <pxr/usd/usdGeom/primvarsAPI.h>
#endif

namespace geomp {

bool MeshIO::readOBJ(const std::string& fileName, PolygonSoup& soup,
					 std::vector<std::pair<int, int>>& uncuttableEdges,
					 std::string& error)
{
	std::ifstream in(fileName.c_str());
	// treat a string like a stream 
	// (e.g., parsing a string as if it were input).
	std::istringstream buffer;

	if (!in.is_open()) {
		error = "Unable to open file";
		return false;

	} else {
		buffer.str(std::string(std::istreambuf_iterator<char>(in),
							   std::istreambuf_iterator<char>()));
	}

	std::string line;
	int nVertices = 0;
	bool seenFace = false;
	soup.positions.clear();
	soup.indices.clear();
	uncuttableEdges.clear();

	// extract a line from the file and stores it in the `line` string
	while (std::getline(buffer, line)) {
		std::istringstream ss(line); 
		std::string token; // same as a word
		ss >> token; // read the first word from the line

		if (token == "v") {
			double x, y, z;
			ss >> x >> y >> z;
			soup.positions.emplace_back(Vector(x, y, z));

			if (seenFace) {
				nVertices = 0;
				seenFace = false;
			}
			nVertices++;

		} else if (token == "f") {
			seenFace = true;
			int vertexCount = 0;
			int rootIndex = 0;
			int prevIndex = 0;

			// read the next word/token from the line stream, token can be 1 or 1/2
			while (ss >> token) { 
				std::istringstream indexStream(token); // create a new stream
				std::string indexString;

				// Split `x/y` using '/'
				if (std::getline(indexStream, indexString, '/')) {
					int index = std::stoi(indexString); 
					// get the index of the vertex
					if (index < 0) index = nVertices + index;
					else index -= 1; // 1-based to 0-based

					if (vertexCount == 0) rootIndex = index;

					vertexCount++;
					// Triangulation of Faces with More Than 3 Vertices
					if (vertexCount > 3) {
						// triangulate polygon if vertexCount > 3
						soup.indices.emplace_back(rootIndex);
						soup.indices.emplace_back(prevIndex);
						soup.indices.emplace_back(index);

						// tag edge as uncuttable
						int i = rootIndex;
						int j = prevIndex;
						std::pair<int, int> edge(i, j);
						uncuttableEdges.emplace_back(edge);

					} else {
						soup.indices.emplace_back(index);
					}

					prevIndex = index;
				}
			}
		}
	}

	in.close();
	return true;
}


bool MeshIO::writeMeshOBJ(const std::string& fileName, const Mesh& mesh)
{
	std::ofstream out(fileName.c_str());
	if (!out.is_open()) {
		return false;
	}

	// write vertices
	for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
		out << "v " << v->position.x << " " << v->position.y << " " << v->position.z << std::endl;
	}

	// write faces
	for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
		out << "f";
		HalfEdgeCIter he = f->halfEdge();
		do {
			out << " " << he->vertex()->index + 1;
			he = he->next();
		} while (he != f->halfEdge());
		out << std::endl;
	}

	out.close();
	return true;
}

bool MeshIO::writeMeshOBJ_uv(
	const VertexData<std::pair<double, double>>& uv_positions, 
	const std::string& fileName, 
	const Mesh& mesh)
{
	std::ofstream out(fileName.c_str());
	if (!out.is_open()) {
		return false;
	}

	// write uv coordinates
	for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
		out << "vt " << uv_positions[v].first << " " << uv_positions[v].second << std::endl;
	}

	// write faces
	for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
		out << "f";
		HalfEdgeCIter he = f->halfEdge();
		do {
			VertexCIter v = he->vertex();
			out << " " << v->index << "/" << v->referenceIndex;
			he = he->next();
		} while (he != f->halfEdge());
		out << std::endl;
	}

	out.close();
	return true;

}


#ifdef USE_USD
bool MeshIO::readUSD(const std::string& fileName, PolygonSoup& soup,
					 std::vector<std::pair<int, int>>& uncuttableEdges,
					 std::string& error)
{
	// open USD file
	auto stage = pxr::UsdStage::Open(fileName);
	if (!stage) {
		error = "Unable to open file";
		return false;
	}

	// iterate over all the meshes in the USD file and get
	// their vertex positions and face indices
	for (auto prim : stage->Traverse()) {
		if (prim.IsA<pxr::UsdGeomMesh>()) {
			pxr::UsdGeomMesh mesh(prim);

			// get vertex positions and face indices
			pxr::VtVec3fArray points;
			mesh.GetPointsAttr().Get(&points);
			pxr::VtIntArray faceVertexCounts;
			mesh.GetFaceVertexCountsAttr().Get(&faceVertexCounts);
			pxr::VtIntArray faceVertexIndices;
			mesh.GetFaceVertexIndicesAttr().Get(&faceVertexIndices);
			if (points.size() == 0 || faceVertexCounts.size() == 0 || faceVertexIndices.size() == 0) {
				continue;
			}

			// add vertex positions to soup
			int globalVertexCount = (int)soup.positions.size();
			for (int i = 0; i < (int)points.size(); i++) {
				soup.positions.emplace_back(Vector(points[i][0], points[i][1], points[i][2]));
			}

			// add face indices to soup
			int index = 0;
			for (int i = 0; i < (int)faceVertexCounts.size(); i++) {
				int vertexCount = faceVertexCounts[i];
				int rootIndex = globalVertexCount + faceVertexIndices[index];
				int prevIndex = globalVertexCount + faceVertexIndices[index];

				for (int j = 0; j < vertexCount; j++) {
					int currentIndex = globalVertexCount + faceVertexIndices[index + j];

					if (j == 0) rootIndex = currentIndex;
					if (j > 2) {
						// triangulate polygon if vertexCount > 3
						soup.indices.emplace_back(rootIndex);
						soup.indices.emplace_back(prevIndex);
						soup.indices.emplace_back(currentIndex);

						// tag edge as uncuttable
						int i = rootIndex;
						int j = prevIndex;
						std::pair<int, int> edge(i, j);
						uncuttableEdges.emplace_back(edge);

					} else {
						soup.indices.emplace_back(currentIndex);
					}

					prevIndex = currentIndex;
				}

				index += vertexCount;
			}
		}
	}

	return true;
}
#endif

void MeshIO::separateComponents(const PolygonSoup& soup, int nComponents,
								const std::vector<int>& faceComponent,
								const std::vector<uint8_t>& isCuttableModelEdge,
								std::vector<PolygonSoup>& soups,
								std::vector<std::vector<uint8_t>>& isCuttableSoupEdge,
								std::vector<std::pair<int, int>>& modelToMeshMap,
								std::vector<std::vector<int>>& meshToModelMap)
{
	int nIndices = (int)soup.indices.size();
	int nVertices = (int)soup.positions.size();
	meshToModelMap.resize(nComponents);
	modelToMeshMap.resize(nVertices);

	if (nComponents == 1) {
		soups.emplace_back(std::move(soup));
		isCuttableSoupEdge.emplace_back(std::move(isCuttableModelEdge));
		meshToModelMap[0].resize(nVertices);

		for (int i = 0; i < nVertices; i++) {
			meshToModelMap[0][i] = i;
			modelToMeshMap[i] = std::make_pair(0, i);
		}

	} else {
		// create soups
		soups.resize(nComponents);
		isCuttableSoupEdge.resize(nComponents);
		std::vector<int> vertexIndexMap(nVertices, -1);

		for (int I = 0; I < nIndices; I += 3) {
			int c = faceComponent[I];

			for (int J = 0; J < 3; J++) {
				int i = soup.indices[I + J];

				// insert vertex if it hasn't been seen
				if (vertexIndexMap[i] == -1) {
					int index = (int)soups[c].positions.size();
					vertexIndexMap[i] = index;
					soups[c].positions.emplace_back(soup.positions[i]);
					meshToModelMap[c].emplace_back(i);
					modelToMeshMap[i] = std::make_pair(c, index);
				}

				// add index
				soups[c].indices.emplace_back(vertexIndexMap[i]);
			}
		}

		// construct tables
		for (int c = 0; c < nComponents; c++) {
			soups[c].vertexAdjacency.construct(soups[c].positions.size(), soups[c].indices);
			isCuttableSoupEdge[c].resize(soups[c].vertexAdjacency.getEdgeCount(), 1);
		}

		// mark cuttable edges for each soup
		for (int I = 0; I < nIndices; I += 3) {
			int c = faceComponent[I];

			for (int J = 0; J < 3; J++) {
				int K = (J + 1) % 3;
				int i = soup.indices[I + J];
				int j = soup.indices[I + K];

				int eIndex = soup.vertexAdjacency.getEdgeIndex(i, j);

				// add edge if uncuttable
				if (!isCuttableModelEdge[eIndex]) {
					int ii = vertexIndexMap[i];
					int jj = vertexIndexMap[j];

					int eIndexSoup = soups[c].vertexAdjacency.getEdgeIndex(ii, jj);
					isCuttableSoupEdge[c][eIndexSoup] = 0;
				}
			}
		}
	}
}

void MeshIO::preallocateElements(const PolygonSoup& soup, Mesh& mesh)
{
	// clear arrays
	mesh.vertices.clear();
	mesh.edges.clear();
	mesh.faces.clear();
	mesh.corners.clear();
	mesh.halfEdges.clear();
	mesh.boundaries.clear();

	// reserve space (reserving extra for hole filling and generators)
	int nVertices = (int)soup.positions.size();
	int nEdges = soup.vertexAdjacency.getEdgeCount();
	int nFaces = (int)soup.indices.size()/3;
	int nCorners = 3*nFaces;
	int nHalfedges = 2*nEdges;

	mesh.vertices.reserve(2*nVertices);
	mesh.edges.reserve(2*nEdges);
	mesh.faces.reserve(2*nFaces);
	mesh.corners.reserve(2*nCorners);
	mesh.halfEdges.reserve(2*nHalfedges);
}

bool MeshIO::hasIsolatedVertices(const Mesh& mesh)
{
	for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
		if (v->isIsolated()) {
			return true;
		}
	}

	return false;
}

bool MeshIO::hasNonManifoldVertices(const Mesh& mesh)
{
	VertexData<int> adjacentFaces(mesh, 0);
	for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
		HalfEdgeCIter h = f->halfEdge();
		do {
			adjacentFaces[h->vertex()]++;

			h = h->next();
		} while (h != f->halfEdge());
	}

	for (BoundaryCIter b = mesh.boundaries.begin(); b != mesh.boundaries.end(); b++) {
		HalfEdgeCIter h = b->halfEdge();
		do {
			adjacentFaces[h->vertex()]++;

			h = h->next();
		} while (h != b->halfEdge());
	}

	for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
		if (adjacentFaces[v] != v->degree()) {
			return true;
		}
	}

	return false;
}

bool MeshIO::buildMesh(const PolygonSoup& soup,
					   const std::vector<uint8_t>& isCuttableEdge,
					   Mesh& mesh, std::string& error)
{
	// preallocate elements
	preallocateElements(soup, mesh);

	// create and insert vertices
	int nVertices = (int)soup.positions.size();
	std::vector<VertexIter> indexToVertex(nVertices);
	for (int i = 0; i < nVertices; i++) {
		VertexIter v = mesh.vertices.emplace(mesh.vertices.end(), Vertex(&mesh));
		v->position = soup.positions[i];
		v->index = (int)mesh.vertices.size() - 1;
		indexToVertex[i] = v;
	}

	// create and insert halfedges, edges and "real" faces
	int nEdges = soup.vertexAdjacency.getEdgeCount();
	std::vector<int> edgeCount(nEdges, 0);
	std::vector<HalfEdgeIter> existingHalfEdges(nEdges);
	std::vector<int> hasFlipEdge(2*nEdges, 0);

	for (int I = 0; I < (int)soup.indices.size(); I += 3) {
		// create new face
		FaceIter f = mesh.faces.emplace(mesh.faces.end(), Face(&mesh));
		f->index = (int)mesh.faces.size() - 1;

		// create a halfedge for each edge of the newly created face
		std::vector<HalfEdgeIter> halfEdges(3);
		for (int J = 0; J < 3; J++) {
			halfEdges[J] = mesh.halfEdges.emplace(mesh.halfEdges.end(), HalfEdge(&mesh));
			halfEdges[J]->index = (int)mesh.halfEdges.size() - 1;
		}

		// initialize the newly created halfedges
		// order is (0 1), (1 2), (2 0)
		for (int J = 0; J < 3; J++) {
			// current halfedge goes from vertex i to vertex j
			int K = (J + 1) % 3;
			int i = soup.indices[I + J];
			int j = soup.indices[I + K];

			// set the current halfedge's attributes
			HalfEdgeIter h = halfEdges[J];
			h->setNext(halfEdges[K]);
			h->setPrev(halfEdges[(J + 3 - 1) % 3]);
			h->onBoundary = false;
			hasFlipEdge[h->index] = 0;

			// point the new halfedge and vertex i to each other
			VertexIter v = indexToVertex[i];
			h->setVertex(v);
			v->setHalfEdge(h);

			// point the new halfedge and face to each other
			h->setFace(f);
			f->setHalfEdge(h);

			int eIndex = soup.vertexAdjacency.getEdgeIndex(i, j);
			if (edgeCount[eIndex] > 0) {
				// if a halfedge between vertex i and j has been created in the past,
				// then it is the flip halfedge of the current halfedge
				HalfEdgeIter flip = existingHalfEdges[eIndex];
				h->setFlip(flip);
				flip->setFlip(h);
				h->setEdge(flip->edge());

				hasFlipEdge[h->index] = 1;
				hasFlipEdge[flip->index] = 1;
				edgeCount[eIndex]++;

			} else {
				// create an edge and set its halfedge
				EdgeIter e = mesh.edges.emplace(mesh.edges.end(), Edge(&mesh));
				e->index = (int)mesh.edges.size() - 1;
				h->setEdge(e);
				e->setHalfEdge(h);
				e->isCuttable = isCuttableEdge[eIndex];

				// record the newly created edge and halfedge from vertex i to j
				existingHalfEdges[eIndex] = h;
				edgeCount[eIndex] = 1;
			}

			// check for non-manifold edges
			if (edgeCount[eIndex] > 2) {
				mesh.status = Mesh::ErrorCode::nonManifoldEdges;
				error = "Mesh has non-manifold edges.";
				return false;
			}
		}
	}

	// create and insert boundary halfedges and "imaginary" faces for boundary cycles
	// also create and insert corners
	for (HalfEdgeIter h = mesh.halfEdges.begin(); h != mesh.halfEdges.end(); h++) {
		// if a halfedge has no flip halfedge, create a new face and
		// link it the corresponding boundary cycle
		if (!hasFlipEdge[h->index]) {
			// create new face
			FaceIter f = mesh.boundaries.emplace(mesh.boundaries.end(), Face(&mesh));
			f->index = (int)mesh.boundaries.size() - 1;

			// walk along boundary cycle
			std::vector<HalfEdgeIter> boundaryCycle;
			boundaryCycle.reserve(2*mesh.edges.size() - mesh.halfEdges.size());
			HalfEdgeIter he = h;
			do {
				// create a new halfedge
				HalfEdgeIter bH = mesh.halfEdges.emplace(mesh.halfEdges.end(), HalfEdge(&mesh));
				bH->index = (int)mesh.halfEdges.size() - 1;
				boundaryCycle.emplace_back(bH);

				// grab the next halfedge along the boundary that does not have a
				// flip halfedge
				HalfEdgeIter nextHe = he->next();
				while (hasFlipEdge[nextHe->index]) {
					nextHe = nextHe->flip()->next();
				}

				// set the current halfedge's attributes
				bH->setVertex(nextHe->vertex());
				bH->setEdge(he->edge());
				bH->onBoundary = true;

				// point the new halfedge and face to each other
				bH->setFace(f);
				f->setHalfEdge(bH);

				// point the new halfedge and he to each other
				bH->setFlip(he);
				he->setFlip(bH);

				// continue walk
				he = nextHe;
			} while (he != h);

			// link the cycle of boundary halfedges together
			int n = (int)boundaryCycle.size();
			for (int i = 0; i < n; i++) {
				boundaryCycle[i]->setNext(boundaryCycle[(i + n - 1) % n]); // boundary halfedges are linked in clockwise order
				boundaryCycle[i]->setPrev(boundaryCycle[(i + 1) % n]);
				hasFlipEdge[boundaryCycle[i]->index] = 1;
				hasFlipEdge[boundaryCycle[i]->flip()->index] = 1;
			}
		}

		// point the newly created corner and its halfedge to each other
		if (!h->onBoundary) {
			CornerIter c = mesh.corners.emplace(mesh.corners.end(), Corner(&mesh));
			c->index = (int)mesh.corners.size() - 1;
			c->setHalfEdge(h);
			h->setCorner(c);
		}
	}

	// check if mesh has isolated vertices
	if (hasIsolatedVertices(mesh)) {
		mesh.status = Mesh::ErrorCode::isolatedVertices;
		error = "Mesh has isolated vertices.";
		return false;
	}

	// check if mesh has non-manifold vertices
	if (hasNonManifoldVertices(mesh)) {
		mesh.status = Mesh::ErrorCode::nonManifoldVertices;
		error = "Mesh has non manifold vertices.";
		return false;
	}

	return true;
}


void MeshIO::BarycenterDivideMesh(
	const Mesh& mesh, Mesh& newMesh)
{
	// // assume no boundaries
	// if (mesh.boundaries.size() > 0) 
    //     throw std::runtime_error("BarycenterDivideMesh: only support mesh without boundaries");

	// reserve enough space for the new mesh
	size_t nVertices_old = mesh.vertices.size();
	size_t nEdges_old = mesh.edges.size();
	size_t nHalfEdges_old = mesh.halfEdges.size();
	size_t nFaces_old = mesh.faces.size();
	
	// each edge contributes 1 new vertex, so does each face
	size_t nVertices = nVertices_old + nEdges_old + nFaces_old;
	// each edge becomes 2, and each face gives 6 new edges
	size_t nEdges = 2*nEdges_old + 6*nFaces_old;
	size_t nHalfEdges = 2*nEdges;
	// each face gives 6 new faces
	size_t nFaces = 6*nFaces_old;
	// each face still has 3 corners
	size_t nCorners = 3*nFaces;
	newMesh.vertices.reserve(2*nVertices);
	newMesh.edges.reserve(2*nEdges);
	newMesh.faces.reserve(2*nFaces);
	newMesh.corners.reserve(2*nCorners);
	newMesh.halfEdges.reserve(2*nHalfEdges);
	newMesh.boundaries.reserve(2*mesh.boundaries.size());
	

	// index of new items in the new mesh
	int nV = 0;  
	int nE = 0;
	int nHe = 0;
	int nF = 0;

	// insert vertices from the old mesh
	for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
		VertexIter newV = newMesh.vertices.emplace(newMesh.vertices.end(), Vertex(&newMesh));
		newV->position = v->position;
		newV->index = nV++;
	}

	// store the map from an edge in the old mesh to its mid vertex in the new mesh
	EdgeData<VertexIter> midVertex(mesh);
	// map a halfedge to a pair of halfedges that are created by subdividing the edge
	// e.g. (u, v) -> (u, w) and (w, v)
	HalfEdgeData<std::pair<HalfEdgeIter, HalfEdgeIter>> halfEdgesToNewHEs(mesh);

	// iterate over edges of the old mesh
	// insert mid vertices, edges and halfedges from the old mesh
	for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
		HalfEdgeCIter he = e->halfEdge();
		VertexCIter u = he->vertex();
		VertexCIter v = he->flip()->vertex();

		// insert mid vertex w
		VertexIter w = newMesh.vertices.emplace(newMesh.vertices.end(), Vertex(&newMesh));
		w->position = 0.5*(u->position + v->position);
		w->index = nV++;

		// store the mid vertex per edge
		midVertex[e] = w;

		// insert new edge u->w and its two halfedges
		EdgeIter e1 = newMesh.edges.emplace(newMesh.edges.end(), Edge(&newMesh));
		e1->index = nE++;
		HalfEdgeIter uw = newMesh.halfEdges.emplace(newMesh.halfEdges.end(), HalfEdge(&newMesh));
		uw->index = nHe++;uw->setVertex(u);uw->setEdge(e1);
		HalfEdgeIter wu = newMesh.halfEdges.emplace(newMesh.halfEdges.end(), HalfEdge(&newMesh));
		wu->index = nHe++;wu->setVertex(w);wu->setEdge(e1);
		uw->setFlip(wu);wu->setFlip(uw);

		// insert new edge w->v and its two halfedges
		EdgeIter e2 = newMesh.edges.emplace(newMesh.edges.end(), Edge(&newMesh));
		e2->index = nE++;
		HalfEdgeIter wv = newMesh.halfEdges.emplace(newMesh.halfEdges.end(), HalfEdge(&newMesh));
		wv->index = nHe++;wv->setVertex(w);wv->setEdge(e2);
		HalfEdgeIter vw = newMesh.halfEdges.emplace(newMesh.halfEdges.end(), HalfEdge(&newMesh));
		vw->index = nHe++;vw->setVertex(v);vw->setEdge(e2);
		wv->setFlip(vw);vw->setFlip(wv);

		// update the connectivity
		e1->setHalfEdge(uw);e2->setHalfEdge(wv);
		w->setHalfEdge(wu);

		// store the new halfedges
		halfEdgesToNewHEs[he] = std::make_pair(uw, wv);
		halfEdgesToNewHEs[he->flip()] = std::make_pair(vw, wu);
	}

	// iterate over faces of the old mesh
	// insert center vertices, 6 faces
	for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
		// create a new vertex at the center of the face
		VertexIter d = newMesh.vertices.emplace(newMesh.vertices.end(), Vertex(&newMesh));
		d->position = f->barycenter(); d->index = nV++;

		// insert edges first 
		HalfEdgeCIter he = f->halfEdge(); // he = (u, v)
		HalfEdgeCIter he1 = he->next(); // he1 = (v, w)
		HalfEdgeCIter he2 = he1->next(); // he2 = (w, u)

		// get the midpoints of the edges and u, v, w
		VertexIter a = midVertex[he->edge()];
		VertexIter b = midVertex[he1->edge()];
		VertexIter c = midVertex[he2->edge()];
		VertexIter u = he->vertex(); VertexIter v = he1->vertex(); VertexIter w = he2->vertex();
		
		// insert 6 edges from d to a, b, c, u, v, w and their halfedges
		// store halfedges starting from d
		std::vector<HalfEdgeIter> sixNewHalfEdges; 
		sixNewHalfEdges.reserve(6);
		std::vector<VertexIter> candidates = {u, v, w, a, b, c};
		for (VertexIter candidate: candidates) {
			// insert edge e=(d, candidate) and its two halfedges
			EdgeIter e = newMesh.edges.emplace(newMesh.edges.end(), Edge(&newMesh));
			e->index = nE++;

			// halfedge d to candidate
			HalfEdgeIter dCand = newMesh.halfEdges.emplace(newMesh.halfEdges.end(), HalfEdge(&newMesh));
			dCand->index = nHe++;dCand->setVertex(d);dCand->setEdge(e);
			
			// halfedge candidate to d
			HalfEdgeIter candd = newMesh.halfEdges.emplace(newMesh.halfEdges.end(), HalfEdge(&newMesh));
			candd->index = nHe++;candd->setVertex(candidate);candd->setEdge(e);

			// update the connectivity
			dCand->setFlip(candd);candd->setFlip(dCand);
			e->setHalfEdge(dCand);
			
			// set one of the halfedges for center
			if (candidate == u) {
				d->setHalfEdge(dCand);
			}

			// store the halfedge for inserting new faces later
			sixNewHalfEdges.emplace_back(dCand);
		}

		// insert 6 new faces
		HalfEdgeIter du = sixNewHalfEdges[0];
		HalfEdgeIter dv = sixNewHalfEdges[1];
		HalfEdgeIter dw = sixNewHalfEdges[2];
		HalfEdgeIter da = sixNewHalfEdges[3];
		HalfEdgeIter db = sixNewHalfEdges[4];
		HalfEdgeIter dc = sixNewHalfEdges[5];

		HalfEdgeIter ua = halfEdgesToNewHEs[he].first;
		HalfEdgeIter av = halfEdgesToNewHEs[he].second;
		HalfEdgeIter vb = halfEdgesToNewHEs[he1].first;
		HalfEdgeIter bw = halfEdgesToNewHEs[he1].second;
		HalfEdgeIter wc = halfEdgesToNewHEs[he2].first;
		HalfEdgeIter cu = halfEdgesToNewHEs[he2].second;

		// f1 = (u, a, d)
		FaceIter f1 = newMesh.faces.emplace(newMesh.faces.end(), Face(&newMesh));
		f1->index = nF++; setNewFace(f1, ua, da->flip(), du);

		// f2 = (a, v, d)
		FaceIter f2 = newMesh.faces.emplace(newMesh.faces.end(), Face(&newMesh));
		f2->index = nF++; setNewFace(f2, av, dv->flip(), da);

		// f3 = (v, b, d)
		FaceIter f3 = newMesh.faces.emplace(newMesh.faces.end(), Face(&newMesh));
		f3->index = nF++; setNewFace(f3, vb, db->flip(), dv);

		// f4 = (b, w, d)
		FaceIter f4 = newMesh.faces.emplace(newMesh.faces.end(), Face(&newMesh));
		f4->index = nF++; setNewFace(f4, bw, dw->flip(), db);

		// f5 = (w, c, d)
		FaceIter f5 = newMesh.faces.emplace(newMesh.faces.end(), Face(&newMesh));
		f5->index = nF++; setNewFace(f5, wc, dc->flip(), dw);

		// f6 = (c, u, d)
		FaceIter f6 = newMesh.faces.emplace(newMesh.faces.end(), Face(&newMesh));
		f6->index = nF++; setNewFace(f6, cu, du->flip(), dc);
	}

	// create and insert corners in the new mesh
	for (HalfEdgeIter h = newMesh.halfEdges.begin(); h != newMesh.halfEdges.end(); h++) {
		// point the newly created corner and its halfedge to each other
		CornerIter c = newMesh.corners.emplace(newMesh.corners.end(), Corner(&newMesh));
		c->index = (int)newMesh.corners.size() - 1;
		c->setHalfEdge(h);
		h->setCorner(c);
	}
}

void MeshIO::setNewFace(
FaceIter f, 
HalfEdgeIter h1,
HalfEdgeIter h2,
HalfEdgeIter h3)
{
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

/*
 * @brief Baricenter divide a mesh, while label edges on Generator 1 and Generator 2.
 * 				Assume the input mesh has label edges on Generator 1 as onCut. 
 * @param mesh input mesh
 * @param newMesh the new mesh to be modified
 * @param dualCycle the second generator stored as faces in the old mesh
 * @return a pair of edges on Generator 1 and Generator 2

*/
std::pair<std::vector<EdgeCIter>, std::vector<EdgeCIter>> MeshIO::BarycenterDivideLabelMesh(const Mesh& mesh, Mesh& newMesh, const std::vector<FaceCIter>& dualCycle){
	// assume no boundaries
	if (mesh.boundaries.size() > 0) 
        throw std::runtime_error("BarycenterDivideMesh: only support mesh without boundaries");

	// reserve enough space for the new mesh
	size_t nVertices_old = mesh.vertices.size();
	size_t nEdges_old = mesh.edges.size();
	size_t nHalfEdges_old = mesh.halfEdges.size();
	size_t nFaces_old = mesh.faces.size();
	
	// each edge contributes 1 new vertex, so does each face
	size_t nVertices = nVertices_old + nEdges_old + nFaces_old;
	// each edge becomes 2, and each face gives 6 new edges
	size_t nEdges = 2*nEdges_old + 6*nFaces_old;
	size_t nHalfEdges = 2*nEdges;
	// each face gives 6 new faces
	size_t nFaces = 6*nFaces_old;
	// each face still has 3 corners
	size_t nCorners = 3*nFaces;
	newMesh.vertices.reserve(2*nVertices);
	newMesh.edges.reserve(2*nEdges);
	newMesh.faces.reserve(2*nFaces);
	newMesh.corners.reserve(2*nCorners);
	newMesh.halfEdges.reserve(2*nHalfEdges);
	newMesh.boundaries.reserve(2*mesh.boundaries.size());
	

	// index of new items in the new mesh
	int nV = 0;  
	int nE = 0;
	int nHe = 0;
	int nF = 0;

	// store the map from an old vertex to a new vertex
	VertexData<VertexIter> o2NV(mesh);
	// insert vertices from the old mesh
	for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
		VertexIter newV = newMesh.vertices.emplace(newMesh.vertices.end(), Vertex(&newMesh));
		newV->position = v->position;
		newV->index = nV++;
		o2NV[v] = newV;
	}

	// store the map from an edge in the old mesh to its mid vertex in the new mesh
	EdgeData<VertexIter> midVertex(mesh);
	// map a halfedge to a pair of halfedges that are created by subdividing the edge
	// e.g. (u, v) -> (u, w) and (w, v)
	HalfEdgeData<std::pair<HalfEdgeIter, HalfEdgeIter>> halfEdgesToNewHEs(mesh);

	// store the first generator edges
	std::vector<EdgeCIter> firstGeneratorEdges;

	// iterate over edges of the old mesh
	// insert mid vertices, edges and halfedges from the old mesh
	for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
		HalfEdgeCIter he = e->halfEdge();
		VertexCIter u = he->vertex(); // from old mesh, need to use o2NV[u]
		VertexCIter v = he->flip()->vertex();

		// insert mid vertex w
		VertexIter w = newMesh.vertices.emplace(newMesh.vertices.end(), Vertex(&newMesh));
		w->position = 0.5*(u->position + v->position);
		w->index = nV++;

		// store the mid vertex per edge
		midVertex[e] = w;

		// insert new edge u->w and its two halfedges
		EdgeIter e1 = newMesh.edges.emplace(newMesh.edges.end(), Edge(&newMesh));
		e1->isCuttable = true; e1->index = nE++;
		HalfEdgeIter uw = newMesh.halfEdges.emplace(newMesh.halfEdges.end(), HalfEdge(&newMesh));
		uw->index = nHe++;uw->setVertex(o2NV[u]);uw->setEdge(e1);
		HalfEdgeIter wu = newMesh.halfEdges.emplace(newMesh.halfEdges.end(), HalfEdge(&newMesh));
		wu->index = nHe++;wu->setVertex(w);wu->setEdge(e1);
		uw->setFlip(wu);wu->setFlip(uw);

		// insert new edge w->v and its two halfedges
		EdgeIter e2 = newMesh.edges.emplace(newMesh.edges.end(), Edge(&newMesh));
		e2->isCuttable = true; e2->index = nE++;
		HalfEdgeIter wv = newMesh.halfEdges.emplace(newMesh.halfEdges.end(), HalfEdge(&newMesh));
		wv->index = nHe++;wv->setVertex(w);wv->setEdge(e2);
		HalfEdgeIter vw = newMesh.halfEdges.emplace(newMesh.halfEdges.end(), HalfEdge(&newMesh));
		vw->index = nHe++;vw->setVertex(o2NV[v]);vw->setEdge(e2);
		wv->setFlip(vw);vw->setFlip(wv);

		// update the connectivity
		e1->setHalfEdge(uw);e2->setHalfEdge(wv);
		w->setHalfEdge(wu);
		o2NV[u]->setHalfEdge(uw);
		o2NV[v]->setHalfEdge(vw);

		// store the new halfedges
		halfEdgesToNewHEs[he] = std::make_pair(uw, wv);
		halfEdgesToNewHEs[he->flip()] = std::make_pair(vw, wu);

		// if edge e is on the cut, then e1 and e2 are on the cut
		if (e->onCut) {
			// std::cout << "edge " << *e << " onCut" << std::endl;
			e1->onCut = true; e1->generatorIndex = 1;
			e2->onCut = true; e2->generatorIndex = 1;
			// add the two new edges to the first generator
			firstGeneratorEdges.emplace_back(e1);
			firstGeneratorEdges.emplace_back(e2);
		}
	}

	// iterate over faces of the old mesh
	// insert center vertices, 6 faces
	for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
		// create a new vertex at the center of the face
		VertexIter d = newMesh.vertices.emplace(newMesh.vertices.end(), Vertex(&newMesh));
		d->position = f->barycenter(); d->index = nV++;

		// insert edges first 
		HalfEdgeCIter he = f->halfEdge(); // he = (u, v)
		HalfEdgeCIter he1 = he->next(); // he1 = (v, w)
		HalfEdgeCIter he2 = he1->next(); // he2 = (w, u)

		// get the midpoints of the edges and u, v, w
		VertexIter a = midVertex[he->edge()];
		VertexIter b = midVertex[he1->edge()];
		VertexIter c = midVertex[he2->edge()];
		VertexIter u = o2NV[he->vertex()]; 
		VertexIter v = o2NV[he1->vertex()]; 
		VertexIter w = o2NV[he2->vertex()];
		
		// insert 6 edges from d to a, b, c, u, v, w and their halfedges
		// store halfedges starting from d
		std::vector<HalfEdgeIter> sixNewHalfEdges; 
		sixNewHalfEdges.reserve(6);
		std::vector<VertexIter> candidates = {u, v, w, a, b, c};
		for (VertexIter candidate: candidates) {
			// insert edge e=(d, candidate) and its two halfedges
			EdgeIter e = newMesh.edges.emplace(newMesh.edges.end(), Edge(&newMesh));
			e->index = nE++; e->isCuttable = true;

			// halfedge d to candidate
			HalfEdgeIter dCand = newMesh.halfEdges.emplace(newMesh.halfEdges.end(), HalfEdge(&newMesh));
			dCand->index = nHe++;dCand->setVertex(d);dCand->setEdge(e);
			
			// halfedge candidate to d
			HalfEdgeIter candd = newMesh.halfEdges.emplace(newMesh.halfEdges.end(), HalfEdge(&newMesh));
			candd->index = nHe++;candd->setVertex(candidate);candd->setEdge(e);

			// update the connectivity
			dCand->setFlip(candd);candd->setFlip(dCand);
			e->setHalfEdge(dCand);
			
			// set one of the halfedges for center
			if (candidate == u) {
				d->setHalfEdge(dCand);
			}

			// store the halfedge for inserting new faces later
			sixNewHalfEdges.emplace_back(dCand);
		}

		// insert 6 new faces
		HalfEdgeIter du = sixNewHalfEdges[0];
		HalfEdgeIter dv = sixNewHalfEdges[1];
		HalfEdgeIter dw = sixNewHalfEdges[2];
		HalfEdgeIter da = sixNewHalfEdges[3];
		HalfEdgeIter db = sixNewHalfEdges[4];
		HalfEdgeIter dc = sixNewHalfEdges[5];

		HalfEdgeIter ua = halfEdgesToNewHEs[he].first;
		HalfEdgeIter av = halfEdgesToNewHEs[he].second;
		HalfEdgeIter vb = halfEdgesToNewHEs[he1].first;
		HalfEdgeIter bw = halfEdgesToNewHEs[he1].second;
		HalfEdgeIter wc = halfEdgesToNewHEs[he2].first;
		HalfEdgeIter cu = halfEdgesToNewHEs[he2].second;

		// f1 = (u, a, d)
		FaceIter f1 = newMesh.faces.emplace(newMesh.faces.end(), Face(&newMesh));
		f1->index = nF++; setNewFace(f1, ua, da->flip(), du);

		// f2 = (a, v, d)
		FaceIter f2 = newMesh.faces.emplace(newMesh.faces.end(), Face(&newMesh));
		f2->index = nF++; setNewFace(f2, av, dv->flip(), da);

		// f3 = (v, b, d)
		FaceIter f3 = newMesh.faces.emplace(newMesh.faces.end(), Face(&newMesh));
		f3->index = nF++; setNewFace(f3, vb, db->flip(), dv);

		// f4 = (b, w, d)
		FaceIter f4 = newMesh.faces.emplace(newMesh.faces.end(), Face(&newMesh));
		f4->index = nF++; setNewFace(f4, bw, dw->flip(), db);

		// f5 = (w, c, d)
		FaceIter f5 = newMesh.faces.emplace(newMesh.faces.end(), Face(&newMesh));
		f5->index = nF++; setNewFace(f5, wc, dc->flip(), dw);

		// f6 = (c, u, d)
		FaceIter f6 = newMesh.faces.emplace(newMesh.faces.end(), Face(&newMesh));
		f6->index = nF++; setNewFace(f6, cu, du->flip(), dc);
	}


	// transform faces in the dual cycle to edges
	// and store the edges in the second generator
	std::vector<EdgeCIter> secondGeneratorEdges;
	for (size_t i = 0; i < dualCycle.size(); i++) {
		FaceCIter f = dualCycle[i];
		FaceCIter g = dualCycle[(i + 1) % dualCycle.size()];

		// find the halfedge that is shared by faces f and g
		EdgeCIter e = sharedEdgeByFaces(f, g);

		// label and two edges as onCut and add them to the second generator
		EdgeIter ef = halfEdgesToNewHEs[e->halfEdge()].first->next()->edge();
		ef->onCut = true; ef->generatorIndex = 2;
		EdgeIter eg = halfEdgesToNewHEs[e->halfEdge()->flip()].first->next()->edge();
		eg->onCut = true; eg->generatorIndex = 2;
		secondGeneratorEdges.emplace_back(ef);
		secondGeneratorEdges.emplace_back(eg);
	}


	// create and insert corners in the new mesh
	for (HalfEdgeIter h = newMesh.halfEdges.begin(); h != newMesh.halfEdges.end(); h++) {
		// point the newly created corner and its halfedge to each other
		CornerIter c = newMesh.corners.emplace(newMesh.corners.end(), Corner(&newMesh));
		c->index = (int)newMesh.corners.size() - 1;
		c->setHalfEdge(h);
		h->setCorner(c);
	}

	// find the shared root of the two generators by counting the apperance of each vertex
	std::unordered_map<int, int> vertexCount;
	for (EdgeCIter e: firstGeneratorEdges) {
		vertexCount[e->halfEdge()->vertex()->index]++;
		vertexCount[e->halfEdge()->flip()->vertex()->index]++;
	}
	for (EdgeCIter e: secondGeneratorEdges) {
		vertexCount[e->halfEdge()->vertex()->index]++;
		vertexCount[e->halfEdge()->flip()->vertex()->index]++;
	}
	// check if all vertices should appear twice except the root appears 4 times
	int rootIndex = -1;
	for (auto& pair: vertexCount) {
		if (pair.second == 4) {
			rootIndex = pair.first;
			// std::cout << "find root with index: " << rootIndex << std::endl;
		} else if (pair.second != 2) {
			throw std::runtime_error("BarycenterDivideLabelMesh: the two generators appear more than once");
		}
	}

	// return the pair of two edge generators
	return std::make_pair(firstGeneratorEdges, secondGeneratorEdges);
}


void MeshIO::normalize_radius(Model& model, bool rescale_to_unit_radius)
{
	// compute center of mass
	Vector cm;
	int nVertices = 0;
	for (int i = 0; i < model.size(); i++) {
		for (VertexCIter v = model[i].vertices.begin(); v != model[i].vertices.end(); v++) {
			cm += v->position;
			nVertices++;
		}
	}
	cm /= nVertices;

	// translate to origin and determine radius
	double radius = 1e-8;
	for (int i = 0; i < model.size(); i++) {
		for (VertexIter v = model[i].vertices.begin(); v != model[i].vertices.end(); v++) {
			v->position -= cm;
			radius = std::max(radius, v->position.norm());
		}
	}

	// rescale to unit radius
	if (rescale_to_unit_radius) {
		for (int i = 0; i < model.size(); i++) {
			for (VertexIter v = model[i].vertices.begin(); v != model[i].vertices.end(); v++) {
				v->position /= radius;
			}

			model[i].radius = radius;
			model[i].cm = cm;
		}
	}
}

void MeshIO::normalize_area(Model& model)
{
	// compute center of mass
	Vector cm;
	int nVertices = 0;
	for (int i = 0; i < model.size(); i++) {
		for (VertexCIter v = model[i].vertices.begin(); v != model[i].vertices.end(); v++) {
			cm += v->position;
			nVertices++;
		}
	}
	cm /= nVertices;

	// translate to origin
	for (int i = 0; i < model.size(); i++) {
		for (VertexIter v = model[i].vertices.begin(); v != model[i].vertices.end(); v++) {
			v->position -= cm;
		}
	}

	// rescale to unit area
	for (int i = 0; i < model.size(); i++) {
		double total_area = 0.0;
		for (FaceCIter f = model[i].faces.begin(); f != model[i].faces.end(); f++) {
			total_area += area(f);
		}
	
		double scale = std::sqrt(1/total_area);
		for (VertexIter v = model[i].vertices.begin(); v != model[i].vertices.end(); v++) {
			v->position *= scale;
		}
	}
}


bool MeshIO::buildModel(const std::vector<std::pair<int, int>>& uncuttableEdges,
						PolygonSoup& soup, Model& model, std::string& error, 
						std::optional<bool> normalize)
{
	// remove isolated vertices
	bool removedIsolatedVertices = soup.removeIsolatedVertices();

	// construct vertex-edge and edge-face adjacency maps
	soup.vertexAdjacency.construct(soup.positions.size(), soup.indices);
	soup.edgeFaceAdjacency.construct(soup.vertexAdjacency, soup.indices);

	// split non-manifold vertices
	std::vector<uint8_t> isCuttableModelEdge(soup.vertexAdjacency.getEdgeCount(), 1);
	bool didSplitNonManifoldVertices = soup.splitNonManifoldVertices();
	if (removedIsolatedVertices || didSplitNonManifoldVertices) {
		// for simplicity, allow uncuttable edges to be cut during flattening since
		// connectivity is being changed in any case when vertices are removed or split
		isCuttableModelEdge.resize(soup.vertexAdjacency.getEdgeCount(), 1);

	} else {
		for (int i = 0; i < (int)uncuttableEdges.size(); i++) {
			int eIndex = soup.vertexAdjacency.getEdgeIndex(uncuttableEdges[i].first, uncuttableEdges[i].second);
			isCuttableModelEdge[eIndex] = 0;
		}
	}

	// orient faces and assign components to them
	std::vector<int> faceComponent;
	int nComponents = soup.orientFacesAndAssignComponents(faceComponent);

	// clear edge-face adjacency map as its no longer needed
	soup.edgeFaceAdjacency.clear();

	// separate model into components
	std::vector<PolygonSoup> soups;
	std::vector<std::vector<uint8_t>> isCuttableSoupEdge;
	separateComponents(soup, nComponents, faceComponent, isCuttableModelEdge,
					   soups, isCuttableSoupEdge, model.modelToMeshMap,
					   model.meshToModelMap);

	// build halfedge meshes
	model.meshes.resize(soups.size());
	for (int i = 0; i < (int)soups.size(); i++) {
		if (!buildMesh(soups[i], isCuttableSoupEdge[i], model[i], error)) {
			return false;
		}
	}

	// normalize model, default is False
	if (normalize.has_value() && normalize.value()) {
		// std::cout << "run normalize"<< std::endl;
		normalize_area(model);
	}
	return true;
}

bool MeshIO::read(const std::string& fileName, Model& model, std::string& error, std::optional<bool> normalize)
{
	// read polygon soup from obj file
	PolygonSoup soup;

	// uncuttable Edges are those from traingulation of faces(general polygons) with more than 3 vertices
	std::vector<std::pair<int, int>> uncuttableEdges;

	// Start to read the file
	if (fileName.find(".obj") != std::string::npos) {
		if (!readOBJ(fileName, soup, uncuttableEdges, error)) {
			return false;
		}
#ifdef USE_USD
	} else if (fileName.find(".usd") != std::string::npos) {
		if (!readUSD(fileName, soup, uncuttableEdges, error)) {
			return false;
		}
#endif
	} else {
		error = "Unable to open file";
		return false;
	}

	// build model
	if (!buildModel(uncuttableEdges, soup, model, error, normalize)) {
		return false;
	}

	return true;
}


// void MeshIO::collectModelUvs(Model& model, bool normalizeUvs,
// 							 const std::vector<uint8_t>& isSurfaceMappedToSphere,
// 							 const std::vector<Vector>& originalUvIslandCenters,
// 							 const std::vector<Vector>& newUvIslandCenters,
// 							 const std::vector<uint8_t>& isUvIslandFlipped,
// 							 const Vector& modelMinBounds,
// 							 const Vector& modelMaxBounds,
// 							 std::vector<Vector>& positions,
// 							 std::vector<Vector>& uvs,
// 							 std::vector<int>& vIndices,
// 							 std::vector<int>& uvIndices,
// 							 std::vector<int>& indicesOffset)
// {
// 	// initialize
// 	int nVertices = model.nVertices();
// 	int nFaces = 0;
// 	for (int i = 0; i < model.size(); i++) {
// 		nFaces += (int)model[i].faces.size();
// 	}

// 	positions.clear(); positions.reserve(nVertices);
// 	uvs.clear(); uvs.reserve(nVertices);
// 	vIndices.clear(); vIndices.reserve(3*nFaces);
// 	uvIndices.clear(); uvIndices.reserve(3*nFaces);
// 	indicesOffset.clear(); indicesOffset.reserve(nFaces + 1);
// 	indicesOffset.emplace_back(0);

// 	// collect vertex positions
// 	for (int i = 0; i < nVertices; i++) {
// 		std::pair<int, int> vData = model.localVertexIndex(i);
// 		const Mesh& mesh = model[vData.first];
// 		VertexCIter v = mesh.vertices.begin() + vData.second;

// 		if (isSurfaceMappedToSphere[vData.first] == 1) {
// 			const Vector& uv = v->wedge()->uv;
// 			positions.emplace_back(uv);

// 		} else {
// 			Vector p = v->position*mesh.radius + mesh.cm;
// 			positions.emplace_back(p);
// 		}
// 	}

// 	// collect UVs and indices
// 	int nV = 0;
// 	int nUvs = 0;
// 	for (int i = 0; i < model.size(); i++) {
// 		// compute uv radius and shift
// 		Vector minExtent(modelMinBounds.x, modelMinBounds.y);
// 		double dx = modelMaxBounds.x - minExtent.x;
// 		double dy = modelMaxBounds.y - minExtent.y;
// 		double extent = std::max(dx, dy);
// 		minExtent.x -= (extent - dx)/2.0;
// 		minExtent.y -= (extent - dy)/2.0;

// 		// compute sphere radius if component has been mapped to a sphere
// 		double sphereRadius = 1.0;
// 		if (isSurfaceMappedToSphere[i] == 1) {
// 			for (WedgeCIter w = model[i].wedges().begin(); w != model[i].wedges().end(); w++) {
// 				sphereRadius = std::max(w->uv.norm(), sphereRadius);
// 			}
// 		}

// 		// collect interior UVs
// 		int uvCount = 0;
// 		HalfEdgeData<int> uvIndexMap(model[i]);

// 		// compute the ratio of the surface areas of the mesh and flattened mesh
// 		double lengthRatio = std::sqrt(model[i].areaRatio());

// 		for (VertexCIter v = model[i].vertices.begin(); v != model[i].vertices.end(); v++) {
// 			if (!v->onBoundary()) {
// 				Vector uv = v->wedge()->uv;
// 				if (isSurfaceMappedToSphere[i] == 1) {
// 					uv /= sphereRadius;
// 					uv.x = 0.5 + atan2(uv.z, uv.x)/(2*M_PI);
// 					uv.y = 0.5 - asin(uv.y)/M_PI;

// 				} else {
// 					uv *= model[i].radius*lengthRatio;
// 				}

// 				uv -= originalUvIslandCenters[i];
// 				if (isUvIslandFlipped[i] == 1) uv = Vector(-uv.y, uv.x);
// 				uv += newUvIslandCenters[i];
// 				uv -= minExtent;
// 				if (normalizeUvs) uv /= extent;
// 				uvs.emplace_back(uv);

// 				HalfEdgeCIter he = v->halfEdge();
// 				do {
// 					uvIndexMap[he->next()] = uvCount;

// 					he = he->flip()->next();
// 				} while (he != v->halfEdge());

// 				uvCount++;
// 			}
// 		}

// 		// collect boundary UVs
// 		for (WedgeCIter w: model[i].cutBoundary()) {
// 			Vector uv = w->uv;
// 			if (isSurfaceMappedToSphere[i] == 1) {
// 				uv /= sphereRadius;
// 				uv.x = 0.5 + atan2(uv.z, uv.x)/(2*M_PI);
// 				uv.y = 0.5 - asin(uv.y)/M_PI;

// 			} else {
// 				uv *= model[i].radius*lengthRatio;
// 			}

// 			uv -= originalUvIslandCenters[i];
// 			if (isUvIslandFlipped[i] == 1) uv = Vector(-uv.y, uv.x);
// 			uv += newUvIslandCenters[i];
// 			uv -= minExtent;
// 			if (normalizeUvs) uv /= extent;
// 			uvs.emplace_back(uv);

// 			HalfEdgeCIter he = w->halfEdge()->prev();
// 			do {
// 				uvIndexMap[he->next()] = uvCount;

// 				if (he->edge()->onCut) break;
// 				he = he->flip()->next();
// 			} while (!he->onBoundary);

// 			uvCount++;
// 		}

// 		// collect indices
// 		int uncuttableEdges = 0;
// 		for (FaceCIter f = model[i].faces.begin(); f != model[i].faces.end(); f++) {
// 			if (!f->fillsHole) {
// 				if (uncuttableEdges > 0) {
// 					uncuttableEdges--;
// 					continue;
// 				}

// 				HalfEdgeCIter he = f->halfEdge()->next();
// 				while (!he->edge()->isCuttable) he = he->next();
// 				HalfEdgeCIter fhe = he;
// 				std::unordered_set<int> seenUncuttableEdges;

// 				do {
// 					VertexCIter v = he->vertex();
// 					int vIndex = v->referenceIndex == -1 ? v->index : v->referenceIndex;
// 					vIndices.emplace_back(model.globalVertexIndex(i, vIndex));
// 					uvIndices.emplace_back(nUvs + uvIndexMap[he->next()]);
// 					nV++;

// 					he = he->next();
// 					while (!he->edge()->isCuttable) {
// 						seenUncuttableEdges.emplace(he->edge()->index);
// 						he = he->flip()->next();
// 					}

// 				} while (he != fhe);

// 				indicesOffset.emplace_back(nV);
// 				uncuttableEdges = (int)seenUncuttableEdges.size();
// 			}
// 		}

// 		nUvs += uvCount;
// 	}
// }

void writeString(std::ofstream& out, const std::string& str)
{
	out.write(str.c_str(), str.size());
}

bool MeshIO::writeOBJ(const std::string& fileName, bool writeOnlyUvs,
					  const std::vector<Vector>& positions,
					  const std::vector<Vector>& uvs,
					  const std::vector<int>& vIndices,
					  const std::vector<int>& uvIndices,
					  const std::vector<int>& indicesOffset)
{
	std::ofstream out(fileName.c_str());
	if (!out.is_open()) {
		return false;
	}

	if (!writeOnlyUvs) {
		for (int i = 0; i < (int)positions.size(); i++) {
			const Vector& p = positions[i];
			writeString(out, "v " + std::to_string(p.x) + " " +
									std::to_string(p.y) + " " +
									std::to_string(p.z) + "\n");
		}
	}

	for (int i = 0; i < (int)uvs.size(); i++) {
		const Vector& uv = uvs[i];
		if (writeOnlyUvs) {
			writeString(out, "v " + std::to_string(uv.x) + " " +
									std::to_string(uv.y) + " 0.0\n");

		} else {
			writeString(out, "vt " + std::to_string(uv.x) + " " +
									 std::to_string(uv.y) + "\n");
		}
	}

	for (int i = 0; i < (int)indicesOffset.size() - 1; i++) {
		writeString(out, "f");

		for (int j = indicesOffset[i]; j < indicesOffset[i + 1]; j++) {
			if (writeOnlyUvs) {
				writeString(out, " " + std::to_string(uvIndices[j] + 1));

			} else {
				writeString(out, " " + std::to_string(vIndices[j] + 1) + "/" +
									   std::to_string(uvIndices[j] + 1));
			}
		}

		writeString(out, "\n");
	}

	out.close();
	return true;
}


bool MeshIO::readEigenMatrices(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces,
                              PolygonSoup& soup, std::vector<std::pair<int, int>>& uncuttableEdges,
                              std::string& error)
{
    // Check input dimensions
    if (vertices.cols() > 3) {
		int n_cols = vertices.cols();
		// print the vertices matrix
		std::cout << "Vertices matrix: " << std::endl;
		std::cout << vertices << std::endl;
		
        error = "# of columns of Vertices matrix is " + std::to_string(n_cols) + ", expected 2 or 3";
        return false;
    }

    // Clear existing data
    soup.positions.clear();
    soup.indices.clear();
    uncuttableEdges.clear();

    // Add vertices
    for (int i = 0; i < vertices.rows(); i++) {
        soup.positions.emplace_back(Vector(vertices(i, 0), vertices(i, 1), vertices(i, 2)));
    }

    // Add faces
    for (int i = 0; i < faces.rows(); i++) {
        // Check if indices are valid
        if (faces(i, 0) < 0 || faces(i, 0) >= vertices.rows() ||
            faces(i, 1) < 0 || faces(i, 1) >= vertices.rows() ||
            faces(i, 2) < 0 || faces(i, 2) >= vertices.rows()) {
            error = "Invalid face indices";
            return false;
        }

        // Add face indices
        soup.indices.emplace_back(faces(i, 0));
        soup.indices.emplace_back(faces(i, 1));
        soup.indices.emplace_back(faces(i, 2));
    }

    return true;
}

} // namespace geomp
