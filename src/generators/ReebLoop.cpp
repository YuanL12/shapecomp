#include "geomp/generators/ReebLoop.h"
#include "geomp/mesh/Mesh.h"



// Forward declarations for ReebHanTun functions
extern void RandomUniqueDirection(_SimpleMesh &mesh, Vector3 &uniDirection);
extern void CycleLocalOptimization(_SimpleMesh &locMesh, psbmReebGraph &reebgraph, std::vector<int> &OrientTriangles,
                                   std::vector<std::set<int> > &out_v_basis_loops,
                                   std::vector<std::set<int> > &out_h_basis_loops,
                                   const float fScaleRatio);
extern void CycleLocalOptimization_bdry(_SimpleMesh &locMesh, psbmReebGraph &reebgraph, std::vector<int> &OrientTriangles,
                                         std::vector<std::set<int> > &out_v_basis_loops,
                                         std::vector<std::set<int> > &out_h_basis_loops,
                                         std::set<int> extraVertices, const float fScaleRatio);
extern void validateDirectionVector(const _SimpleMesh &mesh, const Vector3 &distinctDirection);

namespace geomp {

void convertMeshToSimpleMesh(const Mesh& mesh, _SimpleMesh& simpleMesh, 
                              std::vector<int>& orientTriangles,
                              std::set<int>& extraVertices,
                              std::vector<Vector3>& meshNormal) {
    // Clear the simple mesh
    simpleMesh.vecVertex.clear();
    simpleMesh.vecEdge.clear();
    simpleMesh.vecTriangle.clear();
    orientTriangles.clear();
    extraVertices.clear();
    meshNormal.clear();
    
    // Convert vertices - maintain index alignment
    simpleMesh.vecVertex.reserve(mesh.vertices.size());
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        _SimpleMeshVertex tmpVer;
        tmpVer.x = static_cast<float>(v->position.x);
        tmpVer.y = static_cast<float>(v->position.y);
        tmpVer.z = static_cast<float>(v->position.z);
        simpleMesh.vecVertex.push_back(tmpVer);
    }
    
    // Convert faces to triangles and build edges
    // Use the same approach as LoadMeshInOFFformat
    std::map<std::pair<int, int>, int, myPairCompare> edgeMapping;
    orientTriangles.reserve(mesh.faces.size() * 3);
    meshNormal.reserve(mesh.faces.size());
    
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isReal() || f->fillsHole) {
            continue; // Skip non-real faces and hole fillers
        }
        
        // Get vertex indices from the face
        HalfEdgeCIter he = f->halfEdge();
        int v0 = he->vertex()->index;
        int v1 = he->next()->vertex()->index;
        int v2 = he->next()->next()->vertex()->index;
        
        // Store oriented triangle
        orientTriangles.push_back(v0);
        orientTriangles.push_back(v1);
        orientTriangles.push_back(v2);
        
        // Compute triangle normal BEFORE sorting (same as LoadMeshInOFFformat)
        // Use original vertex order to match the expected normal direction
        Vector3 leftVec, rightVec;
        leftVec[0] = simpleMesh.vecVertex[v2].x - simpleMesh.vecVertex[v1].x;
        leftVec[1] = simpleMesh.vecVertex[v2].y - simpleMesh.vecVertex[v1].y;
        leftVec[2] = simpleMesh.vecVertex[v2].z - simpleMesh.vecVertex[v1].z;
        
        rightVec[0] = simpleMesh.vecVertex[v0].x - simpleMesh.vecVertex[v1].x;
        rightVec[1] = simpleMesh.vecVertex[v0].y - simpleMesh.vecVertex[v1].y;
        rightVec[2] = simpleMesh.vecVertex[v0].z - simpleMesh.vecVertex[v1].z;
        
        leftVec = leftVec ^ rightVec;  // cross product
        unitize(leftVec);  // normalize
        meshNormal.push_back(leftVec);
        
        // Create triangle with sorted vertices (v0 < v1 < v2)
        _SimpleMeshTriangle tmpTri;
        tmpTri.v0 = v0;
        tmpTri.v1 = v1;
        tmpTri.v2 = v2;
        tmpTri.sortVertices();
        
        // Build edges for this triangle
        _SimpleMeshEdge tmpEdge;
        
        // Edge v0-v1
        std::pair<int, int> edgePair(tmpTri.v0, tmpTri.v1);
        auto mIter = edgeMapping.find(edgePair);
        if (mIter == edgeMapping.end()) {
            tmpEdge.v0 = edgePair.first;
            tmpEdge.v1 = edgePair.second;
            tmpEdge.AdjTri[0] = simpleMesh.vecTriangle.size();
            tmpEdge.AdjTriNum = 1;
            simpleMesh.vecEdge.push_back(tmpEdge);
            tmpTri.e01 = simpleMesh.vecEdge.size() - 1;
            edgeMapping[edgePair] = tmpTri.e01;
        } else {
            simpleMesh.vecEdge[mIter->second].AdjTriNum++;
            simpleMesh.vecEdge[mIter->second].AdjTri[1] = simpleMesh.vecTriangle.size();
            tmpTri.e01 = mIter->second;
        }
        
        // Edge v1-v2
        edgePair = std::make_pair(tmpTri.v1, tmpTri.v2);
        mIter = edgeMapping.find(edgePair);
        if (mIter == edgeMapping.end()) {
            tmpEdge.v0 = edgePair.first;
            tmpEdge.v1 = edgePair.second;
            tmpEdge.AdjTri[0] = simpleMesh.vecTriangle.size();
            tmpEdge.AdjTriNum = 1;
            simpleMesh.vecEdge.push_back(tmpEdge);
            tmpTri.e12 = simpleMesh.vecEdge.size() - 1;
            edgeMapping[edgePair] = tmpTri.e12;
        } else {
            simpleMesh.vecEdge[mIter->second].AdjTriNum++;
            simpleMesh.vecEdge[mIter->second].AdjTri[1] = simpleMesh.vecTriangle.size();
            tmpTri.e12 = mIter->second;
        }
        
        // Edge v0-v2
        edgePair = std::make_pair(tmpTri.v0, tmpTri.v2);
        mIter = edgeMapping.find(edgePair);
        if (mIter == edgeMapping.end()) {
            tmpEdge.v0 = edgePair.first;
            tmpEdge.v1 = edgePair.second;
            tmpEdge.AdjTri[0] = simpleMesh.vecTriangle.size();
            tmpEdge.AdjTriNum = 1;
            simpleMesh.vecEdge.push_back(tmpEdge);
            tmpTri.e02 = simpleMesh.vecEdge.size() - 1;
            edgeMapping[edgePair] = tmpTri.e02;
        } else {
            simpleMesh.vecEdge[mIter->second].AdjTriNum++;
            simpleMesh.vecEdge[mIter->second].AdjTri[1] = simpleMesh.vecTriangle.size();
            tmpTri.e02 = mIter->second;
        }
        
        simpleMesh.vecTriangle.push_back(tmpTri);
    }
    
    // Assign incident edges to vertices
    for (size_t i = 0; i < simpleMesh.vecEdge.size(); i++) {
        simpleMesh.vecVertex[simpleMesh.vecEdge[i].v0].adjEdges.push_back(static_cast<int>(i));
        simpleMesh.vecVertex[simpleMesh.vecEdge[i].v1].adjEdges.push_back(static_cast<int>(i));
    }
    
    // Set the mesh normal pointer (critical for ReebHanTun to work)
    simpleMesh.SetMeshNormalPtr(&meshNormal);
}

ReebHanTunLoops::ReebHanTunLoops(const Mesh& mesh) {
    // Convert geomp::Mesh to _SimpleMesh
    convertMeshToSimpleMesh(mesh, simpleMesh, orientTriangles, extraVertices, meshNormal);
    
    // Check for boundaries (edges with only one adjacent triangle)
    for (size_t i = 0; i < simpleMesh.vecEdge.size(); i++) {
        if (simpleMesh.vecEdge[i].AdjTriNum == 1) {
            throw std::runtime_error("Mesh has boundary. The ReebLoop algorithm only works for closed surfaces.");
        }
    }
    // Check genus
    int EulerCharacteristic = static_cast<int>(simpleMesh.vecVertex.size() + simpleMesh.vecTriangle.size() - simpleMesh.vecEdge.size());
    // For meshes with boundary, we'd need to trace boundary loops properly
    // For now, assume closed mesh (torus should be closed anyway)
    int genus = 1 - EulerCharacteristic / 2;
    if (genus == 0) {
        throw std::runtime_error("Mesh has genus 0 - no generators to compute");
    }
    
    // Set up Reeb graph computation
    Vector3 distinctDirection;
    RandomUniqueDirection(simpleMesh, distinctDirection);
    std::cout << "distinctDirection: " << distinctDirection[0] << " " << distinctDirection[1] << " " << distinctDirection[2] << std::endl;
    
    psbmReebGraph reebGraph;
    reebGraph.ReserveSpaceForEdges(simpleMesh.vecEdge.size());
    
    // Compute scalar field
    std::vector<double> scalarField(simpleMesh.vecVertex.size());
    for (size_t i = 0; i < simpleMesh.vecVertex.size(); i++) {
        scalarField[i] = simpleMesh.vecVertex[i].x * distinctDirection[0] +
                         simpleMesh.vecVertex[i].y * distinctDirection[1] +
                         simpleMesh.vecVertex[i].z * distinctDirection[2];
    }
    
    reebGraph.SetHeightDirection(distinctDirection);
    reebGraph.AssignData(&simpleMesh, scalarField.data());
    reebGraph.scalarDir = 1;
    
    // Compute Reeb graph
    reebGraph.ComputeReebGraph();
    
    // Compute cycles
    reebGraph.ComputingCycle_max_tree();
    reebGraph.compute_path_on_mesh_for_each_simplified_arc();
    reebGraph.EmbedCycleAsEdgePathOnMesh();
    reebGraph.LinkNumberMatrixComputing();
    
    // Get handle and tunnel loops
    const float fEnlargeFactor = 10000.f;
    if (extraVertices.empty()) {
        CycleLocalOptimization(simpleMesh, reebGraph, orientTriangles, v_basis_loops, h_basis_loops,
                               1.f / fEnlargeFactor);
    } else {
        CycleLocalOptimization_bdry(simpleMesh, reebGraph, orientTriangles, v_basis_loops, h_basis_loops, 
                                     extraVertices, 1.f / fEnlargeFactor);
    }
}

ReebHanTunLoops::ReebHanTunLoops(const Mesh& mesh, const Vector3& distinctDirection) {
    // Convert geomp::Mesh to _SimpleMesh
    convertMeshToSimpleMesh(mesh, simpleMesh, orientTriangles, extraVertices, meshNormal);
    
    // Check for boundaries (edges with only one adjacent triangle)
    for (size_t i = 0; i < simpleMesh.vecEdge.size(); i++) {
        if (simpleMesh.vecEdge[i].AdjTriNum == 1) {
            throw std::runtime_error("Mesh has boundary. The ReebLoop algorithm only works for closed surfaces.");
        }
    }
    // Check genus
    int EulerCharacteristic = static_cast<int>(simpleMesh.vecVertex.size() + simpleMesh.vecTriangle.size() - simpleMesh.vecEdge.size());
    // For meshes with boundary, we'd need to trace boundary loops properly
    // For now, assume closed mesh (torus should be closed anyway)
    int genus = 1 - EulerCharacteristic / 2;
    if (genus == 0) {
        throw std::runtime_error("Mesh has genus 0 - no generators to compute");
    }
    
    // Normalize the direction vector to ensure it's exactly normalized
    // (even if normalized in Python, there may be floating point precision differences)
    Vector3 normalizedDirection = distinctDirection;
    unitize(normalizedDirection);
    
    // Valid the direction vector (expects normalized vector)
    validateDirectionVector(simpleMesh, normalizedDirection);
    
    // Set up Reeb graph computation
    psbmReebGraph reebGraph;
    reebGraph.ReserveSpaceForEdges(simpleMesh.vecEdge.size());
    
    // Compute scalar field using normalized direction
    std::vector<double> scalarField(simpleMesh.vecVertex.size());
    for (size_t i = 0; i < simpleMesh.vecVertex.size(); i++) {
        scalarField[i] = simpleMesh.vecVertex[i].x * normalizedDirection[0] +
                         simpleMesh.vecVertex[i].y * normalizedDirection[1] +
                         simpleMesh.vecVertex[i].z * normalizedDirection[2];
    }
    
    reebGraph.SetHeightDirection(normalizedDirection);
    reebGraph.AssignData(&simpleMesh, scalarField.data());
    reebGraph.scalarDir = 1;
    
    // Compute Reeb graph
    reebGraph.ComputeReebGraph();
    
    // Compute cycles
    reebGraph.ComputingCycle_max_tree();
    reebGraph.compute_path_on_mesh_for_each_simplified_arc();
    reebGraph.EmbedCycleAsEdgePathOnMesh();
    reebGraph.LinkNumberMatrixComputing();
    
    // Get handle and tunnel loops
    const float fEnlargeFactor = 10000.f;
    if (extraVertices.empty()) {
        CycleLocalOptimization(simpleMesh, reebGraph, orientTriangles, v_basis_loops, h_basis_loops,
                               1.f / fEnlargeFactor);
    } else {
        CycleLocalOptimization_bdry(simpleMesh, reebGraph, orientTriangles, v_basis_loops, h_basis_loops, 
                                     extraVertices, 1.f / fEnlargeFactor);
    }
}

void ReebHanTunLoops::labelReebGraphGeneratorsOnMesh(Mesh& mesh) {   
    // Label handle generator edges (generatorIndex = 1)
    for (size_t i = 0; i < h_basis_loops.size(); i++) {
        for (int edgeIdx : h_basis_loops[i]) {
            if (edgeIdx < 0 || edgeIdx >= static_cast<int>(simpleMesh.vecEdge.size())) {
                continue;
            }
            int v0 = simpleMesh.vecEdge[edgeIdx].v0;
            int v1 = simpleMesh.vecEdge[edgeIdx].v1;
            
            // Find the corresponding edge in geomp::Mesh
            EdgeIter e = mesh.getEdge(v0, v1);
            e->onCut = true;
            e->generatorIndex = 1;
        }
    }
    
    // Label tunnel generator edges (generatorIndex = 2)
    for (size_t i = 0; i < v_basis_loops.size(); i++) {
        for (int edgeIdx : v_basis_loops[i]) {
            if (edgeIdx < 0 || edgeIdx >= static_cast<int>(simpleMesh.vecEdge.size())) {
                continue;
            }
            int v0 = simpleMesh.vecEdge[edgeIdx].v0;
            int v1 = simpleMesh.vecEdge[edgeIdx].v1;
            
            // Find the corresponding edge in geomp::Mesh
            EdgeIter e = mesh.getEdge(v0, v1);
            e->onCut = true;
            e->generatorIndex = 2;
        }
    }
}

std::pair<std::vector<int>, std::vector<int>> ReebHanTunLoops::findReebGraphGenerators() {
    // find the vertices on the handle generator
    std::set<int> handleVertices;
    for (size_t i = 0; i < h_basis_loops.size(); i++) {
        for (int edgeIdx : h_basis_loops[i]) {
            if (edgeIdx < 0 || edgeIdx >= static_cast<int>(simpleMesh.vecEdge.size())) {
                continue;
            }
            int v0 = simpleMesh.vecEdge[edgeIdx].v0;
            int v1 = simpleMesh.vecEdge[edgeIdx].v1;
            
            handleVertices.insert(v0);
            handleVertices.insert(v1);
        }
    }
    
    // find the vertices on the tunnel generator
    std::set<int> tunnelVertices;
    for (size_t i = 0; i < v_basis_loops.size(); i++) {
        for (int edgeIdx : v_basis_loops[i]) {
            if (edgeIdx < 0 || edgeIdx >= static_cast<int>(simpleMesh.vecEdge.size())) {
                continue;
            }
            int v0 = simpleMesh.vecEdge[edgeIdx].v0;
            int v1 = simpleMesh.vecEdge[edgeIdx].v1;
            
            // Find the corresponding edge in geomp::Mesh
            tunnelVertices.insert(v0);
            tunnelVertices.insert(v1);
        }
    }
    return std::make_pair(std::vector<int>(handleVertices.begin(), handleVertices.end()), std::vector<int>(tunnelVertices.begin(), tunnelVertices.end()));
}
}  // namespace geomp
