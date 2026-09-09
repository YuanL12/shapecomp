/*
(c) 2012 Fengtao Fan
*/
#include "psbmReebGraph.h"
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <string>
#include "SimpleMesh.h"
#include "FilesOutputForOptimalCycles.h"


#include <time.h>
#include <sstream>
#include <boost/progress.hpp>
#include <boost/program_options.hpp>

//
bool CheckBoundaries(_SimpleMesh &inMesh, std::vector<std::vector<std::pair<int, int> > > &boundries) {
    bool ret = false;
    std::vector<bool> edge_flag(inMesh.vecEdge.size(), false);
    for (unsigned int i = 0; i < inMesh.vecEdge.size(); i++) {
        if (!edge_flag[i] && inMesh.vecEdge[i].AdjTriNum == 1) {// this is a boundary edge
            ret = true;
            std::vector<std::pair<int, int> > bdrEdges;
            //
            int terminate_vertex = inMesh.vecEdge[i].v0;
            int current_vertex = inMesh.vecEdge[i].v1;
            int current_edge = i;
            int loop_counter = 0;
            while (current_vertex != terminate_vertex) {
                bdrEdges.push_back(std::pair<int, int>(current_edge, inMesh.vecEdge[current_edge].AdjTri[0]));
                edge_flag[current_edge] = true;
                //
                unsigned int evid = 0;
                for (; evid < inMesh.vecVertex[current_vertex].adjEdges.size(); evid++) {
                    int rot_edge_id = inMesh.vecVertex[current_vertex].adjEdges[evid];
                    if (rot_edge_id != current_edge && inMesh.vecEdge[rot_edge_id].AdjTriNum == 1)
                        break;
                }
                if (evid < inMesh.vecVertex[current_vertex].adjEdges.size()) {
                    current_edge = inMesh.vecVertex[current_vertex].adjEdges[evid];
                    current_vertex = inMesh.vecEdge[current_edge].v0 + inMesh.vecEdge[current_edge].v1 - current_vertex;
                } else {
                    std::cout << "READING MESH MODEL ERROR : " << std::endl;
                    std::cout << "--- CAN NOT FIND ADJACENT BOUDNARY EDGE !" << std::endl;
                    std::cout << "--- PLEASE CHECK THE MESH MODEL" << std::endl;
                    exit(0);
                }
                loop_counter++;
                if (loop_counter > inMesh.vecEdge.size()) {
                    std::cout << "READING MESH MODEL ERROR : " << std::endl;
                    std::cout << "--- PLEASE CHECK THE MESH MODEL" << std::endl;
                    exit(0);
                }
            }
            edge_flag[current_edge] = true;
            bdrEdges.push_back(std::pair<int, int>(current_edge, inMesh.vecEdge[current_edge].AdjTri[0]));
            //
            boundries.push_back(bdrEdges);
        }
        edge_flag[i] = true;
    }
    return ret;
}

void FindEdgeOrientation(_SimpleMesh &inMesh, const int fid, std::vector<int> &OrientTriangles,
                         std::pair<int, int> &EdgeOrientation) {
    int start_index = 3 * fid;
    std::vector<int> TriOrientation(OrientTriangles.begin() + start_index, OrientTriangles.begin() + start_index + 3);
    TriOrientation.push_back(TriOrientation.front());
    //
    for (unsigned int i = 0; i < 3; i++) {
        if ((EdgeOrientation.first == TriOrientation[i] && EdgeOrientation.second == TriOrientation[i + 1]) ||
            (EdgeOrientation.second == TriOrientation[i] && EdgeOrientation.first == TriOrientation[i + 1])) {
            EdgeOrientation.first = TriOrientation[i + 1];
            EdgeOrientation.second = TriOrientation[i];
        }
    }
    return;
}

void CloseHoles(_SimpleMesh &inMesh, std::vector<std::vector<std::pair<int, int> > > &meshBoundaries,
                std::vector<Vector3> &inMeshNormal,
                std::vector<int> &OrientTriangles, std::set<int> &extraVertices) {
    for (unsigned int b = 0; b < meshBoundaries.size(); b++) {
        Vector3 centroid(0.0, 0.0, 0.0);
        int curEdge = 0;
        int nextEdge = 0;
        int curVertex = 0;
        for (unsigned int i = 0; i < meshBoundaries[b].size(); i++) {// the edges on the boundary is ordered
            curEdge = meshBoundaries[b][i].first;
            nextEdge = meshBoundaries[b][0].first;
            if (i < meshBoundaries[b].size() - 1)
                nextEdge = meshBoundaries[b][i + 1].first;
            // shared vertex
            curVertex = inMesh.vecEdge[curEdge].v0;
            if (curVertex != inMesh.vecEdge[nextEdge].v0 &&
                curVertex != inMesh.vecEdge[nextEdge].v1)
                curVertex = inMesh.vecEdge[curEdge].v1;
            //
            centroid = centroid + Vector3(inMesh.vecVertex[curVertex].x,
                                          inMesh.vecVertex[curVertex].y,
                                          inMesh.vecVertex[curVertex].z);
            //
        }
        centroid = centroid / meshBoundaries[b].size();
        //
        // add vertex
        int centroid_index = inMesh.vecVertex.size();
        extraVertices.insert(centroid_index);
        _SimpleMeshVertex tmpVer;
        tmpVer.x = centroid[0];
        tmpVer.y = centroid[1];
        tmpVer.z = centroid[2];
        //
        inMesh.vecVertex.push_back(tmpVer);
        //
        // add edges and triangles
        std::map<int, int> edge_mapping;
        for (unsigned int i = 0; i < meshBoundaries[b].size(); i++) {
            // Find orientation of this edge
            curEdge = meshBoundaries[b][i].first;
            std::pair<int, int> EdgeOrientation(inMesh.vecEdge[curEdge].v0, inMesh.vecEdge[curEdge].v1);
            FindEdgeOrientation(inMesh, meshBoundaries[b][i].second, OrientTriangles, EdgeOrientation);
            //
            OrientTriangles.push_back(EdgeOrientation.first);
            OrientTriangles.push_back(EdgeOrientation.second);
            OrientTriangles.push_back(centroid_index);
            //
            _SimpleMeshTriangle tmpTri;
            _SimpleMeshEdge tmpEdge;
            //
            tmpTri.v0 = EdgeOrientation.first;
            tmpTri.v1 = EdgeOrientation.second;
            tmpTri.v2 = centroid_index;
            //
            //
            Vector3 leftVec, rightVec;
            leftVec[0] = inMesh.vecVertex[tmpTri.v2].x - inMesh.vecVertex[tmpTri.v1].x;
            leftVec[1] = inMesh.vecVertex[tmpTri.v2].y - inMesh.vecVertex[tmpTri.v1].y;
            leftVec[2] = inMesh.vecVertex[tmpTri.v2].z - inMesh.vecVertex[tmpTri.v1].z;

            rightVec[0] = inMesh.vecVertex[tmpTri.v0].x - inMesh.vecVertex[tmpTri.v1].x;
            rightVec[1] = inMesh.vecVertex[tmpTri.v0].y - inMesh.vecVertex[tmpTri.v1].y;
            rightVec[2] = inMesh.vecVertex[tmpTri.v0].z - inMesh.vecVertex[tmpTri.v1].z;
            leftVec = leftVec ^ rightVec;
            unitize(leftVec);
            //leftVec = leftVec / norm(leftVec);
            inMeshNormal.push_back(leftVec);
            //
            tmpTri.sortVertices();
            // [v0 v1] existing edge
            inMesh.vecEdge[curEdge].AdjTriNum++;
            inMesh.vecEdge[curEdge].AdjTri[1] = inMesh.vecTriangle.size();
            //
            tmpTri.e01 = curEdge;
            //
            std::map<int, int>::iterator mIter;
            mIter = edge_mapping.find(tmpTri.v1);
            if (mIter == edge_mapping.end()) {// new edge
                tmpEdge.v0 = tmpTri.v1;
                tmpEdge.v1 = tmpTri.v2;
                tmpEdge.AdjTri[0] = inMesh.vecTriangle.size();
                tmpEdge.AdjTriNum = 1;
                //
                inMesh.vecEdge.push_back(tmpEdge);
                tmpTri.e12 = inMesh.vecEdge.size() - 1;
                edge_mapping[tmpTri.v1] = tmpTri.e12;
                //
                inMesh.vecVertex[tmpEdge.v0].adjEdges.push_back(tmpTri.e12);
                inMesh.vecVertex[tmpEdge.v1].adjEdges.push_back(tmpTri.e12);
            } else {// existed edge
                inMesh.vecEdge[mIter->second].AdjTriNum++;
                inMesh.vecEdge[mIter->second].AdjTri[1] = inMesh.vecTriangle.size();
                //
                tmpTri.e12 = mIter->second;
            }
            mIter = edge_mapping.find(tmpTri.v0);
            if (mIter == edge_mapping.end()) {
                // new edge
                tmpEdge.v0 = tmpTri.v0;
                tmpEdge.v1 = tmpTri.v2;
                tmpEdge.AdjTri[0] = inMesh.vecTriangle.size();
                tmpEdge.AdjTriNum = 1;
                //
                inMesh.vecEdge.push_back(tmpEdge);
                tmpTri.e02 = inMesh.vecEdge.size() - 1;
                edge_mapping[tmpTri.v0] = tmpTri.e02;
                //
                inMesh.vecVertex[tmpEdge.v0].adjEdges.push_back(tmpTri.e02);
                inMesh.vecVertex[tmpEdge.v1].adjEdges.push_back(tmpTri.e02);
            } else {// existed already
                inMesh.vecEdge[mIter->second].AdjTriNum++;
                inMesh.vecEdge[mIter->second].AdjTri[1] = inMesh.vecTriangle.size();
                //
                tmpTri.e02 = mIter->second;
            }
            //
            inMesh.vecTriangle.push_back(tmpTri);
        }
    }
    return;
}

int LoadData(_SimpleMesh &mesh, std::vector<Vector3> &meshNormal,
             const char *mesh_file_name, double &BoundingBoxRadius, std::vector<int> &tris,
             std::set<int> &extraVertices, int &orgTriangleSize,
             const float fEnlargeFactor) {
    int genus = 0;
    _SimpleMeshVertex minBd;
    _SimpleMeshVertex maxBd;
    mesh.LoadMeshInOFFformat(minBd, maxBd, meshNormal, mesh_file_name, tris,
                             fEnlargeFactor);//"E:\\RProgramming\\Models\\torus.off");// "D:\\MeshModels\\OFF-models\\HighGenusCubeHC1.off");//
    mesh.SetMeshNormalPtr(&meshNormal);
    // triangles in TRIS are in the same order as triangles in mesh.vecTriangle;
    //
    orgTriangleSize = mesh.vecTriangle.size();
    //
    std::vector<std::vector<std::pair<int, int> > > meshBoundaries;
    if (CheckBoundaries(mesh, meshBoundaries)) {// it is a mesh with bondary
        int EulerCharacteristic = mesh.vecVertex.size() + mesh.vecTriangle.size() - mesh.vecEdge.size();
        genus = 1 - (EulerCharacteristic + meshBoundaries.size()) / 2;
        if (genus) {
            CloseHoles(mesh, meshBoundaries, meshNormal, tris, extraVertices);
        }
    } else {// it is a closed mesh
        int EulerCharacteristic = mesh.vecVertex.size() + mesh.vecTriangle.size() - mesh.vecEdge.size();
        genus = 1 - EulerCharacteristic / 2;
    }
    if (!genus) {
        std::cout << "NOTHING IS COMPUTED : " << std::endl;
        std::cout << " ---- MESH HAS GENUS 0!" << std::endl;
        exit(1);
    }
    std::cout << "Mesh has genus : " << genus << std::endl;
    //

// compute the bounding box

    return genus;
}
/*
Find the opposite vertex of the edge in the triangle
*/
std::pair<int, int> findOppositeVertex(const _SimpleMesh &mesh, const int edge_idx) {
    auto const& e = mesh.vecEdge[edge_idx];
    int e0 = e.v0;
    int e1 = e.v1;
    
    // Find opposite vertex in first triangle
    int tri_idx_0 = e.AdjTri[0];
    int a0 = mesh.vecTriangle[tri_idx_0].v0;
    int b0 = mesh.vecTriangle[tri_idx_0].v1;
    int c0 = mesh.vecTriangle[tri_idx_0].v2;
    int opp_v_0;
    if (a0 != e0 && a0 != e1) {
        opp_v_0 = a0;
    } else if (b0 != e0 && b0 != e1) {
        opp_v_0 = b0;
    } else {
        opp_v_0 = c0;
    }
    
    // Find opposite vertex in second triangle (if edge has 2 adjacent triangles)
    int opp_v_1;
    if (e.AdjTriNum == 2) {
        int tri_idx_1 = e.AdjTri[1];
        int a1 = mesh.vecTriangle[tri_idx_1].v0;
        int b1 = mesh.vecTriangle[tri_idx_1].v1;
        int c1 = mesh.vecTriangle[tri_idx_1].v2;
        if (a1 != e0 && a1 != e1) {
            opp_v_1 = a1;
        } else if (b1 != e0 && b1 != e1) {
            opp_v_1 = b1;
        } else {
            opp_v_1 = c1;
        }
    } else {
        // If only one triangle, use the same vertex (shouldn't happen for interior edges)
        throw std::runtime_error("Edge has only one adjacent triangle");
    }
    
    return std::make_pair(opp_v_0, opp_v_1);
}
/*
Compute a unit vector such that its dot product with none edge and its dual vectors is close to 0. 
It avoids the case when the direction vector is perpendicular to any edge vectors.
@param mesh: the input mesh
@param uniDirection: the output direction
*/
void RandomUniqueDirection(_SimpleMesh &mesh, Vector3 &uniDirection) {
    std::vector<Vector3> vecEdgeDirections(2 * mesh.vecEdge.size());
    // Compute edge direction vectors for each edge
    // all vectors are pointed to positive Z direction
    for (unsigned int i = 0; i < mesh.vecEdge.size(); i++) {
        vecEdgeDirections[i][0] = mesh.vecVertex[mesh.vecEdge[i].v0].x - mesh.vecVertex[mesh.vecEdge[i].v1].x;
        vecEdgeDirections[i][1] = mesh.vecVertex[mesh.vecEdge[i].v0].y - mesh.vecVertex[mesh.vecEdge[i].v1].y;
        vecEdgeDirections[i][2] = mesh.vecVertex[mesh.vecEdge[i].v0].z - mesh.vecVertex[mesh.vecEdge[i].v1].z;
        // reverse the direction if the vector is pointing to negative Z direction 
        if (vecEdgeDirections[i][2] < 0) {
            for (int j = 0; j < 3; j++)
                vecEdgeDirections[i][j] = -vecEdgeDirections[i][j];
        }
        // normalize the vector
        unitize(vecEdgeDirections[i]);
    }
    int mesh_edge_size = mesh.vecEdge.size();
    int opp_v_0;
    int opp_v_1;
    for (unsigned int i = 0; i < mesh.vecEdge.size(); i++) {
        if (mesh.vecEdge[i].AdjTriNum == 2) {
            auto [opp_v_0, opp_v_1] = findOppositeVertex(mesh, i);
            //
            vecEdgeDirections[i + mesh_edge_size][0] = mesh.vecVertex[opp_v_0].x - mesh.vecVertex[opp_v_1].x;
            vecEdgeDirections[i + mesh_edge_size][1] = mesh.vecVertex[opp_v_0].y - mesh.vecVertex[opp_v_1].y;
            vecEdgeDirections[i + mesh_edge_size][2] = mesh.vecVertex[opp_v_0].z - mesh.vecVertex[opp_v_1].z;
            //
            if (vecEdgeDirections[i + mesh_edge_size][2] < 0) {
                for (int j = 0; j < 3; j++)
                    vecEdgeDirections[i + mesh_edge_size][j] = -vecEdgeDirections[i + mesh_edge_size][j];
            }
            //
            unitize(vecEdgeDirections[i + mesh_edge_size]);
        }
    }
    //
    srand(time(NULL));
    //
    int runTimesCounter = 0;
    int halfRandMax = RAND_MAX >> 2;
    bool bFindDirection = false;
    double minAangleValue = 1.0;
    double error_precision = 1e-7;
    double max_running_counter = 2000;
    while (!bFindDirection) {
        runTimesCounter++;
        for (int i = 0; i < 3; i++) {
            uniDirection[i] = rand() - halfRandMax;
        }
        if (uniDirection[2] < 0) {
            for (int i = 0; i < 3; i++)
                uniDirection[i] = -uniDirection[i];
        }
        //
        unitize(uniDirection);
        //
        bFindDirection = true;
        minAangleValue = 1.0;
        for (unsigned int i = 0; i < vecEdgeDirections.size(); i++) {
            double angleValue = uniDirection * vecEdgeDirections[i];
            if (angleValue < error_precision && angleValue > -error_precision) {
                bFindDirection = false;
                break;
            }
            if (std::abs(minAangleValue) > std::abs(angleValue))
                minAangleValue = angleValue;
        }
        if (runTimesCounter > max_running_counter)
            break;
    }
    //
    //std::cout << "run times : " << runTimesCounter << std::endl;
    //std::cout << "min angle : " << minAangleValue << std::endl;
    //
    return;
}


/*
Validate the direction vector
@param mesh: the input mesh
@param distinctDirection: the input direction
*/
void validateDirectionVector(const _SimpleMesh &mesh, const Vector3 &distinctDirection) {
    std::vector<Vector3> vecEdgeDirections(2 * mesh.vecEdge.size());
    // Compute edge direction vectors for each edge
    // all vectors are pointed to positive Z direction
    for (unsigned int i = 0; i < mesh.vecEdge.size(); i++) {
        vecEdgeDirections[i][0] = mesh.vecVertex[mesh.vecEdge[i].v0].x - mesh.vecVertex[mesh.vecEdge[i].v1].x;
        vecEdgeDirections[i][1] = mesh.vecVertex[mesh.vecEdge[i].v0].y - mesh.vecVertex[mesh.vecEdge[i].v1].y;
        vecEdgeDirections[i][2] = mesh.vecVertex[mesh.vecEdge[i].v0].z - mesh.vecVertex[mesh.vecEdge[i].v1].z;
        // reverse the direction if the vector is pointing to negative Z direction 
        if (vecEdgeDirections[i][2] < 0) {
            for (int j = 0; j < 3; j++)
                vecEdgeDirections[i][j] = -vecEdgeDirections[i][j];
        }
        // normalize the vector
        unitize(vecEdgeDirections[i]);
    }
    int mesh_edge_size = mesh.vecEdge.size();
    int opp_v_0;
    int opp_v_1;
    for (unsigned int i = 0; i < mesh.vecEdge.size(); i++) {
        if (mesh.vecEdge[i].AdjTriNum == 2) {
            auto [opp_v_0, opp_v_1] = findOppositeVertex(mesh, i);
            //
            vecEdgeDirections[i + mesh_edge_size][0] = mesh.vecVertex[opp_v_0].x - mesh.vecVertex[opp_v_1].x;
            vecEdgeDirections[i + mesh_edge_size][1] = mesh.vecVertex[opp_v_0].y - mesh.vecVertex[opp_v_1].y;
            vecEdgeDirections[i + mesh_edge_size][2] = mesh.vecVertex[opp_v_0].z - mesh.vecVertex[opp_v_1].z;
            //
            if (vecEdgeDirections[i + mesh_edge_size][2] < 0) {
                for (int j = 0; j < 3; j++)
                    vecEdgeDirections[i + mesh_edge_size][j] = -vecEdgeDirections[i + mesh_edge_size][j];
            }
            //
            unitize(vecEdgeDirections[i + mesh_edge_size]);
        }
    }
    
    // Ensure the direction vector is normalized
    Vector3 normalizedDir = distinctDirection;
    unitize(normalizedDir);
    
    double error_precision = 1e-7;
    double minDotProduct = 1.0;
    int problematicEdgeIdx = -1;
    for (unsigned int i = 0; i < vecEdgeDirections.size(); i++) {
        double angleValue = normalizedDir * vecEdgeDirections[i];
        double absAngleValue = std::abs(angleValue);
        if (absAngleValue < minDotProduct) {
            minDotProduct = absAngleValue;
            problematicEdgeIdx = i;
        }
        if (angleValue < error_precision && angleValue > -error_precision) {
            std::string errorMsg = "The input direction vector is invalid to produce a Reeb graph. ";
            errorMsg += "Direction is nearly perpendicular to edge direction " + std::to_string(i);
            errorMsg += " (dot product: " + std::to_string(angleValue) + "). ";
            errorMsg += "Minimum dot product found: " + std::to_string(minDotProduct);
            throw std::runtime_error(errorMsg);
        }
    }
}


