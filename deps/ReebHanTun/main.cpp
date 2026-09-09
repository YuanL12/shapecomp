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


char *strLicense = "THIS SOFTWARE IS PROVIDED \"AS-IS\". THERE IS NO WARRANTY OF ANY KIND. "
                   "NEITHER THE AUTHORS NOR THE OHIO STATE UNIVERSITY WILL BE LIABLE FOR "
                   "ANY DAMAGES OF ANY KIND, EVEN IF ADVISED OF SUCH POSSIBILITY. \n"
                   "\n"
                   "This software was developed (and is copyrighted by) the Jyamiti group at "
                   "The Ohio State University. Please do not redistribute this software. "
                   "This program is for academic research use only. This software uses the "
                   "CGAL library (www.cgal.org), Boost library (www.boost.org) and Ann library "
                   "(www.cs.umd.edu/~mount/ANN/) which are covered under their own licenses.\n"
                   "\n"
                   "The CGAL library's license "
                   "(which applies to the CGAL library ONLY and NOT to this program itself) is "
                   "as follows:\n"
                   "\n"
                   "LICENSE\n"
                   "---------------------------------------------------------------------------\n"
                   "\n"
                   "The CGAL software consists of several parts, each of which is licensed under "
                   "an open source license. It is also possible to obtain commercial licenses "
                   "from GeometryFactory (www.geometryfactory.com) for all or parts of CGAL. \n"
                   "\n"
                   "The source code of the CGAL library can be found in the directories "
                   "\"src/CGAL\", \"src/CGALQt\" and \"include/CGAL\". It is specified in each file of "
                   "the CGAL library which license applies to it. This is either the GNU Lesser "
                   "General Public License (as published by the Free Software Foundation; "
                   "version 2.1 of the License) or the Q Public License (version 1.0). The texts "
                   "of both licenses can be found in the files LICENSE.LGPL and LICENSE.QPL. \n"
                   "\n"
                   "Distributed along with CGAL (for the users' convenience), but not part of "
                   "CGAL, are the following third-party libraries, available under their own "
                   "licenses: \n"
                   "\n"
                   "- CORE, in the directories \"include/CORE\" and \"src/Core\", is licensed under"
                   " the QPL (see LICENSE.QPL). \n"
                   "- OpenNL, in the directory \"include/OpenNL\", is licensed under the LGPL"
                   " (see include/OpenNL/LICENSE.OPENNL). \n"
                   "- ImageIO, in the directory \"examples/Surface_mesher/ImageIO\", is licensed"
                   " under the LGPL (see LICENSE.LGPL). \n"
                   "\n"
                   "All other files that do not have an explicit copyright notice (e.g., all "
                   "examples and some demos) are licensed under a very permissive license. The "
                   "exact license text can be found in the file LICENSE.FREE_USE. Note that some "
                   "subdirectories have their own copy of LICENSE.FREE_USE. These copies have "
                   "the same license text and differ only in the copyright holder.\n"
                   "---------------------------------------------------------------------------\n"
                   "\n"
                   "The Boost library's license "
                   "(which applies to the Boost library ONLY and NOT to this program itself) is "
                   "as follows:\n"
                   "\n"
                   "LICENSE\n"
                   "---------------------------------------------------------------------------\n"
                   "Boost Software License - Version 1.0 - August 17th, 2003\n"
                   "\n"
                   "Permission is hereby granted, free of charge, to any person or organization "
                   "obtaining a copy of the software and accompanying documentation covered by "
                   "this license (the \"Software\") to use, reproduce, display, distribute, "
                   "execute, and transmit the Software, and to prepare derivative works of the "
                   "Software, and to permit third-parties to whom the Software is furnished to "
                   "do so, all subject to the following: \n"
                   "\n"
                   "The copyright notices in the Software and this entire statement, including "
                   "the above license grant, this restriction and the following disclaimer, "
                   "must be included in all copies of the Software, in whole or in part, and "
                   "all derivative works of the Software, unless such copies or derivative "
                   "works are solely in the form of machine-executable object code generated by "
                   "a source language processor. \n"
                   "\n"
                   "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR "
                   "IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, "
                   "FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT "
                   "SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE "
                   "FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, "
                   "ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER "
                   "DEALINGS IN THE SOFTWARE. \n"
                   "---------------------------------------------------------------------------\n"
                   "\n"
                   "ANN library's license "
                   "(which applies to the ANN library ONLY and NOT to this program itself) is "
                   "as follows: \n"
                   "\n"
                   "LICENSE\n"
                   "---------------------------------------------------------------------------\n"
                   "The ANN Library (all versions) is provided under the terms and "
                   "conditions of the GNU Lesser General Public Library, which is stated "
                   "below.  It can also be found at: \n"
                   "\n"
                   "   http:\//www.gnu.org/copyleft/lesser.html \n"
                   "---------------------------------------------------------------------------\n";

//char * strLicense = "THIS SOFTWARE IS PROVIDED \"AS-IS\". http:\//www.gnu.org/copyleft/lesser.html THERE IS NO WARRANTY OF ANY KIND. \n";
bool ParseCommand(int argc, char **argv, std::string &InputFile, std::string &OutputFile) {
    try {
        /* Define the program options description
        */
        namespace po = boost::program_options;
        po::options_description desc("ReebHanTun Usage");
        desc.add_options()
                (",h", "Help information")
                (",l", "License information")
                (",I", po::value<std::string>(&InputFile)->required(), "Input file name")
                (",O", po::value<std::string>(&OutputFile)->required(), "Output file name prefix");

        // Parser map
        po::variables_map vm;
        try {
            po::store(po::parse_command_line(argc, argv, desc), vm);

            //
            if (vm.count("-h")) {
                std::cout << desc << std::endl;
            }
            //
            if (vm.count("-l")) {
                std::cout << strLicense << std::endl;
            }
            //
            po::notify(vm);
        }
        catch (boost::program_options::required_option &e) {
            std::cerr << "ERROR: " << e.what() << std::endl;
            return false;
        }
        catch (boost::program_options::error &e) {
            std::cerr << "ERROR: " << e.what() << std::endl;
            return false;
        }
        //if (vm.count("-I"))
        //{
        //	std::cout << vm["-I"].as<std::string>() << std::endl;
        //}
        //if (vm.count("-O"))
        //{
        //	std::cout << vm["-O"].as<std::string>() << std::endl;
        //}
    }
    catch (std::exception &e) {
        std::cerr << "Unhandled Exception reached the top of main: "
                  << e.what() << ", application will now exit" << std::endl;
        return false;

    }
    return true;
}

// forward declarations
int LoadData(_SimpleMesh &mesh, std::vector<Vector3> &meshNormal,
             const char *mesh_file_name, double &BoundingBoxRadius, std::vector<int> &tris,
             std::set<int> &extraVertices, int &orgTriangleSize,
             const float fEnlargeFactor);
void RandomUniqueDirection(_SimpleMesh &mesh, Vector3 &uniDirection);

void CycleLocalOptimization(_SimpleMesh &locMesh, psbmReebGraph &reebgraph, std::vector<int> &OrientTriangles,
    std::vector<std::set<int> > &v_basis,
    std::vector<std::set<int> > &h_basis, const float fEnlargeFactor);

void CycleLocalOptimization_bdry(_SimpleMesh &locMesh, psbmReebGraph &reebgraph, std::vector<int> &OrientTriangles,
         std::vector<std::set<int> > &out_v_basis_loops,
         std::vector<std::set<int> > &out_h_basis_loops,
         std::set<int> extraVertices, const float fEnlargeFactor);

int main(int argc, char **argv) {
    std::string InputFileName;
    std::string OutputFileName;
    Vector3 distinctDirection;
    //
    _SimpleMesh mesh;
    const float fEnlargeFactor = 10000.f;
    std::vector<std::set<int> > v_basis_loops;
    std::vector<std::set<int> > h_basis_loops;
    //
    psbmReebGraph reebGraph;
    std::vector<Vector3> meshNormal;
    //
    if (ParseCommand(argc, argv, InputFileName, OutputFileName)) {
        //
        std::vector<int> OrientTriangles;
        std::set<int> extraVertices;
        int nOrgTriangleSize = 0;
        double BoundingBoxRadius;
        LoadData(mesh, meshNormal, InputFileName.c_str(), BoundingBoxRadius, OrientTriangles, extraVertices,
                 nOrgTriangleSize, fEnlargeFactor);
        //////

        RandomUniqueDirection(mesh, distinctDirection);
        //
        //std::cout << distinctDirection[0] << " " <<
        //			 distinctDirection[1] << " " <<
        //			 distinctDirection[2] << std::endl;

        //

        reebGraph.ReserveSpaceForEdges(mesh.vecEdge.size());
        //
        double *scalarField = new double[mesh.vecVertex.size()];
        int perDir = 0;
        int rayDir = 0;
        for (unsigned int i = 0; i < mesh.vecVertex.size(); i++) {
            scalarField[i] = mesh.vecVertex[i].x * distinctDirection[0] +
                             mesh.vecVertex[i].y * distinctDirection[1] +
                             mesh.vecVertex[i].z * distinctDirection[2];
        }
        reebGraph.SetHeightDirection(distinctDirection);
        reebGraph.AssignData(&mesh, scalarField);
        reebGraph.scalarDir = 1;// x=0, y=1, z=2
        //
        //minimumDiffBetweenVertices(scalarField);
        //
        std::cout << std::endl;
        {
            boost::progress_timer t;
            //
            reebGraph.ComputeReebGraph();
            //
        }
        {
        	reebGraph.WriteReebGraphOBJ("casting.obj", mesh.vecVertex);
        	std::cout << "Vertices in RG : " << reebGraph.pVecReebNode->size() << std::endl;
        	std::cout << "Edges in RG : " << reebGraph.pListReebArc->size() << std::endl;
        	std::cout << "Genus is : " << reebGraph.pListReebArc->size() - reebGraph.pVecReebNode->size() + 1  << std::endl;
        }
        {
            boost::progress_timer t;
            ////
            std::cout << "Time for mapping and linking :" << std::endl;
            //reebGraph.ComputeCycleAndPairing();
            reebGraph.ComputingCycle_max_tree();

//            reebGraph.WriteSimplifiedReebGraphOBJ("casting-sim.obj", mesh.vecVertex);

            //std::cout << "mapping" << std::endl;
            //// computing the cycle on surface
            reebGraph.compute_path_on_mesh_for_each_simplified_arc();//pathArcOnMesh, offsetPathArcOnMesh);

            //std::cout << "embed" << std::endl;
            reebGraph.EmbedCycleAsEdgePathOnMesh();
            //std::cout << "linking" << std::endl;
            reebGraph.LinkNumberMatrixComputing();

        }
        //

        //
        if (extraVertices.empty())
            CycleLocalOptimization(mesh, reebGraph, OrientTriangles, v_basis_loops, h_basis_loops,
                                   1.f / fEnlargeFactor);
        else
            CycleLocalOptimization_bdry(mesh, reebGraph, OrientTriangles, v_basis_loops, h_basis_loops, extraVertices,
                                        1.f / fEnlargeFactor);
        //
        std::cout << "Handle and tunnel loops written in files :  \n";
        FilesOutputForOptimalCycles files_out_op;
        files_out_op.InitMeshPtr(&mesh);
        files_out_op.WriteCyclesInformation(OutputFileName.c_str(), v_basis_loops, h_basis_loops);
        int orgVertexSize = mesh.vecVertex.size();
        if (!extraVertices.empty())
            orgVertexSize = *extraVertices.begin();
        files_out_op.WriteGeomviewListFormat(OutputFileName.c_str(), v_basis_loops, h_basis_loops, OrientTriangles,
                                             orgVertexSize, nOrgTriangleSize, 1.f / fEnlargeFactor);

        ///////////////////////////////////////////////////
    }
    return 0;
}