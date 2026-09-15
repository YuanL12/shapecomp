#pragma once

#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <taskflow/algorithm/for_each.hpp>
#include <taskflow/core/flow_builder.hpp>
#include <taskflow/taskflow.hpp>

#include "geomp/mesh/Mesh.h"
#include "geomp/mesh/Types.h"

namespace geomp {

/*
    SphereLocator
    Parameterize a closed genus-zero mesh onto the unit sphere (Tutte + stereographic
    + Crane centering) and locate spherical triangles containing query points.
*/
class SphereLocator {
public:
	SphereLocator();
	explicit SphereLocator(Mesh& mesh, std::optional<std::string> Tutte_embedding_type = std::nullopt);

	void construct(Mesh& mesh, std::optional<std::string> Tutte_embedding_type = std::nullopt, int pole_type = 0);

	Mesh* get_mesh_ptr() const;
	int get_pole_index() const;
	Eigen::MatrixXd get_sphere_vertices() const;

	std::pair<int, Eigen::Vector3d> findBaryCoords(const Eigen::Vector3d& point) const;
	std::pair<int, Eigen::Vector3d> findBaryCoords(const std::vector<double>& point) const;
	std::vector<std::pair<int, Eigen::Vector3d>> findBaryCoordsParallel(const Eigen::MatrixXd& points) const;

private:
	int get_Tutte_embedding_type(const std::string& type) const;
	Vector sphereVertex(int index) const;
	std::pair<int, Eigen::Vector3d> findBaryCoordsBrute(const Eigen::Vector3d& p) const;
	bool barycentricOnFace(FaceCIter f, const Eigen::Vector3d& p, Eigen::Vector3d& bary, Eigen::Vector3d& raw) const;

	Mesh* mesh_ptr;
	Eigen::MatrixXd sphere_vertices;
	int pole_index;
};

}  // namespace geomp
