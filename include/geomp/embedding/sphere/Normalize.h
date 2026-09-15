/* ==NORMALIZE.H =================================================================================
 *
 * Möbius inversion that recenters spherical vertices so their center of mass is the origin.
 *
 * A. Baden, K. Crane, and M. Kazhdan. "Moebius Registration", Eurographics Symp. Geom. Proc.
 * (2018).
 *
 * Operates on an Eigen matrix of sphere coordinates; does not mutate Mesh::position.
 *
 * Author:  Patrice Koehl and Yuan Luo
 =============================================================================================== */

#pragma once

#include <Eigen/Dense>

namespace geomp {

class SphereNormalize {
  public:
    static void normalize(Eigen::MatrixXd& sphereVertices, bool verbose = false);
    static Eigen::Vector3d centerOfMass(const Eigen::MatrixXd& sphereVertices);

  private:
    static void jacobian(const Eigen::MatrixXd& sphereVertices, double* J);
    static void inv3x3(const double* M, double* Minv);
    static void matVec(const double* Mat, const double* X, double* B);
    static void applyInversion(Eigen::MatrixXd& sphereVertices, const double* center);
};

}  // namespace geomp
