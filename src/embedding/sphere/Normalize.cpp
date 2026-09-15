#include "geomp/embedding/sphere/Normalize.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace geomp {

Eigen::Vector3d SphereNormalize::centerOfMass(const Eigen::MatrixXd& sphereVertices) {
	if (sphereVertices.rows() == 0 || sphereVertices.cols() != 3) {
		throw std::runtime_error("SphereNormalize: expected an (N, 3) matrix");
	}
	Eigen::Vector3d cg = Eigen::Vector3d::Zero();
	for (int i = 0; i < sphereVertices.rows(); ++i) {
		cg += sphereVertices.row(i).transpose();
	}
	cg /= static_cast<double>(sphereVertices.rows());
	return cg;
}

void SphereNormalize::jacobian(const Eigen::MatrixXd& sphereVertices, double* J) {
	for (int i = 0; i < 9; ++i) {
		J[i] = 0.0;
	}
	for (int i = 0; i < sphereVertices.rows(); ++i) {
		const double x = sphereVertices(i, 0);
		const double y = sphereVertices(i, 1);
		const double z = sphereVertices(i, 2);
		J[0] -= x * x;
		J[1] -= x * y;
		J[2] -= x * z;
		J[4] -= y * y;
		J[5] -= y * z;
		J[8] -= z * z;
	}
	const int n = static_cast<int>(sphereVertices.rows());
	J[0] = 1.0 + J[0] / n;
	J[4] = 1.0 + J[4] / n;
	J[8] = 1.0 + J[8] / n;
	J[1] = J[1] / n;
	J[2] = J[2] / n;
	J[5] = J[5] / n;
	J[3] = J[1];
	J[6] = J[2];
	J[7] = J[5];
}

void SphereNormalize::inv3x3(const double* Mat, double* Inv) {
	Inv[0] = Mat[4] * Mat[8] - Mat[7] * Mat[5];
	Inv[1] = -(Mat[3] * Mat[8] - Mat[6] * Mat[5]);
	Inv[2] = Mat[3] * Mat[7] - Mat[6] * Mat[4];
	Inv[3] = Inv[1];
	Inv[4] = Mat[0] * Mat[8] - Mat[6] * Mat[2];
	Inv[5] = -(Mat[0] * Mat[7] - Mat[6] * Mat[1]);
	Inv[6] = Inv[2];
	Inv[7] = Inv[5];
	Inv[8] = Mat[0] * Mat[4] - Mat[3] * Mat[1];

	const double det = Mat[0] * Inv[0] + Mat[1] * Inv[1] + Mat[2] * Inv[2];
	if (std::abs(det) < 1e-18) {
		throw std::runtime_error("SphereNormalize: Jacobian is singular");
	}
	for (int i = 0; i < 9; ++i) {
		Inv[i] /= det;
	}
}

void SphereNormalize::matVec(const double* Mat, const double* B, double* C) {
	C[0] = Mat[0] * B[0] + Mat[1] * B[1] + Mat[2] * B[2];
	C[1] = Mat[3] * B[0] + Mat[4] * B[1] + Mat[5] * B[2];
	C[2] = Mat[6] * B[0] + Mat[7] * B[1] + Mat[8] * B[2];
}

void SphereNormalize::applyInversion(Eigen::MatrixXd& sphereVertices, const double* center) {
	const double coef = 1.0 - center[0] * center[0] - center[1] * center[1] - center[2] * center[2];
	for (int i = 0; i < sphereVertices.rows(); ++i) {
		const double x = sphereVertices(i, 0);
		const double y = sphereVertices(i, 1);
		const double z = sphereVertices(i, 2);
		double den = (x + center[0]) * (x + center[0]);
		den += (y + center[1]) * (y + center[1]);
		den += (z + center[2]) * (z + center[2]);
		if (den <= 1e-18) {
			continue;
		}
		const double val = coef / den;
		sphereVertices(i, 0) = val * (x + center[0]) + center[0];
		sphereVertices(i, 1) = val * (y + center[1]) + center[1];
		sphereVertices(i, 2) = val * (z + center[2]) + center[2];
	}
}

void SphereNormalize::normalize(Eigen::MatrixXd& sphereVertices, bool verbose) {
	if (sphereVertices.rows() == 0 || sphereVertices.cols() != 3) {
		throw std::runtime_error("SphereNormalize: expected an (N, 3) matrix");
	}

	const double TOL = 5.e-5;
	int istep = 0;
	double dist = 0.0;
	double ene = 0.0;
	double J[9], Jinv[9], COM[3], Center[3], dC[3];
	Center[0] = 0.0;
	Center[1] = 0.0;
	Center[2] = 0.0;
	double coef = 0.0;

	if (verbose) {
		std::cout << "\nNormalize sphere vertices:\n";
		std::cout << "        ============================================\n";
		std::cout << "               Iter       Step size          d(O,G)\n";
		std::cout << "        ============================================\n";
	}

	while (true) {
		const Eigen::Vector3d cg = centerOfMass(sphereVertices);
		COM[0] = cg[0];
		COM[1] = cg[1];
		COM[2] = cg[2];
		dist = COM[0] * COM[0] + COM[1] * COM[1] + COM[2] * COM[2];
		ene = std::sqrt(dist);

		if (verbose) {
			std::cout << "           " << std::setw(8) << istep << "    " << std::setw(12) << coef;
			std::cout << "    " << std::setw(12) << ene << "\n";
		}

		if (ene < TOL) {
			break;
		}
		++istep;

		jacobian(sphereVertices, J);
		inv3x3(J, Jinv);
		matVec(Jinv, COM, dC);
		Center[0] = -0.5 * dC[0];
		Center[1] = -0.5 * dC[1];
		Center[2] = -0.5 * dC[2];
		dist = std::sqrt(Center[0] * Center[0] + Center[1] * Center[1] + Center[2] * Center[2]);
		if (dist <= 1e-18) {
			break;
		}
		coef = std::tanh(std::min(dist, 0.8)) / dist;
		Center[0] *= coef;
		Center[1] *= coef;
		Center[2] *= coef;
		applyInversion(sphereVertices, Center);

		if (istep > 100) {
			break;
		}
	}

	if (verbose) {
		std::cout << "        ============================================\n";
		std::cout << "Sphere vertices centered      : Success\n";
		std::cout << "Number of steps needed        : " << istep << "\n";
		std::cout << "Final energy ( i.e. d(O, CG)) : " << ene << "\n";
	}

	if (istep > 100) {
		throw std::runtime_error("SphereNormalize failed because the iteration steps exceed 100");
	}

	for (int i = 0; i < sphereVertices.rows(); ++i) {
		const double n = sphereVertices.row(i).norm();
		if (n > 1e-18) {
			sphereVertices.row(i) /= n;
		}
	}
}

}  // namespace geomp
