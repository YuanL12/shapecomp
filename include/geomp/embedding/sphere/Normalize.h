/* ==NORMALIZE.H =================================================================================
 *
 *  Small functions that apply a Mobius transform on the vertices of a mesh embedded on the sphere,
 *  so that the center of gravity of those points sits at the center of the sphere
 *
 *	Apply Method recently proposed by Keenan Crane:
 *	
 *	A. Baden, K. Crane, and M. Kazhdan. "Moebius Registration", Eurographics Symp. Geom. Proc.
 *	(2018).
 *  Author:  Patrice Koehl and Yuan Luo
 *  Date:    02/20/2024
 *  Version: 2
 *
   =============================================================================================== */

#pragma once
#include <algorithm>
#include "Mesh.h"

/* =============================================================================================
Class
=============================================================================================== */

class Normalize {
public:
	// perform normalization
	static void normalize(Mesh& mesh,  bool verbose);
	// Move to center of mass to zero
	static void inversion2zero(Mesh& mesh);
	static Vector find_center_of_mass(Mesh& mesh);

private:
	//find position of center of mass
	static void centerOfMass(Mesh& mesh, double *cg);
	// Compute Jacobian
	static void jacobian(Mesh& mesh, double *J);
	// Inverse a 3x3 matrix
	static void inv3x3(double *M, double *Minv);
	// Multiple matrix by vector
	static void matVec(double *Mat, double *X, double *B);
	// Apply Inversion
	static void applyInversion(Mesh& mesh, double *center);
};

/* ==============================================================================================
Find center of massh of current vertices
=============================================================================================== */

void Normalize::centerOfMass(Mesh& mesh, double *cg){
	cg[0] = 0; cg[1] = 0; cg[2] = 0;
	for(VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++){
		cg[0] += v->position2.x;
		cg[1] += v->position2.y;
		cg[2] += v->position2.z;
	}
	int n = mesh.vertices.size();
	cg[0] /= n; cg[1] /= n; cg[2] /= n;
}

// compute center of mass 
Vector Normalize::find_center_of_mass(Mesh& mesh){
	Vector cg;
	cg[0] = 0; cg[1] = 0; cg[2] = 0;
	for(VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++){
		cg[0] += v->position2.x;
		cg[1] += v->position2.y;
		cg[2] += v->position2.z;
	}

	int n = mesh.vertices.size();
	cg[0] /= n; cg[1] /= n; cg[2] /= n;
	return cg;
}

/* ==============================================================================================
Jacobian for energy function
=============================================================================================== */
void Normalize::jacobian(Mesh& mesh, double *J){
	for(int i = 0; i < 9; i++) J[i] = 0;
	Vector point;
	for(VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++){
		point = v->position2;
		J[0] -= point[0]*point[0]; J[1] -= point[0]*point[1]; J[2] -= point[0]*point[2];
		J[4] -= point[1]*point[1]; J[5] -= point[1]*point[2];
		J[8] -= point[2]*point[2];
	}
	int n = mesh.vertices.size();
	J[0] = 1.0 + J[0]/n; J[4] = 1.0 + J[4]/n; J[8] = 1.0 + J[8]/n;
	J[1] = J[1]/n; J[2] = J[2]/n; J[5] = J[5]/n;
	J[3] = J[1]; J[6] = J[2]; J[7] = J[5];
}

/* ==============================================================================================
Inverse a 3x3 matrix
=============================================================================================== */
void Normalize::inv3x3(double *Mat, double *Inv){
	Inv[0] = Mat[4]*Mat[8] - Mat[7]*Mat[5];
	Inv[1] = -(Mat[3]*Mat[8] - Mat[6]*Mat[5]);
	Inv[2] = Mat[3]*Mat[7] - Mat[6]*Mat[4];
	Inv[3] = Inv[1];
	Inv[4] = Mat[0]*Mat[8] - Mat[6]*Mat[2];
	Inv[5] = -(Mat[0]*Mat[7] - Mat[6]*Mat[1]);
	Inv[6] = Inv[2];
	Inv[7] = Inv[5];
	Inv[8] = Mat[0]*Mat[4] - Mat[3]*Mat[1];

	double det = Mat[0]*Inv[0] + Mat[1]*Inv[1] + Mat[2]*Inv[2];
	for(int i = 0; i < 9; i++){
		Inv[i] = Inv[i]/det;
	}
}

/* ==============================================================================================
Matrix vector multiply
=============================================================================================== */
void Normalize::matVec(double *Mat, double *B, double *C){
	C[0] = Mat[0]*B[0] + Mat[1]*B[1] + Mat[2]*B[2];
	C[1] = Mat[3]*B[0] + Mat[4]*B[1] + Mat[5]*B[2];
	C[2] = Mat[6]*B[0] + Mat[7]*B[1] + Mat[8]*B[2];
}

/*
Apply inversion on all vertices
*/
void Normalize::applyInversion(Mesh& mesh, double *center){
	double coef;
	double val, den;
	coef = 1.0 - center[0]*center[0] - center[1]*center[1] - center[2]*center[2];
	Vector point;

	for(VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++){
		point = v->position2;
		den = (point[0]+center[0])*(point[0]+center[0]);
		den += (point[1]+center[1])*(point[1]+center[1]);
		den += (point[2]+center[2])*(point[2]+center[2]);

		val = coef/den;
		
		v->position2.x = val*(point[0]+center[0]) + center[0];
		v->position2.y = val*(point[1]+center[1]) + center[1];
		v->position2.z = val*(point[2]+center[2]) + center[2];
	}
}


/* 
Make the centor of mass zero and preserve orientation by applying invsersion twice:
	1. first make the center of mass zero
	2. reflection along x-y plane to preserve orientation
*/
void Normalize::inversion2zero(Mesh& mesh){
	double cm_length_square;
	double val, den;

	// compute the center of mass
	Vector cm;
	Vector point;
	size_t points_count = 0;
	for (const auto& v : mesh.vertices){
		point = v.position2;
		cm[0] += point[0];
		cm[1] += point[1];
		cm[2] += point[2];
		points_count ++; 
	}
	if (points_count > 0) {
		cm[0] /= static_cast<double>(points_count);
		cm[1] /= static_cast<double>(points_count);
		cm[2] /= static_cast<double>(points_count);
	}

	cm_length_square = cm[0]*cm[0] + cm[1]*cm[1] + cm[2]*cm[2];

	Vector center; // center of inversion circle
	for (size_t i = 0; i < 3; i++)
		center[i] = (1/cm_length_square) * cm[i];

	for(VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++){
		point = v->position2;
		den  = (point[0]-center[0])*(point[0]-center[0]);
		den += (point[1]-center[1])*(point[1]-center[1]);
		den += (point[2]-center[2])*(point[2]-center[2]);

		val = (1/cm_length_square - 1)/den;
		
		v->position2.x = -val*(point[0]-center[0]) - center[0];
		v->position2.y = -val*(point[1]-center[1]) - center[1];
		v->position2.z = val*(point[2]-center[2]) + center[2];
	}
}

  /* ===============================================================================================
   Normalize
   =============================================================================================== */

  void Normalize::normalize(Mesh& mesh, bool verbose = false){

  /* ===============================================================================================
	Apply Method recently proposed by Crane:
	
	A. Baden, K. Crane, and M. Kazhdan. "Moebius Registration", Eurographics Symp. Geom. Proc.
	(2018).
   =============================================================================================== */

	double TOL = 5.e-5;
	int istep = 0;
	double dist, ene;
	double J[9], Jinv[9], COM[3], Center[3], dC[3];

  
	// Initialize center to center of the ball; corresponding Inversion is Identity
	Center[0] = 0.; Center[1] = 0.; Center[2] = 0.;
	// Gauss - Newton iterations for finding appropriate center of inversion
	double coef=0;
	if (verbose){
		std::cout << " " << std::endl;
		std::cout << "Normalize mesh (center vertices):" << std::endl;
		std::cout << "================================" << std::endl;
		std::cout << " " << std::endl;

		std::cout << "        " << "============================================" << std::endl;
		std::cout << "        " << "       Iter       Step size          d(O,G) " << std::endl;
		std::cout << "        " << "============================================" << std::endl;
	}
	
	while(true){
  /* 	==========================================================================================
		Compute center of mass; exit loop if distance to origin < TOL
     	========================================================================================== */
  		centerOfMass(mesh, COM);
		dist = COM[0]*COM[0] + COM[1]*COM[1] + COM[2]*COM[2];
		ene = std::sqrt(dist);

		if(verbose){
			std::cout << "        " << "   " << std::setw(8)<< istep << "    " << std::setw(12) << coef;
			std::cout << "    " << std::setw(12) << ene << std::endl;
		}
		
		if(std::sqrt(dist) < TOL) break;
		istep++;
  
		// Compute Jacobian and its inverse
		jacobian(mesh, J);
		inv3x3(J, Jinv);
  /* 	==========================================================================================
		Compute new Inversion center
		If new center outside the unit ball, bring it back in it
     	========================================================================================== */

		matVec(Jinv, COM, dC);
		Center[0] = -0.5*dC[0]; Center[1] = -0.5*dC[1]; Center[2] = -0.5*dC[2];
		dist = std::sqrt(Center[0]*Center[0] + Center[1]*Center[1] + Center[2]*Center[2]);
//		coef = 1.0;
//		if(dist >= 1) coef = 0.5/dist;
		coef=std::tanh(std::min(dist,0.8))/dist;
		Center[0] = coef*Center[0]; Center[1] = coef*Center[1]; Center[2]=coef*Center[2];

  /* 	==========================================================================================
			Apply Inversion
     	========================================================================================== */

		applyInversion(mesh, Center);

		if(istep > 100) break;
	}
	if (verbose){
		std::cout << "        " << "============================================" << std::endl;
		std::cout << " " << std::endl;
		std::cout << "Sphere vertices centered      : Success" << std::endl;
		std::cout << "Number of steps needed        : " << istep << std::endl;
		std::cout << "Final energy ( i.e. d(O, CG)) : " << ene << std::endl;
		std::cout << " " << std::endl;	
	}
	
	if (istep>100) {
        throw std::runtime_error("Normalization failed because the iteartion steps exceed 100");
    }

  }
