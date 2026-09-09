#pragma once

#include <cmath>
#include <Eigen/Dense>
namespace geomp {

class Vector {
public:
	// initializes all components to zero
	Vector();

	// initializes with specified components
	Vector(double x, double y, double z = 0.0);

	// copy constructor
	Vector(const Vector& v);

	// constructor by a pair of double
	Vector(const std::pair<double, double>& p);

	// access
	double& operator[](int index);
	const double& operator[](int index) const;

	// print << 
	friend std::ostream& operator<<(std::ostream& os, const Vector& v);

	// math
	Vector operator*(double s) const;
	Vector operator/(double s) const;
	Vector operator+(const Vector& v) const;
	Vector operator-(const Vector& v) const;
	Vector operator-() const;

	Vector& operator*=(double s);
	Vector& operator/=(double s);
	Vector& operator+=(const Vector& v);
	Vector& operator-=(const Vector& v);

	// returns Euclidean length
	double norm() const;

	// returns Euclidean length squared
	double norm2() const;

	// normalizes vector
	void normalize();

	// returns unit vector in the direction of this vector
	Vector unit() const;

	// convert custom Vector to Eigen::Vector2d
	inline Eigen::Vector2d to2DEigenVector() const {
		return Eigen::Vector2d(x, y);
	}
	// members
	double x, y, z;
};

Vector operator*(double s, const Vector& v);
double dot(const Vector& u, const Vector& v);
Vector cross(const Vector& u, const Vector& v);

} // namespace geomp

#include "Vector.inl"
