#pragma once

#include "geomp/mesh/Mesh.h"

namespace geomp {

template<typename Element, std::vector<Element> Mesh::*meshList, typename ElementIter, typename T>
class MeshData {
public:
	// constructor
	MeshData(const Mesh& mesh);

	// constructor
	MeshData(const Mesh& mesh, const T& initVal);

	// operator[]
	T& operator[](ElementIter e);
	const T& operator[](ElementIter e) const;

private:
	// member
	std::vector<T> data;
};

} // namespace geomp

#include "MeshData.inl"
