#pragma once

#include "geomp/mesh/HalfEdge.h"

namespace geomp {

class CutPtrIterator {
public:
	// constructor
	CutPtrIterator(HalfEdgeIter he, bool justStarted);

	// increment, comparison and dereference operators
	const CutPtrIterator& operator++();
	bool operator==(const CutPtrIterator& other) const;
	bool operator!=(const CutPtrIterator& other) const;
	WedgeIter operator*() const;

private:
	// members
	HalfEdgeIter currHe;
	bool justStarted;
};

class CutPtrSet {
public:
	// constructors
	CutPtrSet();
	CutPtrSet(HalfEdgeIter he);

	// begin and end routines
	CutPtrIterator begin();
	CutPtrIterator end();

	// begin and end routines
	CutPtrIterator begin() const;
	CutPtrIterator end() const;

private:
	// members
	HalfEdgeIter firstHe;
	bool isInvalid;
};

} // namespace geomp

#include "CutIterator.inl"
