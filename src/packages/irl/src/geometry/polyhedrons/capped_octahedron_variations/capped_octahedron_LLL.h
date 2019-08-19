// This file is part of the Interface Reconstruction Library (IRL)
// a library for interface reconstruction and computational geometry operations
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_CAPPED_OCTAHEDRON_VARIATIONS_CAPPED_OCTAHEDRON_LLL_H_
#define SRC_GEOMETRY_POLYHEDRONS_CAPPED_OCTAHEDRON_VARIATIONS_CAPPED_OCTAHEDRON_LLL_H_

#include "src/geometry/general/stored_vertex_access.h"
#include "src/geometry/half_edge_structures/half_edge_polyhedron.h"
#include "src/geometry/polyhedrons/general_polyhedron.h"
#include "src/geometry/polyhedrons/tet.h"
#include "src/parameters/defined_types.h"
#include "src/geometry/polyhedrons/capped_octahedron_variations/capped_octahedron_correction_base.h"

namespace IRL{

template <class Derived, class VertexType>
class CappedOctahedron_LLLSpecialization: public CappedOctahedronCorrectionBase<Derived, VertexType> {
public: 

HalfEdgePolyhedron<VertexType> generateHalfEdgeVersion(void) const; 

template<class HalfEdgePolyhedronType>
void setHalfEdgeVersion(HalfEdgePolyhedronType* a_half_edge_version) const;

static constexpr UnsignedIndex_t getNumberOfSimplicesInDecomposition(void);

static constexpr std::array<UnsignedIndex_t, 4> getSimplexIndicesFromDecomposition(const UnsignedIndex_t a_tet);

ProxyTet<Derived> getSimplexFromDecomposition(const UnsignedIndex_t a_tet) const;

}; 

template <class VertexType>
class StoredCappedOctahedron_LLL: public StoredVertexAccess<StoredCappedOctahedron_LLL<VertexType>,VertexType,7>, public CappedOctahedron_LLLSpecialization<StoredCappedOctahedron_LLL<VertexType>, VertexType>{
friend StoredVertexAccess<StoredCappedOctahedron_LLL<VertexType>,VertexType,7>;

public: 

using StoredVertexAccess<StoredCappedOctahedron_LLL<VertexType>,VertexType,7>::StoredVertexAccess;

StoredCappedOctahedron_LLL(void) = default;
}; 


// Predefined types 
using CappedOctahedron_LLL= StoredCappedOctahedron_LLL<Pt>; 

} // namespace IRL 

#include "capped_octahedron_LLL.tpp"
#endif //SRC_GEOMETRY_POLYHEDRONS_CAPPED_OCTAHEDRON_VARIATIONS_CAPPED_OCTAHEDRON_LLL_H_
