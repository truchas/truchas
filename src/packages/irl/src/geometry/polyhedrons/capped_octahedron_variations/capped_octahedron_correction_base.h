// This file is part of the Interface Reconstruction Library (IRL)
// a library for interface reconstruction and computational geometry operations
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_CAPPED_OCTAHEDRON_VARIATIONS_CAPPED_OCTAHEDRON_CORRECTION_BASE_H_
#define SRC_GEOMETRY_POLYHEDRONS_CAPPED_OCTAHEDRON_VARIATIONS_CAPPED_OCTAHEDRON_CORRECTION_BASE_H_

#include "src/geometry/polyhedrons/general_polyhedron.h"
#include "src/geometry/polyhedrons/tet.h"
#include "src/parameters/defined_types.h"

namespace IRL{

template <class Derived, class VertexType>
class CappedOctahedronCorrectionBase: public GeneralPolyhedron<Derived, VertexType, ProxyTet<Derived>> {
  Derived& getDerived(void);
  const Derived& getDerived(void) const;
public:

VertexType& operator[](const UnsignedIndex_t a_index);
const VertexType& operator[](const UnsignedIndex_t a_index) const;

// Adjust the 9th point to match the volume supplied.
inline void adjustCapToMatchVolume(const Volume a_correct_volume);

}; 



} // namespace IRL 

#include "src/geometry/polyhedrons/capped_octahedron_variations/capped_octahedron_correction_base.tpp"
#endif //SRC_GEOMETRY_POLYHEDRONS_CAPPED_OCTAHEDRON_VARIATIONS_CAPPED_OCTAHEDRON_CORRECTION_BASE_H_
