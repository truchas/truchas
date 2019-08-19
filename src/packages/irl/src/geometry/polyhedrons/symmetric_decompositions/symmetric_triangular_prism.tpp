// This file is part of the Interface Reconstruction Library (IRL)
// a library for interface reconstruction and computational geometry operations
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_SYMMETRIC_DECOMPOSITIONS_SYMMETRIC_TRIANGULAR_PRISM_TPP_
#define SRC_GEOMETRY_POLYHEDRONS_SYMMETRIC_DECOMPOSITIONS_SYMMETRIC_TRIANGULAR_PRISM_TPP_

#include <cassert>

#include "src/geometry/general/moment_calculation_through_simplices.h"
#include "src/geometry/half_edge_structures/half_edge.h"

namespace IRL {

namespace symmetric_triangular_prism_triangulation {
static constexpr UnsignedIndex_t datum_index = 0;
static constexpr std::array<std::array<UnsignedIndex_t, 3>, 12>
    face_triangle_decomposition{{{6, 1, 2},
                                 {7, 4, 1},
                                 {7, 3, 4},
                                 {8, 2, 1},
                                 {8, 1, 4},
                                 {8, 4, 5},
                                 {8, 5, 2},
                                 {9, 5, 3},
                                 {9, 2, 5},
                                 {10, 4, 3},
                                 {10, 3, 5},
                                 {10, 5, 4}}};
}  // namespace symmetric_triangular_prism_triangulation

template <class Derived, class VertexType>
HalfEdgePolyhedron<VertexType> SymmetricTriangularPrismSpecialization<
    Derived, VertexType>::generateHalfEdgeVersion(void) const {
  HalfEdgePolyhedron<VertexType> half_edge_version;
  this->setHalfEdgeVersion(&half_edge_version);
  return half_edge_version;
}

template <class Derived, class VertexType>
template <class HalfEdgePolyhedronType>
void SymmetricTriangularPrismSpecialization<Derived, VertexType>::
    setHalfEdgeVersion(HalfEdgePolyhedronType* a_half_edge_version) const {
  using HalfEdgeType = typename HalfEdgePolyhedronType::half_edge_type;

  a_half_edge_version->resize(54, 11, 18);

  for (UnsignedIndex_t v = 0; v < 11; ++v) {
    a_half_edge_version->getVertex(v).setLocation((*this)[v]);
  }

  static constexpr std::array<UnsignedIndex_t, 54> ending_vertex_mapping{
      {0, 1, 6, 1, 2, 6, 2, 0, 6, 1, 0, 7,  4, 1, 7,  3, 4, 7,
       0, 3, 7, 2, 1, 8, 1, 4, 8, 4, 5, 8,  5, 2, 8,  0, 2, 9,
       3, 0, 9, 5, 3, 9, 2, 5, 9, 4, 3, 10, 3, 5, 10, 5, 4, 10}};
  static constexpr std::array<UnsignedIndex_t, 54> previous_half_edge_mapping{
      {2,  0,  1,  5,  3,  4,  8,  6,  7,  11, 9,  10, 14, 12, 13, 17, 15, 16,
       20, 18, 19, 23, 21, 22, 26, 24, 25, 29, 27, 28, 32, 30, 31, 35, 33, 34,
       38, 36, 37, 41, 39, 40, 44, 42, 43, 47, 45, 46, 50, 48, 49, 53, 51, 52}};
  static constexpr std::array<UnsignedIndex_t, 54> next_half_edge_mapping{
      {1,  2,  0,  4,  5,  3,  7,  8,  6,  10, 11, 9,  13, 14, 12, 16, 17, 15,
       19, 20, 18, 22, 23, 21, 25, 26, 24, 28, 29, 27, 31, 32, 30, 34, 35, 33,
       37, 38, 36, 40, 41, 39, 43, 44, 42, 46, 47, 45, 49, 50, 48, 52, 53, 51}};
  static constexpr std::array<UnsignedIndex_t, 54> face_mapping{
      {0,  0,  0,  1,  1,  1,  2,  2,  2,  3,  3,  3,  4,  4,  4,  5,  5,  5,
       6,  6,  6,  7,  7,  7,  8,  8,  8,  9,  9,  9,  10, 10, 10, 11, 11, 11,
       12, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15, 16, 16, 16, 17, 17, 17}};
  static constexpr std::array<UnsignedIndex_t, 54> opposite_half_edge_mapping{
      {8,  10, 3,  2,  22, 6,  5,  34, 0,  14, 1,  18, 17, 25, 9,  20, 46, 12,
       11, 37, 15, 32, 4,  24, 23, 13, 27, 26, 52, 30, 29, 43, 21, 38, 7,  42,
       41, 19, 33, 44, 49, 36, 35, 31, 39, 53, 16, 48, 47, 40, 51, 50, 28, 45}};
  for (UnsignedIndex_t n = 0;
       n < static_cast<UnsignedIndex_t>(ending_vertex_mapping.size()); ++n) {
    HalfEdgeType& current_half_edge = a_half_edge_version->getHalfEdge(n);
    current_half_edge = HalfEdgeType(
        &a_half_edge_version->getVertex(ending_vertex_mapping[n]),
        &a_half_edge_version->getHalfEdge(previous_half_edge_mapping[n]),
        &a_half_edge_version->getHalfEdge(next_half_edge_mapping[n]),
        &a_half_edge_version->getFace(face_mapping[n]));
    current_half_edge.setOppositeHalfEdge(
        &a_half_edge_version->getHalfEdge(opposite_half_edge_mapping[n]));
    current_half_edge.getFace()->setStartingHalfEdge(&current_half_edge);
    current_half_edge.getVertex()->setHalfEdge(&current_half_edge);
  }
}

template <class Derived, class VertexType>
constexpr UnsignedIndex_t SymmetricTriangularPrismSpecialization<
    Derived, VertexType>::getNumberOfSimplicesInDecomposition(void) {
  return static_cast<UnsignedIndex_t>(
      symmetric_triangular_prism_triangulation::face_triangle_decomposition
          .size());
}

template <class Derived, class VertexType>
constexpr std::array<UnsignedIndex_t, 4>
SymmetricTriangularPrismSpecialization<Derived, VertexType>::
    getSimplexIndicesFromDecomposition(const UnsignedIndex_t a_tet) {
  assert(a_tet < SymmetricTriangularPrismSpecialization::
                     getNumberOfSimplicesInDecomposition());
  return {symmetric_triangular_prism_triangulation::face_triangle_decomposition
              [a_tet][0],
          symmetric_triangular_prism_triangulation::face_triangle_decomposition
              [a_tet][1],
          symmetric_triangular_prism_triangulation::face_triangle_decomposition
              [a_tet][2],
          symmetric_triangular_prism_triangulation::datum_index};
}

template <class Derived, class VertexType>
ProxyTet<Derived> SymmetricTriangularPrismSpecialization<Derived, VertexType>::
    getSimplexFromDecomposition(const UnsignedIndex_t a_tet) const {
  assert(a_tet < this->getNumberOfSimplicesInDecomposition());
  return {static_cast<const Derived&>(*this),
          this->getSimplexIndicesFromDecomposition(a_tet)};
}

}  // namespace IRL
#endif  // SRC_GEOMETRY_POLYHEDRONS_SYMMETRIC_DECOMPOSITIONS_SYMMETRIC_TRIANGULAR_PRISM_TPP_
