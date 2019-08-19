// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_POLYHEDRON_TPP_
#define SRC_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_POLYHEDRON_TPP_

namespace IRL {

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
SegmentedHalfEdgePolyhedron<FaceType, VertexType>
HalfEdgePolyhedron<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                   kMaxVertices, kMaxFaces>::generateSegmentedPolyhedron(void) {
  SegmentedHalfEdgePolyhedron<FaceType, VertexType> a_polyhedron;
  this->setSegmentedPolyhedron(&a_polyhedron);
  return a_polyhedron;
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
template <class SegmentedType>
void HalfEdgePolyhedron<
    PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges, kMaxVertices,
    kMaxFaces>::setSegmentedPolyhedron(SegmentedType* a_polytope) {
  a_polytope->setNumberOfFaces(0);
  a_polytope->setNumberOfVertices(0);
  for (UnsignedIndex_t n = 0; n < this->getNumberOfFaces(); ++n) {
    a_polytope->addFace(&this->getFace(n));
  }
  for (UnsignedIndex_t n = 0; n < this->getNumberOfVertices(); ++n) {
    a_polytope->addVertex(&this->getVertex(n));
  }
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_POLYHEDRON_TPP_
