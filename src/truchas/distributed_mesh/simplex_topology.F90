!!
!! SIMPLEX_TOPOLOGY
!!
!! This module defines a natural convention for labeling the facets of the
!! triangle and tetrahedron cell that is peculiar to simplices.  The module
!! SIMPLEX_GEOMETRY is the companion to this module.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! August 2015.  Adapted and refactored from earlier code.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! USAGE NOTES
!!
!!  An oriented k-simplex (or cell) is an ordered list [n1 n2 ... nk] of
!!  integer (node) indices.  Here we are interested in edges (k=2), triangles
!!  (k=3), and tetrahedrons (k=4).  The ith oriented "face" of a k-simplex is
!!  the (k-1)-simplex that results from dropping the ith index.  Thus the
!!  oriented "faces" of the oriented triangle [n1 n2 n3] are the edges [n2 n3],
!!  [n1 n3], and [n1 n2] -- the ith edge is the edge opposite the ith vertex.
!!  Similarly, the oriented faces of the oriented tetrahedron [n1 n2 n3 n4] are
!!  the triangles [n2 n3 n4], [n1 n3 n4], [n1 n2 n4], and [n1 n2 n3] -- the ith
!!  face is the face opposite the ith vertex.
!!
!!  The parameter arrays TRI_EDGE_VERT and TET_FACE_VERT provide the indices
!!  needed to extract the oriented faces from a triangle and tetrahedron,
!!  respectively.  For example if the array tet(:) = [n1, n2, n3, n4] then
!!  tet(TET_FACE_VERT(:,i)) gives the oriented ith face of the tetrahedron.
!!  NB: THE ORDER OF THE INDICES IN THESE ARRAYS IS DELIBERATE AND ESSENTIAL.
!!
!!  An oriented triangle [n1 n2 n3] in R^3 has an associated unit normal vector
!!  according to the right-hand-rule. Relative to this direction, edges 1 and 3
!!  are oriented in the counter-clockwise direction while edge 2 is oriented in
!!  the clockwise direction.
!!
!!  If the volume of the oriented tetrahedron is positive, then oriented faces
!!  1 and 3 are directed outward while faces 2 and 4 are directed inward.  If
!!  the volume is negative the face orientations are reversed.
!!
!!  Similarly, if the area of an oriented triangle in R^2 is positive, then
!!  edges 1 and 3 are oriented counter-clockwise and edge 2 clockwise.  If the
!!  area is negative the orientations are reversed.
!!
!!  With this convention note that a shared edge of two triangle faces of a
!!  tetrahedron has the same orientation relative to either triangle. Namely,
!!  if edge i of face k is the same as edge i' of face k' then
!!  TET_FACE_VERT(TRI_EDGE_VERT(:,i),k) = TET_FACE_VERT(TRI_EDGE_VERT(:,i'),k')
!!  This establishes a unique orientation of the edges of a tetrahedron.
!!
!!  It is also convenient to label the edges of a tetrahedon.  The parameter
!!  array TET_EDGE_VERT defines the labeling and provides the indices needed
!!  to extract the oriented edges from a tetrahedron:  if the array tet(:) =
!!  [n1, n2, n3, n4] then tet(TET_EDGE_VERT(:,i)) gives the oriented ith edge
!!  of the tetrahedron. The labeling of the edges in a tetrahedron is rather
!!  immaterial, but the orientation of each edge is not.
!!
!!  The parameter array TET_FACE_EDGE gives the edge indices for each face of
!!  the tet: TET_FACE_EDGE(i,k) is the index of the edge of face k opposite
!!  vertex i of the face.  It satisfies the following identity:
!!  TET_EDGE_VERT(:,TET_FACE_EDGE(i,k)) = TET_FACE_VERT(TRI_EDGE_VERT(:,i),k)
!!
!!  The labeling convention just described enables a very powerful technique
!!  for handling the orientation issues that arise when working with degrees
!!  of freedom that live on oriented facets of a simplicial mesh.  Given a
!!  tet mesh with an arbitrary numbering of the nodes, orient the edges [n1 n2]
!!  with n1 < n2, orient the faces [n1 n2 n3] with n1 < n2 < n3, and orient the
!!  tetrahedra [n1 n2 n3 n4] with n1 < n2 < n3 < n4.  It follows then that the
!!  local orientation of the faces and edges of a tetrahedron and the local
!!  orientation of the edges of a triangle, as defined above, coincides with
!!  this orientation of the mesh edges and faces.  Consequently no additional
!!  data is needed to determine the relationship between the global orientation
!!  of an edge/face and its local orientation relative to a cell or face --
!!  the local orientation *is* the global orientation.  Furthermore, once the
!!  orientation of edges, faces, and cells has been established (through the
!!  ordered lists of their vertices) the nodes can be renumbered arbitrarily;
!!  this property -- local orientation == global orientation -- persists!
!!

module simplex_topology

  implicit none
  private

  integer, parameter, public :: TRI_EDGE_VERT(2,3) = &
      reshape(source=[2,3, 1,3, 1,2], shape=[2,3])

  integer, parameter, public :: TET_EDGE_VERT(2,6) = &
      reshape(source=[1,2, 1,3, 1,4, 2,3, 2,4, 3,4], shape=[2,6])

  integer, parameter, public :: TET_FACE_VERT(3,4) = &
      reshape(source=[2,3,4, 1,3,4, 1,2,4, 1,2,3], shape=[3,4])

  integer, parameter, public :: TET_FACE_EDGE(3,4) = &
      reshape(source=[6,5,4, 6,3,2, 5,3,1, 4,2,1], shape=[3,4])

end module simplex_topology
