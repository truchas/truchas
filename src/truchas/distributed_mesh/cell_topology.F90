!!
!! CELL_TOPOLOGY
!!
!! This module defines the convention used for labeling the vertices and faces
!! of some 3D cell types: tetrahedron, pyramid, wedge/prism, and hexahedron.
!! The conventions are consistent with the ExodusII format.
!!
!! Neil Carlson <nnc@lanl.gov>
!! 16 Feb 2006; refactored August 2015
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module cell_topology

  implicit none
  private

  public :: get_face_nodes, get_link_face_nodes, link_edges
  public :: num_cell_faces, cell_face_sizes, cell_face_sig, cell_edges
  public :: normalize_facet, reverse_facet, facet_parity

  !! These data arrays define the intrinsic labeling conventions for the common
  !! 3D cells used in Truchas.  The conventions are consistent with those used
  !! by the ExodusII mesh format:
  !!
  !! * <T>_FACES(<T>_XFACE(k):<T>_XFACE(k+1)-1) are the cell vertices that
  !!   define the kth face oriented out of the cell for cell type <T>.
  !!
  !! * <T>_FSIZE(k) is the number of vertices on the kth face of cell type <T>.
  !!
  !! * <T>_EDGES(:,k) are the 2 cell vertices defining edge k of cell type <T>.
  !!
  !! * <T>_FACE_SIG(k) is a bit mask array tagging the vertices on the kth face
  !!   of cell type <T>.  BTEST(<T>_FACE_SIG(k), i-1) is true if vertex i is
  !!   contained in face k of cell type <T>.
  !!
  !! * <T>_VERT_SIG(k) is a bit mask array tagging the faces that contain the
  !!   kth vertex for cell type <T>.  BTEST(<T>_VERT_SIG(k), i-1) is true if
  !!   face i contains in vertex k of cell type <T>.

  integer, target, public :: TET4_XFACE(5), TET4_FACES(12), TET4_FSIZE(4), TET4_EDGES(2,6)
  integer, target, public :: TET4_FACE_SIG(4), TET4_VERT_SIG(4)
  data TET4_XFACE/1,4,7,10,13/
  data TET4_FACES/1,2,4, 2,3,4, 1,4,3, 1,3,2/
  data TET4_FSIZE/3,3,3,3/
  data TET4_EDGES/1,2, 1,3, 1,4, 2,3, 2,4, 3,4/
  data TET4_FACE_SIG/b'1011', b'1110', b'1101', b'0111'/
  data TET4_VERT_SIG/b'1101', b'1011', b'1110', b'0111'/

  integer, target, public :: PYR5_XFACE(6), PYR5_FACES(16), PYR5_FSIZE(5), PYR5_EDGES(2,8)
  integer, target, public :: PYR5_FACE_SIG(5), PYR5_VERT_SIG(5)
  data PYR5_XFACE/1,4,7,10,13,17/
  data PYR5_FACES/1,2,5,  2,3,5,  3,4,5,  1,5,4,  1,4,3,2/
  data PYR5_FSIZE/3,3,3,3,4/
  data PYR5_EDGES/1,2, 1,4, 1,5, 2,3, 2,5, 3,4, 3,5, 4,5/
  data PYR5_FACE_SIG/b'10011', b'10110', b'11100', b'11001', b'01111'/
  data PYR5_VERT_SIG/b'11001', b'10011', b'10110', b'11100', b'01111'/

  integer, target, public :: WED6_XFACE(6), WED6_FACES(18), WED6_FSIZE(5), WED6_EDGES(2,9)
  integer, target, public :: WED6_FACE_SIG(5), WED6_VERT_SIG(6)
  data WED6_XFACE/1,5,9,13,16,19/
  data WED6_FACES/1,2,5,4,  2,3,6,5,  1,4,6,3,  1,3,2,  4,5,6/
  data WED6_FSIZE/4,4,4,3,3/
  data WED6_EDGES/1,2, 1,3, 1,4, 2,3, 2,5, 3,6, 4,5, 4,6, 5,6/
  data WED6_FACE_SIG/b'011011', b'110110', b'101101', b'000111', b'111000'/
  data WED6_VERT_SIG/b'01101', b'01011', b'01110', b'10101', b'10011', b'10110'/

  integer, target, public :: HEX8_XFACE(7), HEX8_FACES(24), HEX8_FSIZE(6), HEX8_EDGES(2,12)
  integer, target, public :: HEX8_FACE_SIG(6), HEX8_VERT_SIG(8)
  data HEX8_XFACE/1,5,9,13,17,21,25/
  data HEX8_FACES/1,2,6,5, 2,3,7,6, 3,4,8,7, 1,5,8,4, 1,4,3,2, 5,6,7,8/
  data HEX8_FSIZE/4,4,4,4,4,4/
  data HEX8_EDGES/1,2, 1,4, 1,5, 2,3, 2,6, 3,4, 3,7, 4,8, 5,6, 5,8, 6,7, 7,8/
  data HEX8_FACE_SIG/b'00110011', b'01100110', b'11001100', b'10011001', b'00001111', b'11110000'/
  data HEX8_VERT_SIG/b'011001', b'010011', b'010110', b'011100', b'101001', b'100011', b'100110', b'101100'/

  !! These provide the list of local indices of faces that contain each of the
  !! cell vertices: *_VERT_FACE(:,k) are the three face indices that contain
  !! vertex k of the cell.  NB: Clients may assume the lists are sorted.
  !! Note that pyramid elements are not supported.

  integer, parameter, public :: HEX8_VERT_FACE(3,8) = &
      reshape(source=[1,4,5, 1,2,5, 2,3,5, 3,4,5, 1,4,6, 1,2,6, 2,3,6, 3,4,6], shape=[3,8])

  integer, parameter, public :: WED6_VERT_FACE(3,6) = &
      reshape(source=[1,3,4, 1,2,4, 2,3,4, 1,3,5, 1,2,5, 2,3,5], shape=[3,6])

  integer, parameter, public :: TET4_VERT_FACE(3,4) = &
      reshape(source=[2,3,4, 1,3,4, 1,2,4, 1,2,3], shape=[3,4])

  !! An interface link appears formally equivalent to a wedge or hex cell, but
  !! its sole purpose is to provide a link between the bottom and top face of
  !! the cell (4 and 5 for the wedge, and 5 and 6 for the hex). These parameters
  !! provide the lists of vertices on those pair of faces with an inward directed
  !! orientation, which would be consistent with an outward orientation of the
  !! matching faces of the neighboring cells that the interface link connects.

  integer, parameter :: TRI_LINK_FACE_VERT(3,2) = reshape(source=[1,2,3, 4,6,5], shape=[3,2])
  integer, parameter :: QUAD_LINK_FACE_VERT(4,2) = reshape(source=[1,2,3,4, 5,8,7,6], shape=[4,2])

  !! These data arrays provide the pair of vertices defining the link edges.
  integer, target :: TRI_LINK_EDGES(2,3), QUAD_LINK_EDGES(2,4)
  data TRI_LINK_EDGES/1,4, 2,5, 3,6/
  data QUAD_LINK_EDGES/1,5, 2,6, 3,7, 4,8/

contains

  !! Returns a pointer to the array of cell edges for the given cell type.
  !! CNODES is the list of nodes defining the cell, but only the size of the
  !! array is significant and is used as a surrogate for the cell type.  This
  !! works for tet, pyramid, wedge, and hex cells; a null pointer is returned
  !! for anything else.

  function cell_edges (cnodes) result (edges)
    integer, intent(in) :: cnodes(:)
    integer, pointer :: edges(:,:)
    select case (size(cnodes))
    case (4)  ! tet
      edges => TET4_EDGES
    case (5)  ! pyramid
      edges => PYR5_EDGES
    case (6)  ! wedge
      edges => WED6_EDGES
    case (8)  ! hex
      edges => HEX8_EDGES
    case default
      edges => null()
    end select
  end function cell_edges

  !! Returns the number of faces the given cell possesses.  CNODES is the list
  !! of nodes defining the cell, but only the size of the array is significant;
  !! the values are ignored.  We are using the size of the cell connectivity as
  !! a surrogate for the element type.  This works for tet, pyramid, wedge, and
  !! hex cells; for anything else 0 is returned.

  pure integer function num_cell_faces (cnodes)
    integer, intent(in) :: cnodes(:)
    select case (size(cnodes))
    case (4)  ! tet
      num_cell_faces = 4
    case (5)  ! pyramid
      num_cell_faces = 5
    case (6)  ! wedge
      num_cell_faces = 5
    case (8)  ! hex
      num_cell_faces = 6
    case default
      num_cell_faces = 0
    end select
  end function num_cell_faces

  !! Returns a pointer to the cell face signature array for the given cell type.
  !! CNODES is the list of nodes defining the cell, but only the size of the
  !! array is significant and is used as a surrogate for the cell type.  This
  !! works for tet, pyramid, wedge, and hex cells; a null pointer is returned
  !! for anything else.

  function cell_face_sig (cnodes) result (face_sig)
    integer, intent(in) :: cnodes(:)
    integer, pointer :: face_sig(:)
    select case (size(cnodes))
    case (4)  ! tet
      face_sig => TET4_FACE_SIG
    case (5)  ! pyramid
      face_sig => PYR5_FACE_SIG
    case (6)  ! wedge
      face_sig => WED6_FACE_SIG
    case (8)  ! hex
      face_sig => HEX8_FACE_SIG
    case default
      face_sig => null()
    end select
  end function cell_face_sig

  !! Returns a pointer to the array of cell face sizes for the given cell type.
  !! CNODES is the list of nodes defining the cell, but only the size of the
  !! array is significant and is used as a surrogate for the cell type.  This
  !! works for tet, pyramid, wedge, and hex cells; a null pointer is returned
  !! for anything else.

  function cell_face_sizes (cnodes) result (fsizes)
    integer, intent(in) :: cnodes(:)
    integer, pointer :: fsizes(:)
    select case (size(cnodes))
    case (4)  ! tet
      fsizes => TET4_FSIZE
    case (5)  ! pyramid
      fsizes => PYR5_FSIZE
    case (6)  ! wedge
      fsizes => WED6_FSIZE
    case (8)  ! hex
      fsizes => HEX8_FSIZE
    case default
      fsizes => null()
    end select
  end function cell_face_sizes

  !! Returns the list of nodes defining the specified face of the given cell.
  !! CNODES is the list of nodes defining the cell and INDEX the face of the
  !! cell.  The result is returned in the allocatable array FNODES and is
  !! oriented outward with respect to the cell.  If NORMALIZE is present with
  !! value true, the list FNODES begins with the smallest node number.  If
  !! REVERSE is present with value true, the returned FNODES list is oriented
  !! inward with respect to the cell.  This works for tet, pyramid, wedge, and
  !! hex cells; FNODES is returned deallocated for anything else.

  pure subroutine get_face_nodes (cnodes, index, fnodes, normalize, reverse)

    integer, intent(in) :: cnodes(:), index
    integer, allocatable, intent(inout) :: fnodes(:)
    logical, intent(in), optional :: normalize, reverse

    select case (size(cnodes))
    case (4)  ! tet
      fnodes = cnodes(TET4_FACES(TET4_XFACE(index):TET4_XFACE(index+1)-1))
    case (5)  ! pyramid
      fnodes = cnodes(PYR5_FACES(PYR5_XFACE(index):PYR5_XFACE(index+1)-1))
    case (6)  ! wedge
      fnodes = cnodes(WED6_FACES(WED6_XFACE(index):WED6_XFACE(index+1)-1))
    case (8)  ! hex
      fnodes = cnodes(HEX8_FACES(HEX8_XFACE(index):HEX8_XFACE(index+1)-1))
    case default
      if (allocated(fnodes)) deallocate(fnodes)
      return
    end select

    if (present(normalize)) then
      if (normalize) fnodes = cshift(fnodes, shift=minloc(fnodes,dim=1)-1)
    end if

    if (present(reverse)) then
      if (reverse) call reverse_facet (fnodes)
    end if

  end subroutine get_face_nodes

  !! This function return a pointer to the list of nodes defining the specified
  !! oriented link face K in {1,2}. The faces are oriented outward with respect
  !! to the cells they belongs to; if we regard the link nodes as defining a
  !! wedge or hex element, the faces are inward oriented with respect to it.
  !! The pointer target is allocated by the function and the caller is
  !! responsible for deallocating it when it is no longer needed.

  pure subroutine get_link_face_nodes (lnodes, k, fnodes, normalize, reverse)

    integer, intent(in) :: lnodes(:), k
    integer, allocatable, intent(inout) :: fnodes(:)
    logical, intent(in), optional :: normalize, reverse

    select case (size(lnodes))
    case (6) ! triangular link faces
      fnodes = lnodes(TRI_LINK_FACE_VERT(:,k))
    case (8) ! quadrilateral link faces
      fnodes = lnodes(QUAD_LINK_FACE_VERT(:,k))
    case default
      if (allocated(fnodes)) deallocate(fnodes)
      return
    end select

    !! If specified, rotate the smallest value to the initial position.
    if (present(normalize)) then
      if (normalize) call normalize_facet (fnodes)
    end if

    !! If specified, reverse the orientation of the face.
    if (present(reverse)) then
      if (reverse) call reverse_facet (fnodes)
    end if

  end subroutine get_link_face_nodes

  !! Returns a pointer to the array of link edges for the given link type.
  !! LNODES is the list of nodes defining the link, but only the size of the
  !! array is significant and is used as a surrogate for the link type.

  function link_edges (lnodes)
    integer, intent(in) :: lnodes(:)
    integer, pointer :: link_edges(:,:)
    select case (size(lnodes))
    case (6)
      link_edges => TRI_LINK_EDGES
    case (8)
      link_edges => QUAD_LINK_EDGES
    case default
      link_edges => null()
    end select
  end function link_edges

  !! An facet is an edge in R^2, or triangle, quadrilateral, or general polygon
  !! in R^3.  It is described by a list of the distinct node numbers that define
  !! it, with the order of the nodes defining an orientation for the face.
  !! In R^3 the nodes are ordered counter-clockwise about the face with respect
  !! to the choice of normal vector to the surface.  In R^2 there is a unique
  !! description of an oriented edge, but in R^3 there are multiple equivalent
  !! descriptions of an oriented face, depending on the initial node.  A
  !! description of a face is normalized when the initial node is the one having
  !! the smallest number.

  !! This subroutine normalizes the LIST of node numbers defining a face,
  !! preserving its orientation.  It is a no-op for an edge in R^2.

  pure subroutine normalize_facet (list)
    integer, intent(inout) :: list(:)
    if (size(list) > 2) list = cshift(list, shift=minloc(list,dim=1)-1)
  end subroutine normalize_facet

  !! Thus subroutine reverses the orientation of the given facet.  In R^3,
  !! the returned facet remains normalized if the input facet is normalized.

  pure subroutine reverse_facet (list)
    integer, intent(inout) :: list(:)
    integer :: j1, j2
    select case (size(list))
    case (2)
      call swap (list(1), list(2))
    case (3)
      call swap (list(2), list(3))
    case (4)
      call swap (list(2), list(4))
    case (5:)
      j1 = 2
      j2 = size(list)
      do while (j1 < j2)
        call swap (list(j1), list(j2))
        j1 = j1 + 1
        j2 = j2 - 1
      end do
    end select
  contains
    pure subroutine swap (a, b)
      integer, intent(inout) :: a, b
      integer :: n
      n = a; a = b; b = n
    end subroutine
  end subroutine reverse_facet

  !! This function returns the relative parity of the oriented facets F1 and F2:
  !! it returns 1 if the facets are the same and have the same orientation, -1 if
  !! the facets are the same but have opposite orientations, and 0 otherwise.
  !! Faces in R^3 must be normalized.

  pure integer function facet_parity (f1, f2) result (parity)

    integer, intent(in) :: f1(:), f2(:)

    integer :: n, j

    parity = 0

    if (size(f1) /= size(f2)) return

    select case (size(f1))

    case (2)  ! edge facets

      if (f1(1) == f2(1)) then
        if (f1(2) == f2(2)) parity = 1
      else if (f1(1) == f2(2)) then
        if (f1(2) == f2(1)) parity = -1
      end if

    case (3)  ! tri facets

      if (f1(1) == f2(1)) then
        if (f1(2) == f2(2)) then
          if (f1(3) == f2(3)) parity = 1
        else if (f1(2) == f2(3)) then
          if (f1(3) == f2(2)) parity = -1
        end if
      end if

    case (4)  ! quad facets

      if (f1(1) == f2(1)) then
        if (f1(2) == f2(2)) then
          if (f1(3) /= f2(3)) return
          if (f1(4) /= f2(4)) return
          parity = 1
        else if (f1(2) == f2(4)) then
          if (f1(3) /= f2(3)) return
          if (f1(4) /= f2(2)) return
          parity = -1
        end if
      end if

    case (5:) ! more general facets

      n = size(f1)
      if (f1(1) == f2(1)) then
        if (f1(2) == f2(2)) then
          do j = 3, n
            if (f1(j) /= f2(j)) return
          end do
          parity = 1
        else if (f1(2) == f2(n)) then
          do j = 3, n
            if (f1(j) /= f2(n+2-j)) return
          end do
          parity = -1
        end if
      end if

    end select

  end function facet_parity

end module cell_topology
