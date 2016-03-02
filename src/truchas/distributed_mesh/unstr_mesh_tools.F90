!!
!! UNSTR_MESH_TOOLS
!!
!! This module provides a collection of procedures for generating higher-order
!! topological information for unstructured mixed-element meshes.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Adapted from existing code, February 2014
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module unstr_mesh_tools

  implicit none
  private

  public :: get_cell_neighbor_array, label_mesh_faces

  interface get_cell_neighbor_array
    procedure get_cell_neighbor_array_1, get_cell_neighbor_array_2
  end interface

contains

  subroutine get_cell_neighbor_array_1 (xcnode, cnode, xcnhbr, cnhbr, stat)

    use face_neighbor_table_type

    integer, intent(in) :: xcnode(:), cnode(:)
    integer, allocatable, intent(out) :: xcnhbr(:), cnhbr(:)
    integer, intent(out) :: stat

    type(face_neighbor_table) :: nhbr_table

    call nhbr_table%init (xcnode, cnode)
    call get_cnhbr_aux (nhbr_table, xcnode, cnode, xcnhbr, cnhbr, stat)

  end subroutine get_cell_neighbor_array_1

  subroutine get_cell_neighbor_array_2 (xcnode, cnode, xlnode, lnode, xcnhbr, cnhbr, lnhbr, stat)

    use face_neighbor_table_type
    use cell_topology, only: get_link_face_nodes

    integer, intent(in) :: xcnode(:), cnode(:)
    integer, intent(in) :: xlnode(:), lnode(:)
    integer, allocatable, intent(out) :: xcnhbr(:), cnhbr(:), lnhbr(:,:)
    integer, intent(out) :: stat

    integer :: i, j, k, nlink
    integer, allocatable :: face(:)
    type(face_neighbor_table) :: nhbr_table
    type(face_neighbor), allocatable :: nhbrs(:)

    call nhbr_table%init (xcnode, cnode)
    call get_cnhbr_aux (nhbr_table, xcnode, cnode, xcnhbr, cnhbr, stat)
    if (stat /= 0) return

    !! Lookup the two link faces in the table to get the neighbor cell numbers.
    !! The face node list must be arranged to match that stored in the table.
    nlink = size(xlnode) - 1
    allocate(lnhbr(2,nlink))
    do j = 1, nlink
      associate (link => lnode(xlnode(j):xlnode(j+1)-1))
        do k = 1, 2
          call get_link_face_nodes (link, k, face, normalize=.true.)
          call nhbr_table%get_neighbors (face, nhbrs)
          lnhbr(k,j) = 0  ! default value
          do i = size(nhbrs), 1, -1
            if (nhbrs(i)%p == 1) exit
          end do
          if (i > 0) then  ! found a match ...
            if (cnhbr(nhbrs(i)%k) == 0) then ! it is a boundary face as expected
              lnhbr(k,j) = nhbrs(i)%j
            else  ! not a boundary face; something is wrong
              stat = -1
            end if
          else  ! no match found; something is wrong
            stat = -1
          end if
        end do
      end associate
    end do

  end subroutine get_cell_neighbor_array_2

  subroutine get_cnhbr_aux (nhbr_table, xcnode, cnode, xcnhbr, cnhbr, stat)

    use face_neighbor_table_type
    use cell_topology, only: get_face_nodes

    type(face_neighbor_table), intent(in) :: nhbr_table
    integer, intent(in) :: xcnode(:), cnode(:)
    integer, allocatable, intent(out) :: xcnhbr(:), cnhbr(:)
    integer, intent(out) :: stat

    integer :: i, j, k, jj, kk, nmatch, bad_faces, ncell, offset
    integer, allocatable :: face(:)
    type(face_neighbor), allocatable :: nhbrs(:)

    ncell = size(xcnode) - 1

    !! Generate XCNHBR: CNHBR(XCNHBR(J):XCNHBR(J+1)-1) will store the face
    !! neighbors of cell J.  We infer the cell type from the number of nodes.
    allocate(xcnhbr(size(xcnode)))
    xcnhbr(1) = 1
    do j = 1, ncell
      select case (xcnode(j+1)-xcnode(j))
      case (4)  ! tet - 4 faces
        xcnhbr(j+1) = xcnhbr(j) + 4
      case (5)  ! pyramid - 5 faces
        xcnhbr(j+1) = xcnhbr(j) + 5
      case (6)  ! wedge - 5 faces
        xcnhbr(j+1) = xcnhbr(j) + 5
      case (8)  ! hex - 6 faces
        xcnhbr(j+1) = xcnhbr(j) + 6
      case default
        INSIST(.false.)
      end select
    end do

    allocate(cnhbr(xcnhbr(ncell+1)-1))

    cnhbr = 0
    bad_faces = 0
    stat = 0

    do j = 1, ncell
      associate (cell => cnode(xcnode(j):xcnode(j+1)-1))
        offset = xcnhbr(j)-1  ! for getting local face index
        do k = xcnhbr(j), xcnhbr(j+1)-1
          if (cnhbr(k) /= 0) cycle  ! info already assigned
          !! Get a face and the list of its neighbor cells.
          call get_face_nodes (cell, k-offset, face, normalize=.true.)
          call nhbr_table%get_neighbors (face, nhbrs)
          !! Locate the cell neighbor, but scan all for valid topology.
          jj = 0
          kk = 0
          nmatch = 0
          do i = 1, size(nhbrs)
            if (nhbrs(i)%j == j) cycle  ! found myself
            nmatch = nmatch + 1
            if (nhbrs(i)%p == -1) then  ! a good match (if only one)
              jj = nhbrs(i)%j
              kk = nhbrs(i)%k
            end if
          end do
          !! Store the neighbor information.
          if (nmatch == 1 .and. jj /= 0) then
            !! Found a unique neighbor; assign the neighbor data.
            cnhbr(k) = jj ! my neighbor, and
            cnhbr(kk) = j ! my neighbor's neighbor (me!)
          else if (nmatch /= 0) then
            !! Bad mesh topology; tag the faces involved.
            bad_faces = 1 + bad_faces
            do i = 1, size(nhbrs)
              cnhbr(nhbrs(i)%k) = -bad_faces
            end do
            stat = -1
          end if
        end do
      end associate
    end do

  end subroutine get_cnhbr_aux


  subroutine label_mesh_faces (xcnode, cnode, xlnode, lnode, nface, xcface, cface, lface)

    use cell_topology
    use facet_table_type

    integer, intent(in)  :: xcnode(:), cnode(:)
    integer, intent(in)  :: xlnode(:), lnode(:)
    integer, intent(out) :: nface
    integer, allocatable, intent(out) :: xcface(:), cface(:), lface(:,:)

    integer :: j, k, n, offset, ncell, nlink, max_face, node_max
    integer, allocatable :: face(:)
    type(facet_table) :: table

    ncell = size(xcnode) - 1

    !! Generate XCFACE: CFACE(XCFACE(J):XCFACE(J+1)-1) will store the face
    !! numbers of cell J.  We infer the cell type from the number of nodes.
    allocate(xcface(size(xcnode)))
    xcface(1) = 1
    do j = 1, ncell
      select case (xcnode(j+1)-xcnode(j))
      case (4)  ! tet - 4 faces
        xcface(j+1) = xcface(j) + 4
      case (5)  ! pyramid - 5 faces
        xcface(j+1) = xcface(j) + 5
      case (6)  ! wedge - 5 faces
        xcface(j+1) = xcface(j) + 5
      case (8)  ! hex - 6 faces
        xcface(j+1) = xcface(j) + 6
      case default
        INSIST(.false.)
      end select
    end do
    allocate(cface(xcface(ncell+1)-1))

    ASSERT(minval(cnode) > 0)

    max_face = size(cface)  ! worst case; realistically, closer to half this
    node_max = maxval(cnode)
    call table%init (max_face, node_max)

    do j = 1, ncell
      associate (cell => cnode(xcnode(j):xcnode(j+1)-1))
        offset = xcface(j)-1  ! for getting local face index
        do k = xcface(j), xcface(j+1)-1
          call get_face_nodes (cell, k-offset, face)
          call table%get_facet_label (face, cface(k), insert=.true.)
        end do
      end associate
    end do
    nface = table%number_of_facets()
    INSIST(all(cface /= 0))

    nlink = size(xlnode) - 1
    allocate(lface(2,nlink))
    do j = 1, nlink
      associate (link => lnode(xlnode(j):xlnode(j+1)-1))
        do k = 1, 2
          call get_link_face_nodes (link, k, face)
          call table%get_facet_label (face, lface(k,j), insert=.false.)
        end do
      end associate
    end do
    INSIST(all(lface >= 1))
    INSIST(all(lface <= nface))

  end subroutine label_mesh_faces

end module unstr_mesh_tools
