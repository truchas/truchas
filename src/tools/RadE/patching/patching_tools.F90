!TODO: Finish documentation
!!
!! PATCHING_TOOLS
!!
!! This module provides utility functions used by the various patching
!! algorithms.
!!
!! David Neill-Asanza <dhna@lanl.gov>
!! 9 May 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!


#include "f90_assert.fpp"

module patching_tools

  use kinds, only: r8
  implicit none

contains


  ! TODO: merge with 3D code (unstr_mesh_tools?)
  subroutine get_face_neighbor_array (xface, fnode, xfnhbr, fnhbr, stat, normals, max_angle)

    use edge_neighbor_table_type

    integer, intent(in) :: xface(:), fnode(:)
    integer, allocatable, intent(out) :: xfnhbr(:), fnhbr(:)
    integer, intent(out) :: stat
    real(r8), intent(in), optional :: normals(:,:), max_angle ! max_angle is in radians

    type(edge_neighbor_table) :: nhbr_table

    call nhbr_table%init (xface, fnode)
    call get_fnhbr_aux (nhbr_table, xface, fnode, xfnhbr, fnhbr, stat, normals, max_angle)

  end subroutine get_face_neighbor_array


  ! TODO: merge with 3D code (unstr_mesh_tools?)
  subroutine get_fnhbr_aux (nhbr_table, xface, fnode, xfnhbr, fnhbr, stat, normals, max_angle)

    use edge_neighbor_table_type
    use cell_topology, only: get_edge_nodes

    type(edge_neighbor_table), intent(in) :: nhbr_table
    integer, intent(in) :: xface(:), fnode(:)
    integer, allocatable, intent(out) :: xfnhbr(:), fnhbr(:)
    integer, intent(out) :: stat
    real(r8), intent(in), optional :: normals(:,:), max_angle

    integer :: i, j, k, jj, kk, nmatch, bad_edges, nface, offset
    integer, allocatable :: edge(:)
    real(r8) :: angle
    type(edge_neighbor), allocatable :: nhbrs(:)

    nface = size(xface) - 1

    !! Generate XFNHBR: FNHBR(XFNHBR(J):XFNHBR(J+1)-1) will store the edge
    !! neighbors of face J.
    allocate(xfnhbr(size(xface)))
    xfnhbr(1) = 1
    do j = 1, nface
      ! Each face has one edge per node
      xfnhbr(j+1) = xfnhbr(j) + (xface(j+1)-xface(j))
    end do

    allocate(fnhbr(xfnhbr(nface+1)-1))

    fnhbr = 0
    bad_edges = 0
    stat = 0

    do j = 1, nface
      associate (face => fnode(xface(j):xface(j+1)-1))
        offset = xfnhbr(j)-1  ! for getting local edge index
        do k = xfnhbr(j), xfnhbr(j+1)-1
          if (fnhbr(k) /= 0) cycle  ! info already assigned
          !! Get a edge and the list of its neighbor faces.
          call get_edge_nodes (face, k-offset, edge, normalize=.true.)
          call nhbr_table%get_neighbors (edge, nhbrs)
          !! Locate the face neighbor, but scan all for valid topology.
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
            !! Ignore faces with large angles
            if (present(max_angle)) then
              angle = dot_product(normals(:,j), normals(:,jj))
              if (angle > 1) angle = 1  ! Fix floating point errors
              angle = acos( angle )
              if (angle >= max_angle) cycle
            end if
            !! Found a unique neighbor; assign the neighbor data.
            fnhbr(k) = jj ! my neighbor, and
            fnhbr(kk) = j ! my neighbor's neighbor (me!)
          else if (nmatch /= 0) then
            !! Bad mesh topology; tag the edges involved.
            bad_edges = 1 + bad_edges
            do i = 1, size(nhbrs)
              fnhbr(nhbrs(i)%k) = -bad_edges
            end do
            stat = -1
          end if
        end do
      end associate
    end do

  end subroutine get_fnhbr_aux

  ! TODO: this is implemented in 'truchas/utilities/graph_type.F90'.
  !! Return the connected components of the enclosure dual graph
  subroutine get_connected_faces(nface, xfnhbr, fnhbr, ncomp, xcomp, comp)

    integer, intent(in) :: nface
    integer, allocatable, intent(in) :: xfnhbr(:), fnhbr(:)
    integer, intent(out) :: ncomp
    integer, allocatable, intent(out) :: xcomp(:), comp(:)

    integer :: tag(nface)
    integer :: i

    ncomp = 0
    tag = 0

    !! Tag faces with component ID
    do i = 1, nface
      if (tag(i) /= 0) cycle
      ncomp = ncomp + 1
      call tag_component (i)
    end do
    allocate(xcomp(ncomp+1),comp(nface))

    !! Prepare XCOMP; nodes in component N will be COMP(XCOMP(N):XCOMP(N+1)-1)
    xcomp = 0
    do i = 1, nface
      xcomp(1+tag(i)) = xcomp(1+tag(i)) + 1
    end do
    xcomp(1) = 1
    do i = 2, size(xcomp)
      xcomp(i) = xcomp(i) + xcomp(i-1)
    end do

    !! Fill COMP; XCOMP(N) stores the next free location for component N.
    do i = 1, nface
      comp(xcomp(tag(i))) = i
      xcomp(tag(i)) = xcomp(tag(i)) + 1
    end do

    !! Restore XCOMP; XCOMP(N) is now the start of component N+1 instead of N
    do i = size(xcomp), 2, -1
      xcomp(i) = xcomp(i-1)
    end do
    xcomp(1) = 1
  contains
    !! Tag all the nodes connected to ROOT with the current component number
    recursive subroutine tag_component (root)
      integer, intent(in) :: root
      integer :: k, f
      tag(root) = ncomp
      do k = xfnhbr(root), xfnhbr(root+1)-1
        f = fnhbr(k)
        ! Skip missing neighbors (i.e. f is on the mesh boundary)
        if ( f <= 0) cycle
        if (tag(f) == 0) call tag_component (f)
      end do
    end subroutine
  end subroutine get_connected_faces


end module patching_tools
