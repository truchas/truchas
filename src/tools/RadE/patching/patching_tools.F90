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

#include "f90_assert.fpp"

module patching_tools

  use,intrinsic :: iso_fortran_env, only: r8 => real64, i8 => int64
  implicit none

  real(r8), parameter :: PI = 3.1415926535897932_r8

contains

  !TODO: refactor to more general module?
  subroutine init_random_seed(seed)

    integer, intent(in) :: seed

    integer, allocatable :: put(:)
    integer :: n, i
    integer(i8) :: s

    call random_seed(size=n)
    allocate(put(n))
    s = seed
    do i = 1, n
      put(i) = lcg(s)
    end do
    call random_seed(put=put)

    contains
      !! SOURCE: https://gcc.gnu.org/onlinedocs/gcc-4.9.1/gfortran/RANDOM_005fSEED.html
      !! This simple PRNG is seeded by a single integer.  This PRNG is used to
      !! generate a list of numbers with "high entropy" that will serve as the
      !! actual seed for RANDOM_SEED.
      function lcg(s)
        integer :: lcg
        integer(i8) :: s
        if (s == 0) then
           s = 104729
        else
           s = mod(s, 4294967296_i8)
        end if
        s = mod(s * 279470273_i8, 4294967291_i8)
        lcg = int(mod(s, int(huge(0), i8)), kind(0))
      end function lcg
  end subroutine init_random_seed


  ! TODO: merge with 3D code (unstr_mesh_tools?)
  subroutine get_face_neighbor_array(xface, fnode, xfnhbr, fnhbr, stat, errmsg, normals, max_angle)

    use edge_neighbor_table_type

    integer, intent(in) :: xface(:), fnode(:)
    integer, allocatable, intent(out) :: xfnhbr(:), fnhbr(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    real(r8), intent(in), optional :: normals(:,:), max_angle  ! max_angle is in degrees

    type(edge_neighbor_table) :: nhbr_table

    call nhbr_table%init(xface, fnode)
    call get_fnhbr_aux(nhbr_table, xface, fnode, xfnhbr, fnhbr, stat, errmsg, normals, max_angle)

  end subroutine get_face_neighbor_array


  ! TODO: merge with 3D code (unstr_mesh_tools?)
  subroutine get_fnhbr_aux(nhbr_table, xface, fnode, xfnhbr, fnhbr, stat, errmsg, normals, max_angle)

    use edge_neighbor_table_type
    use cell_topology, only: get_edge_nodes

    type(edge_neighbor_table), intent(in) :: nhbr_table
    integer, intent(in) :: xface(:), fnode(:)
    integer, allocatable, intent(out) :: xfnhbr(:), fnhbr(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    real(r8), intent(in), optional :: normals(:,:), max_angle  ! max_angle is in degrees

    integer :: i, j, k, jj, kk, nmatch, bad_edges, nface, offset
    integer, allocatable :: edge(:)
    real(r8) :: angle, max_angle_rad
    type(edge_neighbor), allocatable :: nhbrs(:)

    nface = size(xface) - 1

    if (present(max_angle)) then
      max_angle_rad = PI*max_angle/180.0_r8
    end if

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
          call get_edge_nodes(face, k-offset, edge, normalize=.true.)
          call nhbr_table%get_neighbors(edge, nhbrs)
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
              if (angle >= max_angle_rad) cycle
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

    if (stat/=0) errmsg = 'bad mesh topology: more than two faces share an edge'

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
      call tag_component(i)
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
    recursive subroutine tag_component(root)
      integer, intent(in) :: root
      integer :: k, f
      tag(root) = ncomp
      do k = xfnhbr(root), xfnhbr(root+1)-1
        f = fnhbr(k)
        ! Skip missing neighbors (i.e. f is on the mesh boundary)
        if ( f <= 0) cycle
        if (tag(f) == 0) call tag_component(f)
      end do
    end subroutine
  end subroutine get_connected_faces


  !! Finds the connected components of a subset of faces of the enclosure dual graph
  subroutine get_connected_faces_subset(xfnhbr, fnhbr, faces, tag, ncomp)

    integer, allocatable, intent(in) :: xfnhbr(:), fnhbr(:)
    integer, intent(in) :: faces(:)
    integer, allocatable, intent(out) :: tag(:)
    integer, intent(out) :: ncomp

    integer :: i

    allocate(tag(size(faces)))
    ncomp = 0
    tag = 0

    !! Find connected components
    do i = 1, size(faces)
      if (tag(i) /= 0) cycle
      ncomp = ncomp + 1
      call tag_component(i)
    end do

  contains
    !! Tag all the faces connected to ROOT with the current component number
    recursive subroutine tag_component(root)
#ifdef NO_2008_FINDLOC
      use f08_intrinsics, only: findloc
#endif
      integer, intent(in) :: root
      integer :: k, f, n, nid
      tag(root) = ncomp
      f = faces(root)
      do k = xfnhbr(f), xfnhbr(f+1)-1
        n = fnhbr(k)
        !! Skip missing neighbors (i.e. f is on the mesh boundary)
        if ( n <= 0) cycle
        nid = findloc(faces, n, dim=1)
        !! Ignore neighbors not in face list
        if (nid == 0) cycle
        if (tag(nid) == 0) call tag_component(nid)
      end do
    end subroutine tag_component
  end subroutine get_connected_faces_subset


  !! Finds the faces of a node
  ! TODO: based on cell_neighboring_vertices in mesh_geom_type.F90
  function faces_neighboring_vertices(e) result(ret)

    use re_encl_type

    type(encl), intent(in) :: e
    integer, allocatable :: ret(:,:)

    integer :: f,k,nid,j(e%nnode)

    !! Get maximum number of faces attached to a node
    j = 0
    do f = 1, e%nface
      do k = e%xface(f), e%xface(f+1)-1
        nid = e%fnode(k)
        j(nid) = j(nid) + 1
      end do
    end do

    allocate(ret(maxval(j),e%nnode))
    ret = -1

    !! j(nid) is the next unwritten index of ret(:,nid)
    j = 1
    do f = 1,e%nface
      do k = e%xface(f), e%xface(f+1)-1
        nid = e%fnode(k)
        ret(j(nid),nid) = f

        j(nid) = j(nid) + 1
      end do
    end do

  end function faces_neighboring_vertices


end module patching_tools
