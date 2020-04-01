!!
!! EDGE_NEIGHBOR_TABLE_TYPE
!!
!! The module defines a structure for efficiently looking up the enclosure faces
!! that contain a given edge. The edge is specified by its ordered list of node
!! numbers, and the neighboring faces are identified by face number, the side
!! index of the edge corresponding to that face, and the orientation of the edge
!! relative to the face. This is intended to be a short-lived structure used to
!! generate ultimate enclosure data.
!!
!! David Neill-Asanza <dhna@lanl.gov>
!! 6 May 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  The module defines the EDGE_NEIGHBOR_TABLE derived type. It has the
!!  following type bound procedures.
!!
!!  INIT (XFACE, FNODE) initializes the neighbor table. The mixed element
!!    enclosure is defined by the integer arrays XFACE and FNODE:
!!    FNODE(XFACE(j):XFACE(j+1)-1) are the nodes defining face j. The convention
!!    for how the nodes and faces of a cell are ordered is defined by the
!!    CELL_TOPOLOGY module.
!!
!!  GET_NEIGHBORS (EDGE, NHBRS) returns the neighbors of the given EDGE in the
!!    rank-1 allocatable array NHBRS.  EDGE is a rank-1 integer array storing the
!!    node numbers of the edge vertices in order (the order defines the
!!    orientation of the edge).  NHBRS is of derived type EDGE_NEIGHBOR which is
!!    also defined by this module.  It is allocated to the proper size by this
!!    procedure.  The type has three public integer components %J is the face
!!    number of the neighbor, %K is the side number of the face that corresponds
!!    to the given edge, and %P is the parity of the edge: %P == 1 if the edge is
!!    oriented outward with respect to the cell, and %P == -1 if it is oriented
!!    inward.
!!    NB: The "side number" is the sequential index of the edge within the entire
!!    enclosure, and not the local side number within the face.  The edges of
!!    face 1 are numbered, in order, followed by those of face 2, etc.  This
!!    turns out to be more convenient for client code than the more intuitive
!!    local side number.
!!

#include "f90_assert.fpp"

! TODO: merge with 3D code (face_neighbor_table_type)
module edge_neighbor_table_type

  use edge_hash_type
  implicit none
  private

  type :: table_record
    integer :: j, k ! face and edge indices
    integer, allocatable :: edge(:) ! edge node list (outward and normalized)
  end type table_record

  type, public :: edge_neighbor_table
    private
    type(edge_hash) :: edge_hash
    integer, allocatable :: xbin(:)
    type(table_record), allocatable :: bin_table(:)
  contains
    procedure :: init
    procedure :: get_neighbors
  end type edge_neighbor_table

  type, public :: edge_neighbor
    integer :: j, k, p
  end type edge_neighbor

contains

  subroutine init (this, xface, fnode)

    use cell_topology, only: get_edge_nodes

    class(edge_neighbor_table), intent(out) :: this
    integer, intent(in) :: xface(:), fnode(:)

    integer :: i, j, k, n, nface, offset
    integer, allocatable :: ecount(:), edge(:)

    nface = size(xface) - 1

    !! Precompute the number of edges for each face.
    !! It must be the number of nodes of the face.
    allocate(ecount(nface))
    do j = 1, nface
      ecount(j) = xface(j+1)-xface(j)
    end do

    !! Set-up the hash function.  It will return an address (or bin number)
    !! in the interval [0, N-1], where N is adjusted upward to a power of 2.
    !! With an ideal hash function, the number of bins is between 1/2 and 1
    !! times the number of edges, but is generally much closer to the
    !! lower bound.  Setting N to the number of edges is generous.
    n = sum(ecount)
    allocate(this%bin_table(n))
    call this%edge_hash%init (n)  ! N is modified
    allocate(this%xbin(0:n))

    !! Count the number of hits to each bin; count for bin N stored in XBIN(N+1).
    this%xbin = 0
    do j = 1, nface
      associate (face => fnode(xface(j):xface(j+1)-1))
        do k = 1, ecount(j)
          call get_edge_nodes (face, k, edge)
          call this%edge_hash%hash (edge, n)
          this%xbin(n+1) = 1 + this%xbin(n+1)
        end do
      end associate
    end do

    !! Prepare XBIN: bin J will be BIN_TABLE(XBIN(J):XBIN(J+1)-1)
    this%xbin(0) = 1
    do j = 1, ubound(this%xbin,1)
      this%xbin(j) = this%xbin(j-1) + this%xbin(j)
    end do

    !! Fill the bin table; use XBIN as a temporary to hold the next free
    !! location for each bin.
    offset = 0
    do j = 1, nface
      associate (face => fnode(xface(j):xface(j+1)-1))
        do k = 1, ecount(j)
          call get_edge_nodes (face, k, edge, normalize=.true.)
          call this%edge_hash%hash (edge, n)
          i = this%xbin(n)
          this%bin_table(i)%j = j
          this%bin_table(i)%k = k + offset
          call move_alloc (edge, this%bin_table(i)%edge)
          this%xbin(n) = i + 1
        end do
        offset = offset + ecount(j)
      end associate
    end do

    !! Restore XBIN: the index of the first element of bin J is now XBIN(J-1)
    !! instead of XBIN(J) as it should be -- fix this.
    do j = ubound(this%xbin,1), 1, -1
      this%xbin(j) = this%xbin(j-1)
    end do
    this%xbin(0) = 1

  end subroutine init


  subroutine get_neighbors (this, edge, nhbrs)

    use cell_topology, only: facet_parity

    class(edge_neighbor_table), intent(in) :: this
    integer, intent(in) :: edge(:)
    type(edge_neighbor), allocatable, intent(out) :: nhbrs(:)

    integer :: i, n
    integer, allocatable :: p(:)

    call this%edge_hash%hash (edge, n)
    associate (bin => this%bin_table(this%xbin(n):this%xbin(n+1)-1))
      allocate(p(size(bin)))
      do i = 1, size(bin)
        p(i) = facet_parity(bin(i)%edge, edge)
      end do
      n = count(p /= 0)
      allocate(nhbrs(n))
      n = 0
      do i = 1, size(bin)
        if (p(i) /= 0) then
          n = n + 1
          nhbrs(n)%j = bin(i)%j
          nhbrs(n)%k = bin(i)%k
          nhbrs(n)%p = p(i)
        end if
      end do
    end associate

  end subroutine get_neighbors

end module edge_neighbor_table_type
