!!
!! FACE_NEIGHBOR_TABLE_TYPE
!!
!! The module defines a structure for efficiently looking up the mesh cells that
!! contain a given face.  The oriented face is specified by its ordered list of
!! node numbers, and the neighboring cells are identified by cell number, the
!! side index of the cell that corresponds to the face, and the orientation of
!! the face relative to the cell. This is intended to be a short-lived structure
!! used to generate ultimate mesh data.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!!
!! PROGRAMMING INTERFACE
!!
!!  The module defines the FACE_NEIGHBOR_TABLE derived type.  It has the
!!  following type bound procedures.
!!
!!  INIT (XCNODE, CNODE) initializes the neighbor table.  The mixed element
!!    mesh is defined by the rank-1 integer arrays XCNODE and CNODE:
!!    CNODE(XCNODE(j):XCNODE(j+1)-1) are the nodes defining cell j.  The cell
!!    type is inferred from the number of nodes, and the convention for how the
!!    nodes and faces of a cell are ordered is defined by the CELL_TOPOLOGY
!!    module.
!!
!!  GET_NEIGHBORS (FACE, NHBRS) returns the neighbors of the given FACE in the
!!    rank-1 allocatable array NHBRS.  FACE is a rank-1 integer array storing
!!    the node numbers of the face vertices in order (the order defines the
!!    orientation of the face).  NHBRS is of derived type FACE_NEIGHBOR which
!!    is also defined by this module.  It is allocated to the proper size by
!!    this procedure.  The type has three public integer components %J is the
!!    cell number of the neighbor, %K is the side number of the cell that is
!!    the given face, and %P is the parity of the face: %P == 1 if the face
!!    is oriented outward with respect to the cell, and %P == -1 if it is
!!    oriented inward.
!!    NB: The "side number" is the sequential index of the side within the
!!    entire mesh, and not the local side number within the cell.  The sides
!!    of cell 1 are numbered, in order, followed by those of cell 2, etc.
!!    This turns out to be more convenient for client code than the more
!!    intuitive local side number.
!!
!! IMPLEMENTATION NOTES
!!
!!  The data structure uses a hash table of the nodes of the sides of the
!!  cells to quickly lookup the neighbors of a given face.
!!
!!  It would require very few changes to make this data structure able to
!!  handle the more general case of facet neighbors: cells containing an
!!  oriented edge, faces containing an oriented edge, etc.  Even cells,
!!  faces, or edges containing a node, but in this case the parity data
!!  is meaningless because a node has no orientation.
!!
!!  An second type of INIT procedure applicable to the case of single
!!  element type meshes may be useful.
!!

#include "f90_assert.fpp"

module face_neighbor_table_type

  use facet_hash_type
  implicit none
  private

  type :: table_record
    integer :: j, k ! cell and face indices
    integer, allocatable :: face(:) ! face node list (outward and normalized)
  end type table_record

  type, public :: face_neighbor_table
    private
    type(facet_hash) :: face_hash
    integer, allocatable :: xbin(:)
    type(table_record), allocatable :: bin_table(:)
  contains
    procedure :: init
    procedure :: get_neighbors
  end type face_neighbor_table

  type, public :: face_neighbor
    integer :: j, k, p
  end type face_neighbor

contains

  subroutine init (this, xcnode, cnode)

    use cell_topology, only: get_face_nodes

    class(face_neighbor_table), intent(out) :: this
    integer, intent(in) :: xcnode(:), cnode(:)

    integer :: i, j, k, n, ncell, offset
    integer, allocatable :: fcount(:), face(:)

    ncell = size(xcnode) - 1

    !! Precompute the number of faces for each cell.
    !! We infer the cell type from its number of nodes.
    allocate(fcount(ncell))
    do j = 1, ncell
      select case (xcnode(j+1)-xcnode(j))
      case (4)  ! tet - 4 faces
        fcount(j) = 4
      case (5)  ! pyramid - 5 faces
        fcount(j) = 5
      case (6)  ! wedge - 5 faces
        fcount(j) = 5
      case (8)  ! hex - 6 faces
        fcount(j) = 6
      case default
        INSIST(.false.)
      end select
    end do

    !! Set-up the hash function.  It will return an address (or bin number)
    !! in the interval [0, N-1], where N is adjusted upward to a power of 2.
    !! With an ideal hash function, the number of bins is between 1/2 and 1
    !! times the number of faces, but is generally much closer to the
    !! lower bound.  Setting N to the number of faces is generous.
    n = sum(fcount)
    allocate(this%bin_table(n))
    call this%face_hash%init (n, maxval(cnode))  ! N is modified
    allocate(this%xbin(0:n))

    !! Count the number of hits to each bin; count for bin N stored in XBIN(N+1).
    this%xbin = 0
    do j = 1, ncell
      associate (cell => cnode(xcnode(j):xcnode(j+1)-1))
        do k = 1, fcount(j)
          call get_face_nodes (cell, k, face)
          call this%face_hash%hash (face, n)
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
    do j = 1, ncell
      associate (cell => cnode(xcnode(j):xcnode(j+1)-1))
        do k = 1, fcount(j)
          call get_face_nodes (cell, k, face, normalize=.true.)
          call this%face_hash%hash (face, n)
          i = this%xbin(n)
          this%bin_table(i)%j = j
          this%bin_table(i)%k = k + offset
          call move_alloc (face, this%bin_table(i)%face)
          this%xbin(n) = i + 1
        end do
        offset = offset + fcount(j)
      end associate
    end do

    !! Restore XBIN: the index of the first element of bin J is now XBIN(J-1)
    !! instead of XBIN(J) as it should be -- fix this.
    do j = ubound(this%xbin,1), 1, -1
      this%xbin(j) = this%xbin(j-1)
    end do
    this%xbin(0) = 1

  end subroutine init


  subroutine get_neighbors (this, face, nhbrs)

    use cell_topology, only: facet_parity

    class(face_neighbor_table), intent(in) :: this
    integer, intent(in) :: face(:)
    type(face_neighbor), allocatable, intent(out) :: nhbrs(:)

    integer :: i, n
    integer, allocatable :: p(:)

    call this%face_hash%hash (face, n)
    associate (bin => this%bin_table(this%xbin(n):this%xbin(n+1)-1))
      allocate(p(size(bin)))
      do i = 1, size(bin)
        p(i) = facet_parity(bin(i)%face, face)
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

end module face_neighbor_table_type
