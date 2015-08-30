!!
!! UNSTR_MESH_TOOLS
!!
!! This module provides a collection of procedures for generating higher-order
!! topological information for unstructured mixed-element meshes.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Adapted from existing code, February 2014
!!

#include "f90_assert.fpp"

module unstr_mesh_tools

  implicit none
  private

  public :: get_cell_neighbor_array, label_mesh_faces

contains

  subroutine get_cell_neighbor_array (xcnode, cnode, xlnode, lnode, xcnhbr, cnhbr, lnhbr, stat)

    use facet_hash_type
    use cell_topology, only: get_face_nodes, get_link_face_nodes, facet_parity

    integer, intent(in) :: xcnode(:), cnode(:)
    integer, intent(in) :: xlnode(:), lnode(:)
    integer, allocatable, intent(out) :: xcnhbr(:), cnhbr(:), lnhbr(:,:)
    integer, intent(out) :: stat

    integer :: i, j, k, n, jj, kk, max_bin_size, nmatch, bad_faces, ncell, nlink, offset
    integer, allocatable :: face(:), xbin(:), p(:)
    type(facet_hash) :: hpar

    type :: table_entry
      integer :: j, k ! cell and face indices
      integer, allocatable :: face(:) ! face node list (outward and normalized)
    end type table_entry
    type(table_entry), allocatable :: bin_table(:)

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

    n = size(cnhbr)
    allocate(bin_table(n))

    !! Set-up the hash function.  It will return an address (or bin number)
    !! in the interval [0, N-1], where N is adjusted upward to a power of 2.
    !! With an ideal hash function, the number of bins is between 1/2 and 1
    !! times the number of faces, but is generally much closer to the
    !! lower bound.  Setting N to the number of faces is generous.
    call hpar%init (n, maxval(cnode))

    allocate(xbin(0:n))

    !! Count the number of hits to each bin; count for bin N stored in XBIN(N+1).
    xbin = 0
    do j = 1, ncell
      associate (cell => cnode(xcnode(j):xcnode(j+1)-1))
        offset = xcnhbr(j)-1  ! for getting local face index
        do k = xcnhbr(j), xcnhbr(j+1)-1
          call get_face_nodes (cell, k-offset, face)
          call hpar%hash (face, n)
          xbin(n+1) = 1 + xbin(n+1)
        end do
      end associate
    end do
    max_bin_size = maxval(xbin)

    !! Prepare XBIN: bin J will be BIN_TABLE(XBIN(J):XBIN(J+1)-1)
    xbin(0) = 1
    do j = 1, ubound(xbin,1)
      xbin(j) = xbin(j-1) + xbin(j)
    end do

    !! Fill the bin table; use XBIN as a temporary to hold the next free
    !! location for each bin.
    do j = 1, ncell
      associate (cell => cnode(xcnode(j):xcnode(j+1)-1))
        offset = xcnhbr(j)-1  ! for getting local face index
        do k = xcnhbr(j), xcnhbr(j+1)-1
          call get_face_nodes (cell, k-offset, face, normalize=.true.)
          call hpar%hash (face, n)
          i = xbin(n)
          bin_table(i)%j = j
          bin_table(i)%k = k
          call move_alloc (face, bin_table(i)%face)
          xbin(n) = i + 1
        end do
      end associate
    end do

    !! Restore XBIN: the index of the first element of bin J is now XBIN(J-1)
    !! instead of XBIN(J) as it should be -- fix this.
    do j = ubound(xbin,1), 1, -1
      xbin(j) = xbin(j-1)
    end do
    xbin(0) = 1

    allocate(p(max_bin_size))

    cnhbr = 0
    bad_faces = 0
    stat = 0

    do j = 1, ncell
      associate (cell => cnode(xcnode(j):xcnode(j+1)-1))
        offset = xcnhbr(j)-1  ! for getting local face index
        do k = xcnhbr(j), xcnhbr(j+1)-1
          if (cnhbr(k) /= 0) cycle  ! info already assigned
          !! Get a face and its corresponding bin.
          call get_face_nodes (cell, k-offset, face, normalize=.true.)
          call hpar%hash (face, n)
          associate (bin => bin_table(xbin(n):xbin(n+1)-1))
            !! Scan bin for *all* matching faces.
            jj = 0
            kk = 0
            nmatch = 0
            do i = 1, size(bin)
              if (bin(i)%j == j) then ! found myself
                p(i) = 1
              else
                p(i) = facet_parity(bin(i)%face, face)
                select case (p(i))
                case (-1) ! a good match (if only one)
                  nmatch = 1 + nmatch
                  jj = bin(i)%j
                  kk = bin(i)%k
                case (1)  ! a bad match (wrong orientation)
                  nmatch = 1 + nmatch
                end select
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
              do i = 1, size(bin)
                if (p(i) /= 0) cnhbr(bin(i)%k) = -bad_faces
              end do
              stat = -1
            end if
          end associate
        end do
      end associate
    end do

    !! Lookup the two link faces in the table to get the neighbor cell numbers.
    !! The face node list must be arranged to match that stored in the table.
    nlink = size(xlnode) - 1
    allocate(lnhbr(2,nlink))
    do j = 1, nlink
      associate (link => lnode(xlnode(j):xlnode(j+1)-1))
        do k = 1, 2
          call get_link_face_nodes (link, k, face, normalize=.true.)
          call hpar%hash (face, n)
          associate (bin => bin_table(xbin(n):xbin(n+1)-1))
            lnhbr(k,j) = 0  ! default value
            do i = size(bin), 1, -1
              if (all(bin(i)%face == face)) exit
            end do
            if (i > 0) then  ! found a match ...
              if (cnhbr(bin(i)%k) == 0) then ! and it's a boundary face as expected
                lnhbr(k,j) = bin(i)%j
              else  ! not a boundary face; something is wrong
                stat = -1
              end if
            else  ! no match found; something is wrong
              stat = -1
            end if
          end associate
        end do
      end associate
    end do

  end subroutine get_cell_neighbor_array


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
