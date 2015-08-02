#include "f90_assert.fpp"

module simpl_mesh_tools

  implicit none
  private

  public :: get_cell_neighbor_array
  public :: label_mesh_faces
  public :: label_mesh_edges

contains

  subroutine get_tet_neighbor_array (cnode, cnhbr, stat)

    integer, intent(in)  :: cnode(:,:)
    integer, intent(out) :: cnhbr(:,:)
    integer, intent(out) :: stat
  
  end subroutine get_tet_neighbor_array

  subroutine get_cell_neighbor_array (cnode, get_face, cnhbr, stat)

    use facet_hash_type
    use cell_topology, only: parity

    integer, intent(in)  :: cnode(:,:)
    integer, intent(out) :: cnhbr(:,:)
    integer, intent(out) :: stat

    interface
      pure subroutine get_face (cell, face_index, face, normalize)
        integer, intent(in) :: cell(:), face_index
        integer, allocatable, intent(inout) :: face(:)
        logical, intent(in), optional :: normalize
      end subroutine
    end interface

    integer :: i, j, k, n, jj, kk, max_bin_size, nmatch, bad_faces
    integer, allocatable :: face(:), xbin(:), p(:)
    type(facet_hash) :: hpar

    type :: table_entry
      integer :: j, k                 ! cell and side indices
      integer, allocatable :: face(:) ! face node list (outward and normalized)
    end type table_entry
    type(table_entry), allocatable :: bin_table(:)

    ASSERT(size(cnode,dim=2) == size(cnhbr,dim=2))

    n = size(cnhbr)
    allocate(bin_table(n))

    !! Set-up the hash function.  It will return an address (or bin number)
    !! in the interval [0, N-1], where N is adjusted upward to a power of 2.
    !! With an ideal hash function, the number of bins is between 1/2 and 1
    !! times the number of faces, but is generally much closer to the
    !! lower bound.  Setting N to the number of faces is generous.
    call hpar%init (n, maxval(cnode))

    !! Count the number of hits to each bin; count for bin N stored in XBIN(N+1).
    allocate(xbin(0:n))
    xbin = 0
    do j = 1, size(cnode,dim=2)
      do k = 1, size(cnhbr,dim=1)
        call get_face (cnode(:,j), k, face)
        call hpar%hash (face, n)
        xbin(n+1) = 1 + xbin(n+1)
      end do
    end do
    max_bin_size = maxval(xbin)

    !! Prepare XBIN: bin J will be BIN_TABLE(XBIN(J):XBIN(J+1)-1)
    xbin(0) = 1
    do j = 1, ubound(xbin,1)
      xbin(j) = xbin(j-1) + xbin(j)
    end do

    !! Fill the bin table; use XBIN as a temporary to hold the next free
    !! location for each bin.
    do j = 1, size(cnode,dim=2)
      do k = 1, size(cnhbr,dim=1)
        call get_face (cnode(:,j), k, face, normalize=.true.)
        call hpar%hash (face, n)
        i = xbin(n)
        bin_table(i)%j = j
        bin_table(i)%k = k
        call move_alloc (face, bin_table(i)%face)
        xbin(n) = i + 1
      end do
    end do

    !! Restore XBIN: the index of the first element of bin J is now XBIN(J-1)
    !! instead of XBIN(J) as it should be -- fix this.
    do j = ubound(xbin,1), 1, -1
      xbin(j) = xbin(j-1)
    end do
    xbin(0) = 1

    cnhbr = 0
    bad_faces = 0
    stat = 0
    allocate(p(max_bin_size)) ! temp that stores parity info while searching a bin
    do j = 1, size(cnode,dim=2)
      do k = 1, size(cnhbr,dim=1)
        if (cnhbr(k,j) /= 0) cycle  ! info already assigned
        !! Get a face and its corresponding bin.
        call get_face (cnode(:,j), k, face, normalize=.true.)
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
              p(i) = parity(bin(i)%face, face)
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
            cnhbr(k,j) = jj   ! my neighbor, and
            cnhbr(kk,jj) = j  ! my neighbor's neighbor (me!)
          else if (nmatch /= 0) then
            !! Bad mesh topology; tag the faces involved.
            bad_faces = 1 + bad_faces
            do i = 1, size(bin)
              if (p(i) /= 0) cnhbr(bin(i)%k,bin(i)%j) = -bad_faces
            end do
            stat = -1
          end if
        end associate
      end do
    end do

  end subroutine get_cell_neighbor_array
  
  subroutine label_mesh_faces (cnode, cface, nface)
  
    use simplex_topology, only: TET_FACE_VERT
  
    integer, intent(in)  :: cnode(:,:)
    integer, intent(out) :: cface(:,:)
    integer, intent(out) :: nface
    
    integer :: max_face
    
    max_face = size(cface)  ! worst case; realistically much closer to half this
    call label_mesh_facets (cnode, TET_FACE_VERT, max_face, cface, nface)
    
  end subroutine label_mesh_faces
  
  subroutine label_mesh_edges (cnode, cedge, nedge)
  
    use simplex_topology, only: TET_EDGE_VERT
  
    integer, intent(in)  :: cnode(:,:)
    integer, intent(out) :: cedge(:,:)
    integer, intent(out) :: nedge
    
    integer :: max_edge
    
    max_edge = 2.5*size(cnode,2)
    call label_mesh_facets (cnode, TET_EDGE_VERT, max_edge, cedge, nedge)
    
  end subroutine label_mesh_edges
  
  subroutine label_mesh_facets (cnode, facet_vert, max_facet, cfacet, nfacet)
  
    use facet_table_type
    use simplex_topology, only:  TET_EDGE_VERT
  
    integer, intent(in)  :: cnode(:,:)
    integer, intent(in)  :: facet_vert(:,:)
    integer, intent(in)  :: max_facet
    integer, intent(out) :: cfacet(:,:)
    integer, intent(out) :: nfacet
    
    integer :: j, k, node_max
    type(facet_table) :: table
    
    ASSERT(size(cnode,dim=1) == 4)
    ASSERT(size(cnode,dim=2) == size(cfacet,dim=2))
    ASSERT(size(cfacet,dim=1) == size(facet_vert,dim=2))
    ASSERT(minval(cnode) > 0)
    
    node_max = maxval(cnode)
    call table%init (max_facet, node_max)
    
    do j = 1, size(cfacet,dim=2)
      do k = 1, size(cfacet,dim=1)
        call table%get_facet_label (cnode(facet_vert(:,k),j), cfacet(k,j), insert=.true.)
        ASSERT(cfacet(k,j) /= 0)
        ASSERT(cfacet(k,j) > 0)
      end do
    end do
    nfacet = table%number_of_facets()
    
    call table%write_hash_performance (0)
    
  end subroutine label_mesh_facets

end module simpl_mesh_tools
