!!
!! UNSTR_MESH_TOOLS
!!
!! This module provides a collection of procedures for generating higher-order
!! topological information for unstructured meshes.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Adapted from existing code, February 2014
!!

#include "f90_assert.fpp"

module unstr_mesh_tools

  implicit none
  private
  
  public :: get_cell_neighbor_array, label_mesh_faces, get_face_node_array

contains

  subroutine get_cell_neighbor_array (cnode, lnode, cnhbr, lnhbr, stat)
  
    use facet_hash_type
    use cell_topology

    integer, intent(in)  :: cnode(:,:)
    integer, intent(out) :: cnhbr(:,:)
    integer, intent(in)  :: lnode(:,:)
    integer, intent(out) :: lnhbr(:,:)
    integer, intent(out) :: stat

    integer :: i, j, k, n, jj, kk, max_bin_size, nmatch, bad_faces
    integer, allocatable :: xbin(:), p(:)
    integer, pointer :: f(:) => null()
    type(facet_hash) :: hpar

    type :: table_entry
      integer :: j, k                     ! cell and local face indices
      integer, pointer :: f(:) => null()  ! face node list (outward and normalized)
    end type table_entry
    type(table_entry), allocatable :: bin_table(:)

    procedure(hex_face_nodes), pointer :: face_nodes

    !! Infer the type of mesh we are working with from the shape of CNODE and CNHBR.
    ASSERT(size(cnode,dim=2) == size(cnhbr,dim=2))
    ASSERT(size(lnode,dim=2) == size(lnhbr,dim=2))
    ASSERT(size(lnhbr,dim=1) == 2)
    if (size(cnode,dim=1) == 4) then ! tetrahedral cells
      ASSERT(size(cnhbr,dim=1) == 4)
      ASSERT(size(lnode,dim=1) == 6)
      face_nodes => tet_face_nodes
    else if (size(cnode,dim=1) == 8) then ! hexahedral cells
      ASSERT(size(cnhbr,dim=1) == 6)
      ASSERT(size(lnode,dim=1) == 8)
      face_nodes => hex_face_nodes
    else
      INSIST(.false.)
    end if
    
    !!! Check that each cell is non-degenerate
    
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
    do j = 1, size(cnode,dim=2)
      do k = 1, size(cnhbr,dim=1)
        f => face_nodes(cnode(:,j), k)
        call hpar%hash (f, n)
        xbin(n+1) = 1 + xbin(n+1)
        deallocate(f)
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
        f => face_nodes(cnode(:,j), k, normalize=.true.)
        call hpar%hash (f, n)
        i = xbin(n)
        bin_table(i)%f => f
        bin_table(i)%j = j
        bin_table(i)%k = k
        xbin(n) = i + 1
      end do
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

    do j = 1, size(cnhbr,dim=2)
      do k = 1, size(cnhbr,dim=1)
      
        if (cnhbr(k,j) /= 0) cycle  ! info already assigned
        
        !! Get a face and its corresponding bin.
        f => face_nodes(cnode(:,j), k, normalize=.true.)
        call hpar%hash (f, n)
        associate (bin => bin_table(xbin(n):xbin(n+1)-1))
        
          !! Scan bin for *all* matching faces.
          jj = 0
          kk = 0
          nmatch = 0
          do i = 1, size(bin)
            if (bin(i)%j == j) then ! found myself
              p(i) = 1
            else
              p(i) = parity(bin(i)%f, f)
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
          deallocate(f)

          if (nmatch == 1 .and. jj /= 0) then
            !! Found a unique neighbor; assign the neighbor data.
            cnhbr(k,j) = jj   ! my neighbor, and
            cnhbr(kk,jj) = j  ! my neighbor's neighbor (me!)
          else if (nmatch /= 0) then
            !! Bad mesh topology; tag the faces involved.
            bad_faces = 1 + bad_faces
            do i = 1, size(bin)
              if (p(i) /= 0) cnhbr(bin(i)%j, bin(i)%k) = -bad_faces
            end do
            stat = -1
          end if
        end associate
        
      end do
    end do
    
    !! Lookup the two link faces in the table to get the neighbor cell numbers.
    !! The face node list must be arranged to match that stored in the table.
    do j = 1, size(lnode,dim=2)
      do k = 1, 2
        f => link_face_nodes(lnode(:,j), k, normalize=.true.)
        call hpar%hash (f, n)
        associate (bin => bin_table(xbin(n):xbin(n+1)-1))
          lnhbr(k,j) = 0  ! default value
          do i = size(bin), 1, -1
            if (all(bin(i)%f == f)) exit
          end do
          if (i > 0) then  ! found a match ...
            if (cnhbr(bin(i)%k,bin(i)%j) == 0) then ! and it's a boundary face as expected
              lnhbr(k,j) = bin(i)%j
            else  ! not a boundary face; something is wrong
              stat = -1
            end if
          else  ! no match found; something is wrong
            stat = -1
          end if
          deallocate(f)
        end associate
      end do
    end do
    
    !! Clean-up
    do j = lbound(bin_table,1), ubound(bin_table,1)
      if (associated(bin_table(j)%f)) deallocate(bin_table(j)%f)
    end do
    deallocate(p, xbin, bin_table)
    
  end subroutine get_cell_neighbor_array


  subroutine label_mesh_faces (cnode, nface, cface, lnode, lface)

    use cell_topology
    use facet_table_type

    integer, intent(in)  :: cnode(:,:)
    integer, intent(out) :: nface
    integer, intent(out) :: cface(:,:)
    integer, intent(in)  :: lnode(:,:)
    integer, intent(out) :: lface(:,:)

    integer :: j, k, n, max_face, node_max
    integer, pointer :: f(:) => null()
    type(facet_table) :: table

    procedure(hex_face_nodes), pointer :: face_nodes

    !! Infer the type of mesh we are working with from the shape of CNODE and CFACE.
    ASSERT(size(cnode,dim=2) == size(cface,dim=2))
    ASSERT(size(lnode,dim=2) == size(lface,dim=2))
    ASSERT(size(lface,dim=1) == 2)
    if (size(cnode,dim=1) == 4) then ! tetrahedral cells
      ASSERT(size(cface,dim=1) == 4)
      ASSERT(size(lnode,dim=1) == 6)
      face_nodes => tet_face_nodes
    else if (size(cnode,dim=1) == 8) then ! hexahedral cells
      ASSERT(size(cface,dim=1) == 6)
      ASSERT(size(lnode,dim=1) == 8)
      face_nodes => hex_face_nodes
    else
      INSIST(.false.)
    end if

    ASSERT(minval(cnode) > 0)

    max_face = product(shape(cface))  ! worst case; realistically, closer to half this
    node_max = maxval(cnode)
    call table%init (max_face, node_max)

    do j = 1, size(cface,dim=2)
      do k = 1, size(cface,dim=1)
        f => face_nodes(cnode(:,j), k)
        call table%get_facet_label (f, cface(k,j))
      end do
    end do
    nface = table%number_of_facets()
    INSIST(all(cface /= 0))

    do j = 1, size(lnode,dim=2)
      do k = 1, 2
        f => link_face_nodes(lnode(:,j), k, normalize=.true.)
        call table%get_facet_label2 (f, lface(k,j))
        deallocate(f)
      end do
    end do
    INSIST(all(lface >= 1))
    INSIST(all(lface <= nface))

  end subroutine label_mesh_faces

  !! Generate the face-node connectivity array FNODE.  The procedure assumes
  !! every face is cooriented with exactly one local cell face (which are
  !! outward-oriented).  The CFACE/CFPAR arrays produced by LABEL_MESH_FACES
  !! have this property.  The returned face-node lists are normalized (the
  !! smallest node begins the list).

  subroutine get_face_node_array (cnode, cface, cfpar, fnode)

    use cell_topology

    integer, intent(in)  :: cnode(:,:)  ! cell-node connectivity array
    integer, intent(in)  :: cface(:,:)  ! cell-face connectivity array
    integer, intent(in)  :: cfpar(:)    ! cell-face orientation parity bit mask
    integer, intent(out) :: fnode(:,:)  ! face-node connectivity array

    integer :: j, k

    procedure(get_hex_face_nodes), pointer :: get_face_nodes

    !! Infer the type of mesh we are working with from the shape of CNODE and CFACE.
    ASSERT(size(cface,dim=2) == size(cnode,dim=2))
    ASSERT(size(cface,dim=2) == size(cfpar))
    if (size(cnode,dim=1) == 4) then ! tetrahedral cells
      ASSERT(size(cface,dim=1) == 4)
      ASSERT(size(fnode,dim=1) == 3)
      get_face_nodes => get_tet_face_nodes
    else if (size(cnode,dim=1) == 8) then ! hexahedral cells
      ASSERT(size(cface,dim=1) == 6)
      ASSERT(size(fnode,dim=1) == 4)
      get_face_nodes => get_hex_face_nodes
    else
      INSIST(.false.)
    end if

    ASSERT(minval(cnode) > 0)
    ASSERT(minval(cface) > 0)
    ASSERT(maxval(cface) <= size(fnode,dim=2))

    fnode = 0
    do j = 1, size(cface,dim=2)
      do k = 1, size(cface,dim=1)
        if (fnode(1,cface(k,j)) == 0) then
          call get_face_nodes (fnode(:,cface(k,j)), cnode(:,j), k, normalize=.true.)
          if (btest(cfpar(j),pos=k)) call reverse_facet (fnode(:,cface(k,j)))
        end if
      end do
    end do

    ASSERT(minval(fnode) > 0)

  end subroutine get_face_node_array

end module unstr_mesh_tools
