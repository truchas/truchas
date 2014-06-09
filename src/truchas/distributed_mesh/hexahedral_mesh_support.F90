#include "f90_assert.fpp"

module hexahedral_mesh_support

  implicit none
  private
  
  public :: label_cell_edges, label_cell_faces
  public :: assemble_face_node_list
  public :: assemble_cell_neighbor_list
  public :: cell_neighbor_info
  public :: assemble_link_neighbor_list
  
contains

  subroutine assemble_cell_neighbor_list (cnode, cnhbr, stat)
  
    use hashing
    use cell_topology

    integer, intent(in)  :: cnode(:,:)
    integer, intent(out) :: cnhbr(:,:)
    integer, intent(out) :: stat

    integer :: i, j, k, n, jj, kk, max_bin_size, nmatch, bad_faces
    integer, allocatable :: xbin(:), p(:)
    integer, pointer :: f(:) => null()
    type(hash_param) :: hpar

    type :: table_entry
      integer :: j, k                     ! cell and local face indices
      integer, pointer :: f(:) => null()  ! face node list (outward and normalized)
    end type table_entry
    type(table_entry), allocatable, target :: bin_table(:)
    type(table_entry), pointer :: bin(:) => null()

    ASSERT( size(cnode,dim=2) == size(cnhbr,dim=2) )
    ASSERT( size(cnode,dim=1) == 8 )
    ASSERT( size(cnhbr,dim=1) == 6 )
    
    !!! Check that each cell is non-degenerate
    
    n = size(cnhbr)
    allocate(bin_table(n))

    !! Set-up the hash function.  It will return an address (or bin number)
    !! in the interval [0, N-1], where N is adjusted upward to a power of 2.
    !! With an ideal hash function, the number of bins is between 1/2 and 1
    !! times the number of faces, but is generally much closer to the
    !! lower bound.  Setting N to the number of faces is generous.
    call initialize_hash_param (hpar, n, maxval(cnode))

    allocate(xbin(0:n))

    !! Count the number of hits to each bin; count for bin N stored in XBIN(N+1).
    xbin = 0
    do j = 1, size(cnode,dim=2)
      do k = 1, size(cnhbr,dim=1)
        f => hex_face_nodes(cnode(:,j), k)
        call hash (hpar, f, n)
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
        f => hex_face_nodes(cnode(:,j), k, normalize=.true.)
        call hash (hpar, f, n)
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
        f => hex_face_nodes(cnode(:,j), k, normalize=.true.)
        call hash (hpar, f, n)
        bin => bin_table(xbin(n):xbin(n+1)-1)
        
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
        
      end do
    end do
    
    !! Clean-up
    do j = lbound(bin_table,1), ubound(bin_table,1)
      if (associated(bin_table(j)%f)) deallocate(bin_table(j)%f)
    end do
    deallocate(p, xbin, bin_table)
    
  end subroutine assemble_cell_neighbor_list
  

  subroutine cell_neighbor_info (cnode, lnode, cnhbr, lnhbr, stat)
  
    use hashing
    use cell_topology

    integer, intent(in)  :: cnode(:,:)
    integer, intent(out) :: cnhbr(:,:)
    integer, intent(in)  :: lnode(:,:)
    integer, intent(out) :: lnhbr(:,:)
    integer, intent(out) :: stat

    integer :: i, j, k, n, jj, kk, max_bin_size, nmatch, bad_faces
    integer, allocatable :: xbin(:), p(:)
    integer, pointer :: f(:) => null()
    type(hash_param) :: hpar

    type :: table_entry
      integer :: j, k                     ! cell and local face indices
      integer, pointer :: f(:) => null()  ! face node list (outward and normalized)
    end type table_entry
    type(table_entry), allocatable, target :: bin_table(:)
    type(table_entry), pointer :: bin(:) => null()

    ASSERT( size(cnode,dim=2) == size(cnhbr,dim=2) )
    ASSERT( size(cnode,dim=1) == 8 )
    ASSERT( size(cnhbr,dim=1) == 6 )
    ASSERT( size(lnode,dim=2) == size(lnhbr,dim=2) )
    ASSERT( size(lnode,dim=1) == 8 )
    ASSERT( size(lnhbr,dim=1) == 2 )
    
    !!! Check that each cell is non-degenerate
    
    n = size(cnhbr)
    allocate(bin_table(n))

    !! Set-up the hash function.  It will return an address (or bin number)
    !! in the interval [0, N-1], where N is adjusted upward to a power of 2.
    !! With an ideal hash function, the number of bins is between 1/2 and 1
    !! times the number of faces, but is generally much closer to the
    !! lower bound.  Setting N to the number of faces is generous.
    call initialize_hash_param (hpar, n, maxval(cnode))

    allocate(xbin(0:n))

    !! Count the number of hits to each bin; count for bin N stored in XBIN(N+1).
    xbin = 0
    do j = 1, size(cnode,dim=2)
      do k = 1, size(cnhbr,dim=1)
        f => hex_face_nodes(cnode(:,j), k)
        call hash (hpar, f, n)
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
        f => hex_face_nodes(cnode(:,j), k, normalize=.true.)
        call hash (hpar, f, n)
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
        f => hex_face_nodes(cnode(:,j), k, normalize=.true.)
        call hash (hpar, f, n)
        bin => bin_table(xbin(n):xbin(n+1)-1)
        
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
        
      end do
    end do
    
    !! Lookup the two link faces in the table to get the neighbor cell numbers.
    !! The face node list must be arranged to match that stored in the table.
    do j = 1, size(lnode,dim=2)
      do k = 1, 2
        !! The link faces are sides 5 and 6 of the 'hex' link cell LNODE(:,j).
        f => hex_face_nodes(lnode(:,j), 4+k, normalize=.true., reverse=.true.)
        call hash (hpar, f, n)
        bin => bin_table(xbin(n):xbin(n+1)-1)
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
      end do
    end do
    
    !! Clean-up
    do j = lbound(bin_table,1), ubound(bin_table,1)
      if (associated(bin_table(j)%f)) deallocate(bin_table(j)%f)
    end do
    deallocate(p, xbin, bin_table)
    
  end subroutine cell_neighbor_info
  
  subroutine assemble_link_neighbor_list (cnode, cnhbr, lnode, lnhbr)
  
    use hashing
    use cell_topology

    integer, intent(in)  :: cnode(:,:)
    integer, intent(in)  :: cnhbr(:,:)
    integer, intent(in)  :: lnode(:,:,:)
    integer, intent(out) :: lnhbr(:,:)
    
    integer :: j, k, tsize, n, inc
    integer, pointer :: f(:)
    type(hash_param) :: hpar
    
    type :: table_entry
      integer, pointer :: f(:) => null() ! the key
      integer :: j  ! the associated value
    end type
    type(table_entry), allocatable :: table(:)

    ASSERT(size(cnode,dim=2) == size(cnhbr,dim=2))
    ASSERT(size(cnode,dim=1) == 8)
    ASSERT(size(cnhbr,dim=1) == 6)
    ASSERT(size(lnode,dim=3) == size(lnhbr,dim=2))
    ASSERT(size(lnode,dim=1) == 4)
    ASSERT(size(lnode,dim=2) == 2)
    ASSERT(size(lnhbr,dim=1) == 2)
    
    !! Create a hash table containing all the boundary faces.
    !! The link faces must match one of these.
    tsize = 1.25 * count(cnhbr == 0)
    call initialize_hash_param (hpar, tsize, maxval(cnode))
    allocate(table(0:tsize-1))
    do j = 1, size(cnhbr,dim=2)
      do k = 1, size(cnhbr,dim=1)
        if (cnhbr(k,j) == 0) then ! this is a boundary face
          f => hex_face_nodes(cnode(:,j), k, normalize=.true.)
          call hash (hpar, f, n, inc)
          do while (associated(table(n)%f))
            n = n - inc
            if (n < 0) n = n + tsize
          end do
          table(n)%f => f
          table(n)%j = j
        end if
      end do
    end do
    
    !! Lookup the two link faces in the table to get the neighbor cell numbers.
    !! The face node list must be arranged to match that stored in the table.
    allocate(f(size(lnode,dim=1)))
    do j = 1, size(lnode,dim=3)
      do k = 1, 2
        f = lnode(:,k,j)
        call normalize_facet(f)
        if (k == 2) call reverse_facet(f) 
        call hash (hpar, f, n, inc)
        do while (associated(table(n)%f))
          if (all(f == table(n)%f)) exit
          n = n - inc
          if (n < 0) n = n + tsize
        end do
        ASSERT(associated(table(n)%f))
        lnhbr(k,j) = table(n)%j
      end do
    end do
    deallocate(f)
    
    !! Free the table.
    do j = 0, tsize-1
      if (associated(table(j)%f)) deallocate(table(j)%f)
    end do
    deallocate(table)
    
  end subroutine assemble_link_neighbor_list


  subroutine label_cell_edges (cnode, cedge, nedge)
  
    use facet_labeling
    use cell_topology
  
    integer, intent(in)  :: cnode(:,:)
    integer, intent(out) :: cedge(:,:)
    integer, intent(out) :: nedge
    
    integer :: j, k, max_edge, node_max
    integer, pointer :: f(:) => null()
    type(facet_table) :: table
    
    ASSERT( size(cnode,dim=1) == 8 )
    ASSERT( size(cnode,dim=2) == size(cedge,dim=2) )
    ASSERT( size(cedge,dim=1) == 12 )
    ASSERT( minval(cnode) == 1 )
    
    max_edge = 12 * size(cnode,dim=2)  ! worst case; realistically, something close to 3*ncell
    node_max = maxval(cnode)
    
    call create_facet_table (table, max_edge, node_max)
    
    do j = 1, size(cedge,dim=2)
      do k = 1, size(cedge,dim=1)
        f => hex_edge_nodes(cnode(:,j), k)
        call get_facet_label (table, f, cedge(k,j))
      end do
    end do
    nedge = number_of_facets (table)
    call destroy_facet_table (table)
    
    INSIST( all(cedge /= 0) )

  end subroutine label_cell_edges

  subroutine label_cell_faces (cnode, cface, nface, lnode, lface)
  
    use facet_labeling
    use cell_topology
  
    integer, intent(in)  :: cnode(:,:)
    integer, intent(out) :: cface(:,:)
    integer, intent(out) :: nface
    integer, intent(in)  :: lnode(:,:)
    integer, intent(out) :: lface(:,:)
    
    integer :: j, k, max_face, node_max
    integer, pointer :: f(:) => null()
    type(facet_table) :: table
    
    ASSERT( size(cnode,dim=1) == 8 )
    ASSERT( size(cnode,dim=2) == size(cface,dim=2) )
    ASSERT( size(cface,dim=1) == 6 )
    ASSERT( minval(cnode) > 0 )
    ASSERT( size(lnode,dim=1) == 8 )
    ASSERT( size(lnode,dim=2) == size(lface,dim=2) )
    ASSERT( size(lface,dim=1) == 2 )
    
    max_face = 6 * size(cnode,dim=2)  ! worst case; realistically, something close to 3*ncell
    node_max = maxval(cnode)
    
    call create_facet_table (table, max_face, node_max)

    do j = 1, size(cface,dim=2)
      do k = 1, size(cface,dim=1)
        f => hex_face_nodes(cnode(:,j), k)
        call get_facet_label (table, f, cface(k,j))
      end do
    end do
    nface = number_of_facets (table)
    INSIST( all(cface /= 0) )

    allocate(f(4))
    do j = 1, size(lnode,dim=2)
      f = lnode(1:4,j)
      call normalize_facet (f)
      call get_facet_label2 (table, f, lface(1,j))
      f = lnode(5:8,j)
      call normalize_facet (f)
      call reverse_facet (f)
      call get_facet_label2 (table, f, lface(2,j))
    end do
    deallocate(f)
    INSIST( all(lface >= 1) )
    INSIST( all(lface <= nface) )

    call destroy_facet_table (table)

  end subroutine label_cell_faces

  
  subroutine assemble_face_node_list (cface, cfpar, cnode, fnode)
  
    use cell_topology
    
    integer, intent(in)  :: cface(:,:)
    integer, intent(in)  :: cfpar(:)
    integer, intent(in)  :: cnode(:,:)
    integer, intent(out) :: fnode(:,:)
    
    integer :: j, k
    
    ASSERT( size(cface,dim=1) == 6 )
    ASSERT( size(cnode,dim=1) == 8 )
    ASSERT( size(cnode,dim=2) == size(cface,dim=2) )
    ASSERT( size(cfpar) == size(cface,dim=2) )
    ASSERT( size(fnode,dim=1) == 4 )
    ASSERT( minval(cnode) > 0 )
    ASSERT( minval(cface) >= 1 )
    ASSERT( maxval(cface) <= size(fnode,dim=2) )

    fnode = 0
    do j = 1, size(cface,dim=2)
      do k = 1, size(cface,dim=1)
        if (fnode(1,cface(k,j)) == 0) then
          fnode(:,cface(k,j)) =  cnode(HEX8_FACE_VERT(:,k),j)
          if (btest(cfpar(j),pos=k)) call reverse_facet (fnode(:,cface(k,j)))
        end if
      end do
    end do
    
    ASSERT( minval(fnode) > 0 )
    
  end subroutine assemble_face_node_list
  
!!  subroutine aux (ab, ac, map, bc)
!!  
!!    integer, intent(in) :: ab(:,:), ac(:,:), map(:,:)
!!    integer, intent(out) :: bc(:,:)
!!  
!!    ASSERT( size(ab,dim=1) == size(map,dim=2) )
!!    ASSERT( size(ab,dim=2) == size(ac,dim=2) )
!!    ASSERT( size(bc,dim=1) == size(map,dim=1) )
!!    ASSERT( minval(map) >= 1 .and. maxval(map) <= size(ac,dim=1) )
!!    ASSERT( minval(ab)  >= 1 .and. maxval(ab) <= size(bc,dim=2) )
!!    
!!    bc = 0
!!    do j = 1, size(ab,dim=2)
!!      do k = 1, size(ab,dim=1)
!!        bc(:,ab(k,j)) = ac(map(:,k),j)
!!      end do
!!    end do
!!
!!  end subroutine aux
  
end module hexahedral_mesh_support
