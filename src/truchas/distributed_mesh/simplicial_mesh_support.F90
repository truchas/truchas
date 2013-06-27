!!
!! SIMPLICIAL_MESH_SUPPORT
!!
!! This module provides some basic procedures for fleshing-out the topological
!! structure of triangular and tetrahedral meshes that start as a bare list
!! of cells.  It also provides some numerically robust procedures for computing
!! geometrical quantities, such as lengths, areas, normals, etc.
!!
!! This module, originally named MESH_SUPPORT, was part of my existing EM solver
!! that was imported into Truchas.  The module has been evolving since the
!! mid-90's and has been used in my finite element and mimetic difference codes.
!!
!! Neil N. Carlson <nnc@newmexico.com>
!! Last revised 4 Mar 2004
!!
!! PROGRAMMING INTERFACE
!!
!!    (TBD)
!!
!! IMPLEMENTATION NOTES
!!
!!    (TBD)
!!

#include "f90_assert.fpp"

module simplicial_mesh_support

  use kinds, only: r8
  implicit none
  private

  public :: label_cell_edges, label_cell_faces
  public :: assemble_cell_neighbor_list
  public :: assemble_edge_node_list, assemble_face_node_list, assemble_face_edge_list
  public :: face_nodes
  
  public :: vector_length, normalize, edge_length, tri_area, signed_tri_area, tet_volume
  public :: face_normal, tet_face_normal
  public :: bc_coord
  
  integer, dimension(2,3), parameter, public :: TRI_EDGE_VERT = &
    reshape(source=(/ 2,3,  1,3,  1,2 /), shape=(/2,3/))

  integer, dimension(2,6), parameter, public :: TET_EDGE_VERT = &
    reshape(source=(/ 1,2,  1,3,  1,4,  2,3,  2,4,  3,4 /), shape=(/2,6/))
    
  integer, dimension(3,4), parameter, public :: TET_FACE_VERT = &
    reshape(source=(/ 2,3,4,  1,3,4,  1,2,4,  1,2,3 /), shape=(/3,4/))
    
  integer, dimension(3,4), parameter, public :: TET_FACE_EDGE = &
    reshape(source=(/ 6,5,4,  6,3,2,  5,3,1,  4,2,1 /), shape=(/3,4/))

  !! Don't mess with this!  The bit sizes are critical to the hash function.
  integer, parameter, private :: i4 = selected_int_kind(9)   ! 4-byte integer, < 2^31
  integer, parameter, private :: i8 = selected_int_kind(18)  ! 8-byte integer, < 2^63

  type, private :: HashTable
    integer(kind=i4) :: tlen = 0  ! Table length
    integer(kind=i4), allocatable :: value(:)   ! Table's key values
    integer(kind=i4), allocatable :: key(:,:)   ! Table's vector keys
    integer(kind=i4) :: n = 0   ! label counter to supply new key values
    !! Hash function parameters
    integer :: tbits  ! Number of table address bits
    integer :: wbits  ! Number of bits in hash multiplier
    integer :: h1bit  ! Start bit for extracting hash address
    integer :: h2bit  ! Start bit for extracting hash address increment
    integer(kind=i8) :: a   ! Hash function multiplier
    integer(kind=i8) :: wmask   ! Multiplier mask
    !! Performance counters
    integer :: nss = 0  ! Number of successful searches
    integer :: nus = 0  ! Number of unsuccessful searches
    integer :: nsp = 0  ! Number of probes in successful searches
    integer :: nup = 0  ! Number of probes in unsuccessful searches
    integer :: msp = 0  ! Maximum probes in a successful probe
    integer :: mup = 0  ! Maximum probes in an unsuccessful probe
  end type HashTable

  interface tri_area
    module procedure tri_area_l, tri_area_x
  end interface
  
  interface tet_face_normal
    module procedure tet_face_normal_one, tet_face_normal_all
  end interface
  
contains

  subroutine label_cell_edges (cnode, cedge, nedge)
  
    integer, intent(inout) :: cnode(:,:)
    integer, intent(out)   :: cedge(:,:), nedge
    
    integer :: j, k, nnode, ncell, tlen
    type(HashTable) :: table
    
    ncell = size(cnode,dim=2)
    nnode = maxval(cnode)
    
    ASSERT( size(cnode,dim=1) == 3 .or. size(cnode,dim=1) == 4 )
    ASSERT( size(cnode,dim=2) == size(cedge,dim=2) )
    ASSERT( minval(cnode) == 1 )
    
    !! 
    call sort_cell_node_lists (cnode)
    
    select case (size(cnode,dim=1))
    case (3)  ! Triangles
      ASSERT( size(cedge,dim=1) == 3 )
      tlen = 4*nnode  ! 2 cells/node, 3 edges/cell double counted -> ~3 faces/node
      call create_hash_table (table, tlen, klen=2, kmax=nnode)
      do j = 1, ncell
        do k = 1, 3
          call get_value (table, cnode(TRI_EDGE_VERT(:,k),j), cedge(k,j))
          if (cedge(k,j) == 0) then
            print *, 'label_cell_edges: PANIC!  hash table is full!'
            stop
          end if
        end do
      end do
      nedge = label_max(table)
      call delete_hash_table (table)
    
    case (4)  ! Tetrahedra
      ASSERT( size(cedge,dim=1) == 6 )
      tlen = 8*nnode  ! ~7 edges/node in subdivided hex mesh; 3 edge, 3 face, 1 interior edges
      call create_hash_table (table, tlen, klen=2, kmax=nnode)
      do j = 1, ncell
        do k = 1, 6
          call get_value (table, cnode(TET_EDGE_VERT(:,k),j), cedge(k,j))
          if (cedge(k,j) == 0) then
            print *, 'label_cell_edges: PANIC!  hash table is full!'
            stop
          end if
        end do
      end do
      nedge = label_max(table)
      call delete_hash_table (table)
    end select
    
    ASSERT( all(cedge > 0) )
    
  end subroutine label_cell_edges
  
  subroutine label_cell_faces (cnode, cface, nface)
  
    integer, intent(inout) :: cnode(:,:)
    integer, intent(out)   :: cface(:,:)
    integer, intent(out)   :: nface
    
    integer :: j, k, nnode, ncell, tlen
    type(HashTable) :: table
    
    nnode = maxval(cnode)
    ncell = size(cnode,dim=2)
    
    ASSERT( size(cnode,dim=1) == 4 )
    ASSERT( size(cnode,dim=2) == size(cface,dim=2) )
    ASSERT( size(cface,dim=1) == 4 )
    ASSERT( minval(cnode) == 1 )
    
    call sort_cell_node_lists (cnode)
    
    tlen = 14*nnode ! 6 cells/node, 4 faces/cell double counted -> ~12 faces/node
    call create_hash_table (table, tlen, klen=3, kmax=nnode)
    
    do j = 1, ncell
      do k = 1, size(cface,dim=1)
        call get_value (table, cnode(TET_FACE_VERT(:,k),j), cface(k,j))
        if (cface(k,j) == 0) then
          print *, 'label_cell_faces: PANIC!  hash table is full!'
          stop
        end if
      end do
    end do
    nface = label_max(table)
    call delete_hash_table (table)
    
    ASSERT( all(cface > 0) )
    
  end subroutine label_cell_faces
  
  subroutine sort_cell_node_lists (cnode)
    integer, intent(inout) :: cnode(:,:)
    integer :: i, j, k, next
    do j = 1, size(cnode,dim=2) ! simple insertion sort of each list
      do k = 2, size(cnode,dim=1)
        i = k
        next = cnode(k,j)
        do while (i > 1)  ! find location i of next value
          if (cnode(i-1,j) <= next) exit
          cnode(i,j) = cnode(i-1,j)
          i = i - 1
        end do
        if (i /= k) cnode(i,j) = next
      end do
    end do
  end subroutine sort_cell_node_lists
  
  subroutine assemble_edge_node_list (cedge, cnode, enode)
  
    integer, intent(in)  :: cedge(:,:)
    integer, intent(in)  :: cnode(:,:)
    integer, intent(out) :: enode(:,:)
    
    integer :: i, j, k
    
    ASSERT( size(cnode,dim=1) == 3 .or. size(cnode,dim=1) == 4 )
    ASSERT( size(cnode,dim=2) == size(cedge,dim=2) )
    ASSERT( size(enode,dim=1) == 2 )
    ASSERT( minval(cnode) > 0 )
    ! This assertion fails for zero-sized arrays.
    !ASSERT( minval(cedge) == 1 .and. maxval(cedge) == size(enode,dim=2) )

    select case (size(cnode,dim=1))
    case (3)  ! Triangles
      ASSERT( size(cedge,dim=1) == 3 )
      enode = 0
      do j = 1, size(cedge,dim=2)
        do k = 1, size(cedge,dim=1)
          i = cedge(k,j)
          if (enode(1,i) <= 0) enode(:,i) = cnode(TRI_EDGE_VERT(:,k),j)
        end do
      end do
    
    case (4)  ! Tetrahedra
      ASSERT( size(cedge,dim=1) == 6 )
      enode = 0
      do j = 1, size(cedge,dim=2)
        do k = 1, size(cedge,dim=1)
          i = cedge(k,j)
          if (enode(1,i) <= 0) enode(:,i) = cnode(TET_EDGE_VERT(:,k),j)
        end do
      end do
    end select
    
    ASSERT( minval(enode) > 0 )
    
  end subroutine assemble_edge_node_list
  
  subroutine assemble_face_node_list (cface, cnode, fnode)
  
    integer, intent(in)  :: cface(:,:)
    integer, intent(in)  :: cnode(:,:)
    integer, intent(out) :: fnode(:,:)
    
    integer :: i, j, k
    
    ASSERT( size(cface,dim=1) == 4 )
    ASSERT( size(cface,dim=2) == size(cnode,dim=2) )
    ASSERT( size(cnode,dim=1) == 4 )
    ASSERT( size(fnode,dim=1) == 3 )
    ASSERT( minval(cnode) > 0 )
    ! This assertion fails for zero-sized arrays.
    !ASSERT( minval(cface) == 1 .and. maxval(cface) == size(fnode,dim=2) )
    
    fnode = 0
    do j = 1, size(cface,dim=2)
      do k = 1, size(cface,dim=1)
        i = cface(k,j)
        if (fnode(1,i) <= 0) fnode(:,i) = cnode(TET_FACE_VERT(:,k),j)
      end do
    end do
    
    ASSERT( minval(fnode) > 0 )
    
  end subroutine assemble_face_node_list

  subroutine assemble_face_edge_list (cface, cedge, fedge)
  
    integer, intent(in)  :: cface(:,:)
    integer, intent(in)  :: cedge(:,:)
    integer, intent(out) :: fedge(:,:)
    
    integer :: i, j, k
    
    ASSERT( size(cface,dim=1) == 4 )
    ASSERT( size(cface,dim=2) == size(cedge,dim=2) )
    ASSERT( size(cedge,dim=1) == 6 )
    ASSERT( size(fedge,dim=1) == 3 )
    ASSERT( minval(cedge) > 0 )
    ! This assertion fails for zero-sized arrays.
    !ASSERT( minval(cface) == 1 .and. maxval(cface) == size(fedge,dim=2) )
    
    fedge = 0
    do j = 1, size(cface,dim=2)
      do k = 1, size(cface,dim=1)
        i = cface(k,j)
        if (fedge(1,i) <= 0) fedge(:,i) = cedge(TET_FACE_EDGE(:,k),j)
      end do
    end do
    
    ASSERT( minval(fedge) > 0 )
    
  end subroutine assemble_face_edge_list
  
  !! This is simply awful :-(
  function face_nodes(fedge, enode) result (fnode)
    integer, intent(in) :: fedge(:)
    integer, intent(in) :: enode(:,:)
    integer :: fnode(3)
    ASSERT( size(fedge) == 3 )
    ASSERT( size(enode,dim=1) == 2 )
    ASSERT( minval(fedge) > 0 )
    ASSERT( maxval(fedge) <= size(enode,dim=2) )
    fnode(1) = enode(1,fedge(3))
    fnode(2) = enode(2,fedge(3))
    fnode(3) = enode(2,fedge(2))
  end function face_nodes
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! ASSEMBLE_CELL_NEIGHBOR_LIST
 !!
 !! This routine is from a much older version of this module, and uses a
 !! different numbering of the triange edge vertices and tet face vertices.
 !! This is irrelevant to the final cell neighbor array; the triangle edges
 !! and tet faces have the same local numbering.  The routine assumes the
 !! cells are positively oriented, in that the standard volume determinant
 !! which depends on the vertex order, is positive.
 !!
  
    
  subroutine assemble_cell_neighbor_list (cnode, cnhbr, stat)

    integer, intent(in)  :: cnode(:,:)
    integer, intent(out) :: cnhbr(:,:)
    integer, intent(out), optional :: stat
    
    integer :: i, j, jj, k, kk, n, ncell, nvert, max_bin_size, nmatch
    integer, allocatable :: bnum(:), xbin(:), p(:), facevert(:,:)

    type cell_face
      integer :: cell, face
    end type cell_face
    type(cell_face), allocatable, target ::  bin_table(:)
    type(cell_face), pointer :: bin(:)
    
    integer, dimension(2,3), parameter :: TriFaceVert = &
      reshape(source=(/ 2,3,  3,1,  1,2 /), shape=(/2,3/))
    integer, dimension(3,4), parameter :: TetFaceVert = &
      reshape(source=(/ 2,3,4,  1,4,3,  1,2,4,  1,3,2 /), shape=(/3,4/))
                                                                                
    nvert = size(cnode,dim=1)
    ncell = size(cnode,dim=2)
    
    ASSERT( nvert >= 3 .and. nvert <= 4 )
    ASSERT( all(shape(cnode) == shape(cnhbr)) )
    ASSERT( minval(cnode) > 0 )  ! assumed in the allocation of XBIN

    stat = 0
    if (ncell == 0) return
    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! BIN TABLE CONSTRUCTION
  !!
  !! We construct a bin table to accelerate the matching of adjacent cells.  The hash function
  !! used is simply the sum of the node numbers of the face vertices.  The cell and local face
  !! numbers are stored rather the face itself.  To find the neighboring cell across a face, only
  !! the appropriate bin needs to be searched.
  !!
  !! NB: It should be easy to find a relatively simple hash function that is far superior,
  !! yielding smaller bins.
  !!
  !! NB: We go to the extra effort of binning every face of every cell so that we may in addition
  !! check for correct mesh topology: unique neighbor having the opposite face orientation.
  !!

    allocate(bnum(nvert), bin_table(size(cnode)), xbin((nvert-1)*maxval(cnode)))
    
    xbin = 0
    do j = 1, ncell
      bnum = sum(cnode(:,j)) - cnode(:,j) ! hash function
      do k = 1, nvert
        xbin(1+bnum(k)) = xbin(1+bnum(k)) + 1
      end do
    end do
    max_bin_size = maxval(xbin)

    !! XBIN(J+1) now contains the size of bin J.

    xbin(1) = 1
    do j = 2, size(xbin)
      xbin(j) = xbin(j) + xbin(j-1)
    end do

    !! XBIN(J) is now the BIN_TABLE index of the first element of bin J;
    !! XBIN(J+1)-1 is the index of the last element.

    do j = 1, ncell
      bnum = sum(cnode(:,j)) - cnode(:,j) ! hash function
      do k = 1, nvert
        n = xbin(bnum(k))
        bin_table(n) % cell = j
        bin_table(n) % face = k
        xbin(bnum(k)) = n + 1
      end do
    end do

    !! The bin table is now filled.  Now, however, XBIN(J) is the
    !! index of the first element of bin J+1 instead of bin J; fix this.

    do j = size(xbin), 2, -1
      xbin(j) = xbin(j-1)
    end do
    xbin(1) = 1
   
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! FIND MATCHING FACES AND CHECK FOR CORRECT MESH TOPOLOGY
  !!

    allocate(facevert(nvert-1,nvert), p(max_bin_size))
    
    select case (nvert)
    case (3)  ! triangular cells
      facevert = TriFaceVert
    case (4)  ! tetrahedral cells
      facevert = TetFaceVert
    end select
    
    cnhbr = 0

    do j = 1, ncell
      bnum = sum(cnode(:,j)) - cnode(:,j)   ! hash function
      do k = 1, nvert
        if (cnhbr(k,j) /= 0) cycle ! info already assigned
        bin => bin_table(xbin(bnum(k)):xbin(bnum(k)+1)-1)
        jj = 0
        nmatch = 0
        do i = 1, size(bin)
          if (bin(i) % cell == j) then  ! found myself
            p(i) = 1
          else
            p(i) = parity(cnode(facevert(:,k),j), cnode(facevert(:,bin(i)%face),bin(i)%cell))
            select case (p(i))
            case (-1)   ! a viable match (should only be one)
              nmatch = nmatch + 1
              jj = bin(i) % cell
              kk = bin(i) % face
            case (1)    ! a match but with the wrong orientation (should not occur)
              nmatch = nmatch + 1
            end select
          end if
        end do
        if (nmatch == 1 .and. jj /= 0) then ! found THE neighbor
          cnhbr(k,j) = jj    ! my neighbor, and
          cnhbr(kk,jj) = j   ! my neighbor's neighbor (me!)
        else if (nmatch > 0) then ! bad mesh topology; tag the cell faces involved
          do i = 1, size(bin)
            if (p(i) /= 0) cnhbr(bin(i)%face,bin(i)%cell) = -(1 + nmatch)
          end do
          if (present(stat)) stat = -1
        end if
      end do
    end do
    
    deallocate(p, facevert, xbin, bin_table, bnum)

  end subroutine assemble_cell_neighbor_list
  
  
  subroutine create_hash_table (table, tlen, klen, kmax)

    type(HashTable), intent(inout) :: table
    integer(kind=i4), intent(in) :: tlen ! Table length (minimum)
    integer(kind=i4), intent(in) :: klen ! Length of integer key vector
    integer(kind=i4), intent(in) :: kmax ! Maximum value of any key component

    integer(kind=i8)             :: ASRC    ! leading 63 bits of 1/GoldenRatio

    data ASRC/o'474335715627751237012'/

    ASSERT( tlen > 0 )
    ASSERT( klen > 0 )
    ASSERT( kmax > 0 )

    table%tbits = bit_width(tlen)                 ! Number of table address bits
    table%tlen  = ibset(0_i4, table%tbits)         ! Table length (power of 2)
    table%wbits = digits(ASRC) - bit_width(kmax)  ! Number of multiplier bits
    table%wmask = ibset(0_i8, table%wbits) - 1_i8 ! Hash product mask
    table%a     = ishft(ASRC, -bit_width(kmax))       ! Hash multiplier
    table % h1bit = table%wbits - table%tbits     ! start bit for extracting hash address
    table % h2bit = table%h1bit - table%tbits     ! start bit for extracting hash address increment
    table % h2bit = max(0_i4, table%h2bit)
    table % h2bit = 0 ! this seems to work better

    if (allocated(table%value)) then
      if (size(table%value) /= table%tlen) deallocate(table%value)
    end if

    if (allocated(table%key)) then
      if (any(shape(table%key) /= (/ klen, table%tlen /))) deallocate(table%key)
    end if

    if (.not.allocated(table%value)) allocate(table%value(0:table%tlen-1))
    if (.not.allocated(table%key)) allocate(table%key(klen,0:table%tlen-1))

    !! Initialize hash table
    table%value = 0
    table%key = 0
    table%n = 0

    !! Initialize performance counters
    table%nss = 0
    table%nus = 0
    table%nsp = 0
    table%nup = 0
    table%msp = 0
    table%mup = 0

  contains

    integer function bit_width (i)
      integer, intent(in) :: i
      bit_width = digits(i)
      do while (.not.btest(i,bit_width-1))
        bit_width = bit_width - 1
        if (bit_width == 0) exit
      end do
    end function bit_width
    
  end subroutine create_hash_table
  
  
  subroutine delete_hash_table (table)
  
    type(HashTable), intent(inout) :: table
    
    if (allocated(table % value)) deallocate(table % value)
    if (allocated(table % key))   deallocate(table % key)
    
    table%tlen = 0
    table%n = 0
    
    table%nss = 0
    table%nus = 0
    table%nsp = 0
    table%nup = 0
    table%msp = 0
    table%mup = 0
    
  end subroutine delete_hash_table
  
  
  integer function label_max (table)
    type(HashTable), intent(in) :: table
    label_max = table%n
  end function label_max
  
  real function fill_factor (table)
    type(HashTable), intent(in) :: table
    fill_factor = real(table%n) / real(table%tlen)
  end function fill_factor
  
  subroutine print_hash_performance (table)
    type(HashTable), intent(in) :: table
    write(unit=*,fmt='(/,a)') 'HASH TABLE PERFORMANCE'
    write(unit=*,fmt='(a)') '-----------------------------------------'
    write(unit=*,fmt=*) '                Fill factor:', real(table%n) / real(table%tlen)
    write(unit=*,fmt=*) '          Number of entries:', table%n
    write(unit=*,fmt=*) '    Avg probes, succ search:', real(table%nsp) / real(table%nss)
    write(unit=*,fmt=*) '  Avg probes, unsucc search:', real(table%nup) / real(table%nus)
    write(unit=*,fmt=*) '    Max probes, succ search:', table%msp
    write(unit=*,fmt=*) '  Max probes, unsucc search:', table%mup
    write(unit=*,fmt=*) 'Succ to unsucc search ratio:', real(table%nss) / real(table%nus)
  end subroutine print_hash_performance
  
    
  subroutine get_value (table, key, value)

    type(HashTable), intent(inout) :: table
    integer(kind=i4), intent(in)   :: key(:)
    integer(kind=i4), intent(out)  :: value
    
    integer :: p, np
    integer(kind=i4) :: i, j
    
    ASSERT( size(key) == size(table%key,dim=1) )
    
    call hash (key, i, j)
    
    np = 0
    value = 0
    
    do
      if (table % value(i) == 0) then ! key not found; insert it
        if (table % n == table % tlen - 1) return ! table is full
        table%n = 1 + table%n     ! next label in sequence
        table%value(i) = table%n  ! set key value
        value = table%value(i)    ! set return value
        table%key(:,i) = key      ! set table key
        !! Update performance counters
        table%nus = table%nus + 1
        table%nup = table%nup + np
        table%mup = max(np, table%mup)
        exit
      else
        np = np + 1 ! Update the number of probes
        p = parity(key, table%key(:,i))  ! Compare keys for a match
        if (p /= 0) then  ! have a match; return value with relative orientation info
          value = p * table%value(i)
          !! Update performance counters
          table%nss = table%nss + 1
          table%nsp = table%nsp + np
          table%msp = max(np, table%msp)
          exit
        end if
        i = i - j   ! Next probe address; wrap to other end of table if necessary.
        if (i < 0) i = i + table%tlen
      end if
    end do
    
  contains
  
    subroutine hash (key, h1, h2)
    
      integer(i4), intent(in) :: key(:)
      integer(i4), intent(out) :: h1, h2
      
      integer :: j
      integer(kind=i8) :: p

      p = key(1) * table%a
      do j = 2, size(key)
        p = ieor(p, key(j) * table%a)
      end do
      p = iand(p, table%wmask)

      h1 = ibits(p, table%h1bit, table%tbits)
      h2 = ior(1_i8, ibits(p, table%h2bit, table%tbits))
      !h2 = 1 ! linear probing -- another disaster!
      !h2 = ior(1_i8, ibits(p, table%h2bit, 5)) ! a disaster!!!
      
      ASSERT( h1 >= 0 .and. h1 < table%tlen )
      ASSERT( h2 >= 1 .and. h2 < table%tlen )
      !print *, 'hash values:', h1, h2
      
    end subroutine hash
    
  end subroutine get_value
  
    
  integer function parity (key1, key2)

    integer, dimension(:), intent(in) :: key1, key2

    ASSERT( size(key1) == size(key2) )

    select case (size(key1))
    case (0)
      parity = 1

    case (1)
      if (key1(1) == key2(1)) then
        parity = 1
      else
        parity = 0
      end if

    case (2)
      parity = parity_pair(key1, key2)

    case (3)
      parity = parity_triple(key1, key2)

    case default
      parity = parity_general(key1, key2)
    end select

  end function parity

  function parity_general (key1, key2) result (p)
    integer, dimension(:), intent(in) :: key1, key2
    integer :: p
    integer, dimension(size(key1)) :: s
    integer :: j, m, n
    ASSERT( size(key1) == size(key2) )
    forall (j=1:size(s)) s(j) = j
    p = 1
    do n = size(key2), 1, -1
      m = n
      do while (key1(s(m)) /= key2(n))
        m = m - 1
        if (m == 0) exit
      end do
      if (m == 0) then
        p = 0
        exit
      else if (m /= n) then
        s(m) = s(n)
        p = -p
      end if
    end do

  end function parity_general

  function parity_pair (key1, key2) result (p)
    integer, dimension(:), intent(in) :: key1, key2
    integer :: p
    p = 0
    if (key1(2) == key2(2)) then
      if (key1(1) == key2(1)) p= 1
    else if (key1(1) == key2(2)) then
      if (key1(2) == key2(1)) p = -1
    end if
  end function parity_pair

  function parity_triple (key1, key2) result (p)
    integer, dimension(:), intent(in) :: key1, key2
    integer :: p
    ASSERT( size(key1) == size(key2) )
    p = 0
    if (key1(3) == key2(3)) then
      if (key1(2) == key2(2)) then
        if (key1(1) == key2(1)) p = 1
      else if (key1(1) == key2(2)) then
        if (key1(2) == key2(1)) p = -1
      end if
    else if (key1(2) == key2(3)) then
      if (key1(3) == key2(2)) then
        if (key1(1) == key2(1)) p = -1
      else if (key1(1) == key2(2)) then
        if (key1(3) == key2(1)) p = 1
      end if
    else if (key1(1) == key2(3)) then
      if (key1(2) == key2(2)) then
        if (key1(3) == key2(1)) p = -1
      else if (key1(3) == key2(2)) then
        if (key1(2) == key2(1)) p = 1
      end if
    end if
  end function parity_triple
  
  logical function valid_mesh (cnode)
  
    integer, intent(in) :: cnode(:,:)
    
    integer :: i, j, k, nvert, ncell, nnode, status
    integer, allocatable :: tag(:), cnhbr(:,:)
    
    nvert = size(cnode,dim=1)
    ncell = size(cnode,dim=2)
    
    valid_mesh = .false.
    
    !! Check for degenerate cells
    do j = 1, ncell
      do k = 1, nvert-1
        do i = k + 1, nvert
          if (cnode(i,j) == cnode(k,j)) return
        end do
      end do
    end do
    
    !! Check for unreferenced nodes
    nnode = maxval(cnode)
    if (minval(cnode) /= 1) return
    allocate(tag(nnode))
    tag = 0
    do j = 1, ncell
      do k = 1, nvert
        tag(cnode(k,j)) = 1
      end do
    end do
    if (any(tag == 0)) return
    deallocate(tag)
    
    !! Check for proper cell connections
    allocate(cnhbr(nvert,ncell))
    call assemble_cell_neighbor_list (cnode, cnhbr, status)
    deallocate(cnhbr)
    if (status /= 0) return
    
    valid_mesh = .true.
    
  end function valid_mesh
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  GEOMETRY FUNCTIONS
!!

  pure function vector_length (x) result (l)
  
    real(kind=r8), intent(in) :: x(:)
    real(kind=r8) :: l
    
    real(kind=r8) :: a, b, c, t
    
    select case (size(x))
    !! TWO-VECTOR CASE
    case (2)
      a = abs(x(1))
      b = abs(x(2))
      if (b > a) then ! reorder
        t = a
        a = b
        b = t
      end if
      l = a * sqrt(1.0_r8 + (b/a)**2)
    !! THREE-VECTOR CASE
    case (3)
      a = abs(x(1))
      b = abs(x(2))
      c = abs(x(3))
      if (b > a) then
        if (c > b) then ! swap largest value to A
          t = a
          a = c
          c = t
        else
          t = a
          a = b
          b = t
        end if
      else if (c > a) then
        t = a
        a = c
        c = t
      end if
      l = a * sqrt(1.0_r8 + ((b/a)**2 + (c/a)**2))
    !! FOR ANYTHING ELSE WE RETURN A BOGUS VALUE
    case default
      l = -huge(1.0_r8)
    end select
    
  end function vector_length
  
  
  pure function normalize (x) result (y)
  
    real(kind=r8), intent(in) :: x(:)
    real(kind=r8) :: y(size(x))
    
    integer :: n
    real(kind=r8) :: u, v, r
    
    select case (size(x))
    !! TWO-VECTOR CASE
    case (2)
      if (abs(x(2)) > abs(x(1))) then
        u = x(1) / abs(x(2))
        r = 1.0_r8 / sqrt(1.0_r8 + u**2)
        y(1) = u * r
        y(2) = sign(r, x(2))
      else
        u = x(2) / abs(x(1))
        r = 1.0_r8 / sqrt(1.0_r8 + u**2)
        y(1) = sign(r, x(1))
        y(2) = u * r
      end if
    !! THREE-VECTOR CASE
    case (3)
      n = 1
      if (abs(x(2)) > abs(x(n))) n = 2
      if (abs(x(3)) > abs(x(n))) n = 3
      select case (n)
      case (1)
        u = x(2) / abs(x(1))
        v = x(3) / abs(x(1))
        r = 1.0_r8 / sqrt(1.0_r8 + (u**2 + v**2))
        y(1) = sign(r, x(1))
        y(2) = u * r
        y(3) = v * r
      case (2)
        u = x(1) / abs(x(2))
        v = x(3) / abs(x(2))
        r = 1.0_r8 / sqrt(1.0_r8 + (u**2 + v**2))
        y(1) = u * r
        y(2) = sign(r, x(2))
        y(3) = v * r
      case (3)
        u = x(1) / abs(x(3))
        v = x(2) / abs(x(3))
        r = 1.0_r8 / sqrt(1.0_r8 + (u**2 + v**2))
        y(1) = u * r
        y(2) = v * r
        y(3) = sign(r, x(3))
      end select
    !! FOR ANYTHING ELSE WE RETURN A BOGUS VALUE
    case default
      y = -huge(1.0_r8)
    end select
    
  end function normalize
  
  
  pure function edge_length (x) result (l)
    real(kind=r8), intent(in) :: x(:,:)
    real(kind=r8) :: l
    l = vector_length (x(:,1)-x(:,2))
  end function edge_length
  

  pure function tri_area_l (l) result (area)
  
    real(kind=r8), intent(in) :: l(:)
    real(kind=r8) :: area, a, b, c, t
    
    a = l(1)
    b = l(2)
    c = l(3)
    
    ! Sort so that a >= b >= c
    if (b > a) then
      if (c > a) then
        t = a
        a = c
        c = t
        if (b > a) then
          t = a
          a = b
          b = t
        end if
      else
        t = a
        a = b
        b = t
      end if
    else
      if (c > b) then
        t = b
        b = c
        c = t
        if (b > a) then
          t = a
          a = b
          b = t
        end if
      end if
    end if
    
    if (c-(a-b) < 0.0_r8) then ! not the lengths of a real triangle
      area = 2.0_r8               ! trigger a floating point exception;
      area = sqrt(1.0_r8 - area)  ! must be a bit obtuse about it.
    else
      area = 0.25_r8 * sqrt((a+(b+c))*(c-(a-b))*(c+(a-b))*(a+(b-c)))
    end if

  end function tri_area_l
  
  
  pure function tri_area_x (x) result (area)
  
    real(kind=r8), intent(in) :: x(:,:)
    real(kind=r8) :: area, l(3)
    
    l(1) = vector_length(x(:,1)-x(:,2))
    l(2) = vector_length(x(:,2)-x(:,3))
    l(3) = vector_length(x(:,3)-x(:,1))
    area = tri_area_l (l)
    
  end function tri_area_x
  
  
  pure function signed_tri_area (x) result (area)
  
    real(kind=r8), intent(in) :: x(:,:)
    real(kind=r8) :: area, l(3), det
    
    !! This should be a plane triangle
    
    l(1) = vector_length(x(:,1)-x(:,2))
    l(2) = vector_length(x(:,2)-x(:,3))
    l(3) = vector_length(x(:,3)-x(:,1))
    area = tri_area_l (l)
    
    !! Naive method of computing (twice) the area; want only its sign.
    det = (x(1,2) - x(1,1))*(x(2,3) - x(2,1)) - (x(1,3) - x(1,1))*(x(2,2) - x(2,1))
    area = sign(area, det)
    area = 0.5_r8 * det
    
  end function signed_tri_area
  
  
  pure function tet_volume (x) result (v)   ! nothing sophisticated here
    real(kind=r8), intent(in) :: x(:,:)
    real(kind=r8) :: v
    v = ((x(1,2)-x(1,1))*((x(2,3)-x(2,1))*(x(3,4)-x(3,1)) - (x(2,4)-x(2,1))*(x(3,3)-x(3,1))) + &
         (x(2,2)-x(2,1))*((x(3,3)-x(3,1))*(x(1,4)-x(1,1)) - (x(3,4)-x(3,1))*(x(1,3)-x(1,1))) + &
         (x(3,2)-x(3,1))*((x(1,3)-x(1,1))*(x(2,4)-x(2,1)) - (x(1,4)-x(1,1))*(x(2,3)-x(2,1)))) &
      / 6.0_r8
  end function tet_volume
  
  
  pure function tet_face_normal_all (x) result (p)
    real(kind=r8), intent(in) :: x(:,:)
    real(kind=r8) :: p(3,4)
    p(:,1) = 0.5_r8 * cross_product(x(:,2) - x(:,4), x(:,3) - x(:,4))
    p(:,2) = 0.5_r8 * cross_product(x(:,3) - x(:,4), x(:,1) - x(:,4))
    p(:,3) = 0.5_r8 * cross_product(x(:,1) - x(:,4), x(:,2) - x(:,4))
    p(:,4) = 0.5_r8 * cross_product(x(:,3) - x(:,1), x(:,2) - x(:,1))
  end function tet_face_normal_all
  
  pure function tet_face_normal_one (x, k) result (p)
    real(kind=r8), intent(in) :: x(:,:)
    integer, intent(in) :: k
    real(kind=r8) :: p(3)
    select case (k)
    case (1)
      p = 0.5_r8 * cross_product(x(:,2) - x(:,4), x(:,3) - x(:,4))
    case (2)
      p = 0.5_r8 * cross_product(x(:,3) - x(:,4), x(:,1) - x(:,4))
    case (3)
      p = 0.5_r8 * cross_product(x(:,1) - x(:,4), x(:,2) - x(:,4))
    case (4)
      p = 0.5_r8 * cross_product(x(:,3) - x(:,1), x(:,2) - x(:,1))
    case default
      p = 0.0_r8
    end select
  end function tet_face_normal_one
    
  pure function face_normal (x) result (p)
    real(kind=r8), intent(in) :: x(:,:)
    real(kind=r8) :: p(3)
    p = 0.5_r8 * cross_product(x(:,2)-x(:,1), x(:,3)-x(:,2))
  end function face_normal
  
  pure function cross_product (a, b) result (axb)
    real(kind=r8), intent(in) :: a(:), b(:)
    real(kind=r8) :: axb(3)
    axb(1) = a(2) * b(3) - a(3) * b(2)
    axb(2) = a(3) * b(1) - a(1) * b(3)
    axb(3) = a(1) * b(2) - a(2) * b(1)
  end function cross_product

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! BARYCENTRIC COORDINATE PROCEDURES
!!
!! Given the vertex coordinates of a 2- or 3-simplex (i.e., triangle or tetrahedron) and a
!! point, this routine computes the barycentric coordinates of the point.  It is assumed
!! that the point lies in the simplex (if not the coordinates will not lie in [0,1]), and
!! that the 2-simplex is in R^2 or 3-simplex is in R^3.  For a simplex embedded in a higher
!! dimension space, (new) routines could be written that instead solve the normal equations
!! for the underdetermined linear system.
!!
  
  function bc_coord (x, p) result (lambda)
  
    real(kind=r8), intent(in) :: x(:,:) ! tetrahedron vertex coordinates
    real(kind=r8), intent(in) :: p(:)   ! point in tetrahedron
    real(kind=r8) :: lambda(size(p)+1)  ! barycentric coordinates of P
    
    integer :: n, j
    real(kind=r8) :: det, y(size(p),size(p)), q(size(p))
    
    n = size(p)
    
    ASSERT( n == 2 .or. n == 3 )
    ASSERT( size(x,dim=1) == n )
    ASSERT( size(x,dim=2) == n+1 )
  
    forall (j = 1:n) y(:,j) = x(:,j) - x(:,n+1)
    q = p - x(:,n+1)
    
    !! Mindless application of Cramer's rule...
    select case (n)
    case (2)
      det = y(1,1)*y(2,2) - y(2,1)*y(1,2)
      lambda(1) = (q(1) * y(2,2) - q(2) * y(1,2)) / det
      lambda(2) = (q(2) * y(1,1) - q(1) * y(2,1)) / det
      lambda(3) = 1.0_r8 - lambda(1) - lambda(2)
    
    case (3)
      det = y(1,1)*(y(2,2)*y(3,3) - y(3,2)*y(2,3)) + &
            y(2,1)*(y(3,2)*y(1,3) - y(1,2)*y(3,3)) + &
            y(3,1)*(y(1,2)*y(2,3) - y(2,2)*y(1,3))
      lambda(1) = (q(1) * (y(2,2)*y(3,3) - y(3,2)*y(2,3)) + &
                   q(2) * (y(3,2)*y(1,3) - y(1,2)*y(3,3)) + &
                   q(3) * (y(1,2)*y(2,3) - y(2,2)*y(1,3))) / det

      lambda(2) = (q(1) * (y(2,3)*y(3,1) - y(3,3)*y(2,1)) + &
                   q(2) * (y(3,3)*y(1,1) - y(1,3)*y(3,1)) + &
                   q(3) * (y(1,3)*y(2,1) - y(2,3)*y(1,1))) / det

      lambda(3) = (q(1) * (y(2,1)*y(3,2) - y(3,1)*y(2,2)) + &
                   q(2) * (y(3,1)*y(1,2) - y(1,1)*y(3,2)) + &
                   q(3) * (y(1,1)*y(2,2) - y(2,1)*y(1,2))) / det

      lambda(4) = 1.0_r8 - lambda(1) - lambda(2) - lambda(3)
    end select
    
    ASSERT( all(lambda >= 0.0_r8) .and. all(lambda <= 1.0_r8) )
    
  end function bc_coord

end module simplicial_mesh_support
