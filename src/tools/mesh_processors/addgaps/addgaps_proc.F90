!!
!! ADDGAPS_PROC
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! 27 Oct 2005
!!
!! This module provides the mesh modification subroutine ADD_GAP_ELEMENTS
!! used by the utility program addgaps.  This version of the module implements
!! a completely new algorithm for adding the gap elements that is more robust
!! and general than the algorithm used in the preceding version.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module addgaps_proc

  use exodus_mesh_type
  use string_utilities, only: i_to_c
  use command_line, only: prog, strict_ss_transf
  implicit none
  private

  public :: add_gap_elements

  integer, parameter, private :: LUN = 10

  !! Data structure for describing the active subset of the input mesh.
  type, private :: active_mesh
    integer :: num_elem
    integer, pointer :: enum(:) => null()
    integer, pointer :: bnum(:) => null()
    integer, pointer :: lnum(:) => null()
    integer, pointer :: nmask(:) => null()
    integer, pointer :: smask(:) => null()
    integer, pointer :: jnbr(:,:) => null()
    integer, pointer :: knbr(:,:) => null()
    integer, pointer :: gside(:,:) => null()
  end type active_mesh

  integer, parameter :: MAX_SIDE = 6  ! works for tet, wedge and hex element types.

  !! HASH FUNCTION OBJECTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Don't mess with these!  The bit sizes are critical to the hash function.
  integer, parameter :: i4 = selected_int_kind(9)   ! 4-byte integer, < 2^31
  integer, parameter :: i8 = selected_int_kind(18)  ! 8-byte integer, < 2^63

  !! Container for the hash function parameters.
  type, private :: hash_param
    integer :: hbits  ! Number of hash address bits
    integer :: wbits  ! Number of bits in hash multiplier
    integer :: h1bit  ! Start bit for extracting hash address
    integer(kind=i8) :: a   ! Hash function multiplier
    integer(kind=i8) :: wmask   ! Multiplier mask
  end type hash_param

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! ADD_GAP_ELEMENTS
 !!
 !! This subroutine adds gap elements to the input mesh INMESH returning the
 !! resulting mesh in OUTMESH.  The gap elements are inserted along the mesh
 !! surfaces described by the side sets whose IDs are given in the list SSID.
 !!
 !! The routine first identifies the active part of the input mesh, and
 !! generates the necessary data for it.  This data is then used to transform
 !! the input mesh to the output mesh containing gap elements.
 !!
 !! Some terminology used throughout the comments.  An active (or gap) node is
 !! a node on one of the gap surfaces.  An element is an active element if it
 !! contains an active node, and a side of an active element is an active side
 !! if the side contains an active node.  An active side is a gap side if it,
 !! or its neighbors corresponding side, was specified by one of the gap side
 !! sets.  The collection of all the active elements forms the active submesh;
 !! all other elements and all inactive sides and nodes can be ignored as they
 !! are not involved in the process of inserting gap elements.
 !!

  subroutine add_gap_elements (inmesh, ssid, outmesh)

    type(exodus_mesh), intent(in)  :: inmesh
    integer,           intent(in)  :: ssid(:)
    type(exodus_mesh), intent(out) :: outmesh

    integer :: i, j, k, n, b, l, bit_mask
    integer, pointer :: list(:), side_sig(:), xside(:), side(:,:), daughter(:)
    logical, allocatable :: active_sset(:), active_node(:), active_elem(:), active_eblk(:)
    type(active_mesh) :: active
    character(len=64) :: log_file

    integer, target, save :: TETRA4_SIDE_SIG(4)
    data TETRA4_SIDE_SIG/b'1011', b'1110', b'1101', b'0111'/

    integer, target, save :: WEDGE6_SIDE_SIG(5)
    data WEDGE6_SIDE_SIG/b'011011', b'110110', b'101101', b'000111', b'111000'/

    integer, target, save :: HEX8_SIDE_SIG(6)
    data HEX8_SIDE_SIG/b'00110011', b'01100110', b'11001100', b'10011001', b'00001111', b'11110000'/

    log_file = trim(prog) // '.log'
    open(LUN,file=trim(log_file),status='replace',action='write',position='rewind')

    !! Generate the active side set mask array ACTIVE_SSET.
    allocate(active_sset(inmesh%num_sset))
    active_sset = .false.
    do i = 1, size(ssid)
      n = inmesh%num_sset
      do while (n > 0)
        if (inmesh%sset(n)%ID == ssid(i)) exit
        n = n - 1
      end do
      if (n == 0) call halt ('no side set with ID ' // i_to_c(ssid(i)))
      active_sset(n) = .true.
    end do

    !! Generate the active node mask array ACTIVE_NODE.
    allocate(active_node(inmesh%num_node))
    active_node = .false.
    do n = 1, inmesh%num_sset
      if (.not.active_sset(n)) cycle
      list => inmesh%side_set_node_list(n)
      do i = 1, size(list)
        active_node(list(i)) = .true.
      end do
      deallocate(list)
    end do

    !! Generate the active element and element block mask arrays ACTIVE_ELEM and ACTIVE_EBLK.
    allocate(active_elem(inmesh%num_elem), active_eblk(inmesh%num_eblk))
    n = 0
    do b = 1, inmesh%num_eblk
      active_eblk(b) = .false.
      do l = 1, inmesh%eblk(b)%num_elem
        n = n + 1
        if (any(active_node(inmesh%eblk(b)%connect(:,l)))) then
          active_elem(n) = .true.
          active_eblk(b) = .true.
        else
          active_elem(n) = .false.
        end if
      end do
    end do

    !! Create the active mesh data structure.
    n = count(active_elem)
    allocate(active%enum(n), active%bnum(n), active%lnum(n), active%nmask(n), active%smask(n))
    active%num_elem = n
    n = 0
    j = 0
    do b = 1, inmesh%num_eblk
      if (active_eblk(b)) then  ! no active elements in this block.

        select case (inmesh%eblk(b)%elem_type)
        case ('TETRA', 'TETRA4')
          side_sig => TETRA4_SIDE_SIG
        case ('WEDGE', 'WEDGE6')
          side_sig => WEDGE6_SIDE_SIG
        case ('HEX', 'HEX8')
          side_sig => HEX8_SIDE_SIG
        case default
          call halt ('unable to handle element type: ' // inmesh%eblk(b)%elem_type)
        end select

        do l = 1, inmesh%eblk(b)%num_elem
          n = n + 1
          if (active_elem(n)) then
            j = j + 1

            active%enum(j) = n
            active%bnum(j) = b
            active%lnum(j) = l

            !! Define NMASK: mark the active nodes on this active element.
            associate(list => inmesh%eblk(b)%connect(:,l))
              bit_mask = 0
              do k = 1, size(list)
                if (active_node(list(k))) bit_mask = ibset(bit_mask, k-1)
              end do
              active%nmask(j) = bit_mask
            end associate

            !! Define SMASK: mark the active sides on this active element.
            bit_mask = 0
            do k = 1, size(side_sig)
              if (iand(active%nmask(j), side_sig(k)) /= 0) bit_mask = ibset(bit_mask, k-1)
            end do
            active%smask(j) = bit_mask
          end if
        end do

      else

        n = n + inmesh%eblk(b)%num_elem

      end if
    end do

    call generate_neighbor_data (inmesh, active)
    call generate_gside_data (inmesh, active_sset, active, xside, side)
    deallocate(active_sset, active_node, active_elem)

    !! Create the output mesh.
    outmesh%title = inmesh%title
    call define_eblk  (inmesh, active, outmesh, daughter) ! element block data
    call define_coord (inmesh, daughter, outmesh)         ! node coordinate data
    call define_nset  (inmesh, daughter, outmesh)         ! node set data
    call define_sset  (inmesh, xside, side, outmesh)      ! side set data
    deallocate(xside, side, daughter)

    deallocate(active%enum, active%bnum, active%lnum, active%nmask, active%smask)
    deallocate(active%jnbr, active%knbr, active%gside)

    close(LUN)
    write(*,'(a)') 'See file "' // trim(log_file) // '" for a summary of the mesh changes.'

  end subroutine add_gap_elements

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! GENERATE_NEIGHBOR_DATA
 !!
 !! This subroutine generates the neighbor data arrays JNBR and KNBR for the
 !! active subset of the input mesh.  When positive, JNBR(K,J) is the index of
 !! the active element that neighbors active element J across local side K, and
 !! KNBR(K,J) is the corresponding local side index of the neighbor.  Boundary
 !! sides are indicated by a zero JNBR value, and inactive/nonexistent sides
 !! to be ignored are indicated by a negative JNBR value.  These arrays
 !! satisfy the reciprocity relations
 !!
 !!      JJ = JNBR(K,J)>0,    J = JNBR(KK,JJ)>0,
 !!      KK = KNBR(K,J),      K = KNBR(KK,JJ).
 !!
 !! The identification of neighbors is accomplished by first binning up the
 !! active sides according to a hash of the nodes on a side.  Once this data
 !! structure has been created, neighbors can be found by a quick search of
 !! a (very) small bin.
 !!
 !! We implicitly assume that the incoming mesh has a valid topology and do
 !! not make the effort here that we could/should to verify a good topology.
 !!

  subroutine generate_neighbor_data (mesh, active)

    type(exodus_mesh), intent(in)    :: mesh
    type(active_mesh), intent(inout) :: active

    integer :: i, j, k, n, jj, kk
    integer, pointer :: xbin(:), key(:)
    type(hash_param) :: hpar

    type :: table_entry
      integer :: j, k
      integer, pointer :: key(:) => null()
    end type table_entry
    type(table_entry), pointer :: bin_table(:), bin(:)

    !! Count the number of active sides N.
    n = 0
    do j = 1, active%num_elem
      n = n + bit_count(active%smask(j))
    end do

    allocate(bin_table(n))

    !! Set-up the hash function.  It will return an address (or bin number)
    !! in the interval [0, N-1], where N is adjusted upward to a power of 2.
    !! With an ideal hash function, the number of bins is between 1/2 and 1
    !! times the number of active sides, but is generally much closer to the
    !! lower bound.  Setting N to the number of active sides is generous.
    call initialize_hash_param (hpar, n, mesh%num_node)

    allocate(xbin(0:n))

    !! Count the number of hits to each bin; store the count for bin N in XBIN(N+1).
    xbin = 0
    do j = 1, active%num_elem
      do k = 1, MAX_SIDE
        if (btest(active%smask(j),k-1)) then
          key => mesh%side_node_list(active%enum(j), k)
          n = hash(hpar, key)
          xbin(n+1) = xbin(n+1) + 1
          deallocate(key)
        end if
      end do
    end do

    !! Prepare XBIN: bin J will be BIN_TABLE(XBIN(J):XBIN(J+1)-1)
    xbin(0) = 1
    do j = 1, ubound(xbin,1)
      xbin(j) = xbin(j) + xbin(j-1)
    end do

    !! Fill the bin table; use XBIN as a temporary to hold the next free
    !! location for each bin.
    do j = 1, active%num_elem
      do k = 1, MAX_SIDE
        if (btest(active%smask(j),k-1)) then
          key => mesh%side_node_list(active%enum(j), k)
          n = hash(hpar, key)
          i = xbin(n)
          bin_table(i)%key => key
          bin_table(i)%j = j
          bin_table(i)%k = k
          xbin(n) = i + 1
        end if
      end do
    end do

    !! Restore XBIN: the index of the first element of bin J is now XBIN(J-1)
    !! instead of XBIN(J) as it should be -- fix this.
    do j = ubound(xbin,1), 1, -1
      xbin(j) = xbin(j-1)
    end do
    xbin(0) = 1

    !! Now find neighbors, filling KNBR and JNBR.
    allocate(active%jnbr(MAX_SIDE,active%num_elem), active%knbr(MAX_SIDE,active%num_elem))
    active%jnbr = -1
    active%knbr = -1
    do j = 1, active%num_elem
      do k = 1, MAX_SIDE
        if (btest(active%smask(j),k-1)) then

          if (active%jnbr(k,j) > 0) cycle  ! info already assigned

          !! The key my neighbor would have for this side.
          key => mesh%side_node_list(active%enum(j), k, reverse=.true.)
          n = hash(hpar, key)

          bin => bin_table(xbin(n):xbin(n+1)-1)

          kk = 0
          jj = 0
          do i = 1, size(bin)
            if (bin(i)%j == j) cycle  ! found myself
            if (size(bin(i)%key) /= size(key)) cycle  ! different type of face
            if (all(key == bin(i)%key)) then  ! found my neighbor
              kk = bin(i)%k
              jj = bin(i)%j
              exit
            end if
          end do

          active%knbr(k,j) = kk
          active%jnbr(k,j) = jj
          if (jj > 0) then  ! my neighbor's neighbor (me)
            active%knbr(kk,jj) = k
            active%jnbr(kk,jj) = j
          end if

        end if
      end do
    end do

    !! Clean-up.
    do j = 1, size(bin_table)
      deallocate(bin_table(j)%key)
    end do
    deallocate(bin_table, xbin)

  end subroutine generate_neighbor_data

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! Hash Function Procedures
 !!
 !! Variable-length key Fibonacci hash function.  See section 6.4 of Knuth's
 !! The Art of Computing Programming, Vol 3.  Note that we require the hash
 !! to be invariant with respect to permutations of the key components.  Thus
 !! we reduce a variable length key to a scalar integer before hashing using
 !! exclusive-or which is commutative and associative.
 !!

  integer function hash (hpar, key)

    type(hash_param), intent(in) :: hpar
    integer, intent(in) :: key(:)

    integer :: j
    integer(kind=i8) :: p

    p = key(1)*hpar%a
    do j = 2, size(key)
      p = ieor(p, key(j)*hpar%a)
    end do
    p = iand(p, hpar%wmask)

    !! Hash address in range [0,2**hpar%hbits-1]
    hash = ibits(p, hpar%h1bit, hpar%hbits)

  end function hash

  subroutine initialize_hash_param (hpar, hsize, kmax)

    type(hash_param), intent(out)   :: hpar
    integer(kind=i4), intent(inout) :: hsize ! Hash address size (minimum)
    integer(kind=i4), intent(in)    :: kmax  ! Max value of any key component

    integer(kind=i8) :: ASRC    ! leading 63 bits of 1/GoldenRatio
    data ASRC/o'474335715627751237012'/

    ASSERT( hsize > 0 )
    ASSERT( kmax > 0 )

    hpar%hbits = bit_width(hsize)                ! Number of hash address bits
    hpar%wbits = digits(ASRC) - bit_width(kmax)  ! Number of multiplier bits
    hpar%wmask = ibset(0_i8, hpar%wbits) - 1_i8  ! Hash product mask
    hpar%a     = ishft(ASRC, -bit_width(kmax))   ! Hash multiplier
    hpar%h1bit = hpar%wbits - hpar%hbits         ! start bit for extracting hash address

    !! Return the hash address size (power of 2)
    hsize = ibset(0_i4, hpar%hbits)

  contains

    integer function bit_width (i)
      integer(kind=i4), intent(in) :: i
      bit_width = digits(i)
      do while (.not.btest(i,bit_width-1))
        bit_width = bit_width - 1
        if (bit_width == 0) exit
      end do
    end function bit_width

  end subroutine initialize_hash_param

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! GENERATE_GSIDE_DATA
 !!
 !! This subroutine generates the GSIDE array which identifies the gap sides
 !! and other active sides of active elements.  When positive, GSIDE(K,J) is
 !! the index of the side set that specified side K of active element J as a
 !! gap side.  GSIDE = 0 for remaining active, non-gap sides.  GSIDE is
 !! negative for inactive or nonexistent sides.  The array satisfies the
 !! reciprocity relation
 !!
 !!   GSIDE(K,J) = GSIDE(KK,JJ) when KK = KNBR(K,J), JJ = JNBR(K,J) > 0.
 !!
 !! This means that adjacent elements have the same GSIDE value for their
 !! shared side.
 !!
 !! Each 'side' in a side set is an element/local-side pair.  The side sets
 !! that specify the gap surface may be one-sided (elements from one side of
 !! the surface), two-sided (elements from both sides of the surface), or
 !! something in between.  It is desirable in the output mesh, however, for
 !! a side set to include 'sides' from both sides of the gap it specified.
 !! This subroutine returns the XSIDE/SIDE arrays that contain any sides
 !! that were missing from the input mesh side set and that should be added
 !! to the output mesh side set.  SIDE(I,XSIDE(N):XSIDE(N+1)-1) are the
 !! missing elements (I=1) and corresponding local sides (I=2) that should
 !! be added to side set N in the output mesh.  Only when a gap side set
 !! is two-sided will there be no missing sides.
 !!

  subroutine generate_gside_data (mesh, active_sset, active, xside, side)

    type(exodus_mesh), intent(in) :: mesh
    logical, intent(in) :: active_sset(:)
    type(active_mesh), intent(inout) :: active
    integer, pointer :: xside(:), side(:,:)

    integer :: i, j, k, n, jj, kk, nn, map(mesh%num_elem)

    allocate(active%gside(MAX_SIDE,active%num_elem))
    active%gside = 0

    !! Generate the active-element-to-element number map.
    map = 0
    do j = 1, active%num_elem
      map(active%enum(j)) = j
    end do

    !! Tag the side-set-specified gap sides with the side set index.
    do n = 1, mesh%num_sset
      if (.not.active_sset(n)) cycle
      do i = 1, mesh%sset(n)%num_side
        k = mesh%sset(n)%face(i)
        j = map(mesh%sset(n)%elem(i))
        if (active%gside(k,j) == 0) then
          active%gside(k,j) = n
        else if (active%gside(k,j) /= n) then
          call halt('overlapping side sets: ID=' // i_to_c(mesh%sset(n)%ID) // ',' // i_to_c(mesh%sset(active%gside(k,j))%ID))
        end if
      end do
    end do

    !! Examine the GSIDE array for consistency.  Gap sides cannot be on the
    !! boundary of the mesh.  If tagged, the neighboring side to a gap side
    !! must be tagged with the same side set index.  If the neighboring side
    !! is not tagged, we tag it with the negative of the side set index to
    !! mark it as a side that needs to be added to the list of new sides for
    !! the side set.  We will unnegative the value later to symmetrize GSIDE.

    allocate(xside(1+mesh%num_sset))
    xside = 0
    do j = 1, active%num_elem
      do k = 1, MAX_SIDE
        n = active%gside(k,j)
        if (n > 0) then
          kk = active%knbr(k,j)
          jj = active%jnbr(k,j)
          if (jj <= 0) call halt ('gap surface (ID=' // i_to_c(mesh%sset(n)%ID) // ') lies on boundary')
          nn = active%gside(kk,jj)
          if (nn == 0) then
            active%gside(kk,jj) = -n
            xside(n+1) = 1 + xside(n+1)
          else if (nn /= n) then
            call halt('overlapping side sets: ID=' // i_to_c(mesh%sset(n)%ID) // ',' // i_to_c(mesh%sset(nn)%ID))
          end if
        end if
      end do
    end do

    !! Setup XSIDE: SIDE(:,XSIDE(N):XSIDE(N+1)-1) are the missing sides for side set N.
    xside(1) = 1
    do j = 2, ubound(xside,1)
      xside(j) = xside(j) + xside(j-1)
    end do
    ASSERT( all(xside(:ubound(xside,1)-1) <= xside(2:)) )

    !! Fill-in the missing sides, correcting GSIDE as we go.  We use XSIDE as a
    !! temporary to hold the next free location in the list for each side set.
    allocate(side(2,xside(ubound(xside,1))-1))
    do j = 1, active%num_elem
      do k = 1, MAX_SIDE
        if (active%gside(k,j) < 0) then
          n = -active%gside(k,j)
          active%gside(k,j) = n
          side(1,xside(n)) = active%enum(j)
          side(2,xside(n)) = k
          xside(n) = 1 + xside(n)
        end if
      end do
    end do

    !! Restore XSIDE: the index of the first element of list J is now
    !! XSIDE(J-1) instead of XSIDE(J) as it should be -- fix this.
    do j = ubound(xside,1), 2, -1
      xside(j) = xside(j-1)
    end do
    xside(1) = 1

    ASSERT( all(active%gside >= 0) )

    !! Tag all inactive/nonexistent sides with a negative value in GSIDE.
    do j = 1, active%num_elem
      do k = 1, MAX_SIDE
        if (.not.btest(active%smask(j),k-1)) active%gside(k,j) = -1
      end do
    end do

  end subroutine generate_gside_data

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! DEFINE_EBLK
 !!
 !! This subroutine defines the element blocks in OUTMESH, adding blocks of
 !! gap elements to the element blocks from INMESH.  The bulk of the mesh
 !! transformation work is done here.
 !!
 !! The fundamental algorithm is simple:
 !!
 !! 1. Each active node of each active element is replaced by a new node.
 !!    This slices the mesh open along all active faces.
 !! 2. A 0-width gap element is inserted between neighboring elements across
 !!    a gap side.  This reconnects the mesh on gap sides.
 !! 3. On the remaining active sides, which are still sliced open, we identify
 !!    matching nodes from either side, and declare them to be equivalent.
 !!    This reconnects the neighbors, repairing the damage done in step 1.
 !!
 !! One virtue of this algorithm is that we don't need to determine a priori
 !! how many daughter nodes a gap node needs to be split into; this happens
 !! automatically during step 3.
 !!
 !! In practice, the algorithm is complicated by the element block structure
 !! of the Exodus mesh, the differing type of gap elements that can arise,
 !! and the need to place the gap elements that are created into specific
 !! element blocks, depending on element type and side set.
 !!
 !! In order to define the OUTMESH element blocks, this subroutine creates
 !! two or more coincident daughter nodes for each gap node.  The relationship
 !! between a gap node and its daughters is described by the returned DAUGHTER
 !! array.  When J <= INMESH%NUM_NODE, DAUGHTER(J)>0 identifies node J as a gap
 !! node.  Its first daughter replaces it as node J, and its next daughter is
 !! node DAUGHTER(J).  The complete set of daughter nodes is obtained by walking
 !! down the linked-list structure until a zero DAUGHTER value is encountered.
 !! All the daughters after the first are new nodes, and are numbered after
 !! the nodes in the input mesh.  DAUGHTER(J)=0 for all original inactive nodes;
 !! these nodes retain their number in the output mesh.
 !!

  subroutine define_eblk (inmesh, active, outmesh, daughter)

    type(exodus_mesh), intent(in)    :: inmesh
    type(active_mesh), intent(in)    :: active
    type(exodus_mesh), intent(inout) :: outmesh
    integer, pointer :: daughter(:)

    integer :: i, j, k, l, n, b, l1, l2, num_vert, num_elem, ID, ios, shift
    integer, allocatable :: parent(:), equiv(:), map(:)
    integer, allocatable :: gap_count(:,:)
    integer, pointer :: list(:), side1(:), side2(:)

    !! Count the gap elements by type and side set.
    allocate(gap_count(3:4,inmesh%num_sset))
    gap_count = 0
    do j = 1, active%num_elem
      list => inmesh%side_size_list(inmesh%eblk(active%bnum(j))%elem_type)
      do k = 1, MAX_SIDE
        n = active%gside(k,j)
        if (n > 0) gap_count(list(k),n) = 1 + gap_count(list(k),n)
      end do
    end do
    gap_count = gap_count/2 ! we double counted above

    !! Allocate the OUTMESH element blocks.
    outmesh%num_elem = inmesh%num_elem + sum(gap_count)
    outmesh%num_eblk = inmesh%num_eblk + count(gap_count > 0)
    allocate(outmesh%eblk(outmesh%num_eblk))

    !! Copy the INMESH element blocks to OUTMESH without change.
    do b = 1, inmesh%num_eblk
      outmesh%eblk(b)%ID = inmesh%eblk(b)%ID
      outmesh%eblk(b)%elem_type = inmesh%eblk(b)%elem_type
      num_vert = size(inmesh%eblk(b)%connect,dim=1)
      num_elem = size(inmesh%eblk(b)%connect,dim=2)
      allocate(outmesh%eblk(b)%connect(num_vert,num_elem))
      outmesh%eblk(b)%connect = inmesh%eblk(b)%connect
      outmesh%eblk(b)%num_elem = num_elem
    end do

    !! Prepare the new gap element blocks in OUTMESH.
    b = inmesh%num_eblk
    do n = 1, size(gap_count,dim=2)
      do l = lbound(gap_count,dim=1), ubound(gap_count,dim=1)
        if (gap_count(l,n) == 0) cycle

        b = b + 1

        !! Create the element block for this group of gap elements.
        select case (l)
        case (3)
          outmesh%eblk(b)%elem_type = 'WEDGE6'
          num_vert = 6
        case (4)
          outmesh%eblk(b)%elem_type = 'HEX8'
          num_vert = 8
        end select
        outmesh%eblk(b)%num_elem = gap_count(l,n)
        allocate(outmesh%eblk(b)%connect(num_vert,outmesh%eblk(b)%num_elem))
        outmesh%eblk(b)%connect = 0

        !! Overwrite the count with the block index for this gap element group.
        gap_count(l,n) = b

        !! Announce the element block.
        write(*,'(a)') 'Creating element block for ' // i_to_c(outmesh%eblk(b)%num_elem) &
          // ' ' // outmesh%eblk(b)%elem_type // ' gap elements from side set ' &
          // i_to_c(inmesh%sset(n)%ID)

        !! Obtain the element block ID from the user.
        outmesh%eblk(b)%ID = 0
        do while (outmesh%eblk(b)%ID == 0)
          write(*,'(4x,a)',advance='no') 'Enter element block ID: '
          read(*,*,iostat=ios) ID
          if (ios /= 0) then
            write(*,'(8x,a)') 'Error reading the ID.'
          else if (ID <= 0) then
            write(*,'(8x,a)') 'Sorry, the ID must be positive.'
          else if (any(ID == outmesh%eblk(:b)%ID)) then
            write(*,'(8x,a)') 'Sorry, that ID is already in use.'
          else
            outmesh%eblk(b)%ID = ID
            write(*,'(4x,a)') 'Setting element block ID to ' // i_to_c(ID)
          end if
        end do

        !! Write info to the log file
        write(LUN,'(a)') 'Creating element block ' // i_to_c(outmesh%eblk(b)%ID) // ' for ' &
          // i_to_c(outmesh%eblk(b)%num_elem) // ' ' // outmesh%eblk(b)%elem_type &
          // ' gap elements from side set ' // i_to_c(inmesh%sset(n)%ID)
      end do
    end do

    !!
    !! STEP 2: SLICE OPEN THE MESH
    !!
    !! Before we can insert gap elements, we need to slice open the output
    !! mesh by replacing each active node of each active element by a new
    !! daughter node.  This disconnects an element from its neighbors across
    !! active sides.  The daughter nodes are numbered consecutively following
    !! the nodes of the starting mesh.  PARENT(J) records the original node
    !! number for daughter node J.  For simplicity we extend the array to all
    !! nodes by setting PARENT(J) = J for original nodes.
    !!

    !! Count the number of daughter nodes that will be created.
    n = 0
    do j = 1, active%num_elem
      n = n + bit_count(active%nmask(j))
    end do

    allocate(parent(n+inmesh%num_node))
    do j = 1, inmesh%num_node
      parent(j) = j
    end do
    parent(1+inmesh%num_node:) = 0

    !! Replace each active node of each active element with a new daughter node.
    n = inmesh%num_node
    do j = 1, active%num_elem
      associate(list => outmesh%eblk(active%bnum(j))%connect(:,active%lnum(j)))
        do k = 1, size(list)
          if (btest(active%nmask(j),k-1)) then
            n = n + 1
            parent(n) = list(k)
            list(k) = n
          end if
        end do
      end associate
    end do

    ASSERT( n == ubound(parent,1) )
    ASSERT( minval(parent) > 0 .and. maxval(parent) <= inmesh%num_node )

    !!
    !! STEP 3: INSERT GAP ELEMENTS AND CLOSE UP THE MESH
    !!
    !! The preceding step sliced open the output mesh along all active sides.
    !! On gap sides we now need to create a gap element that links between the
    !! two elements that were disconnected from each other.  On the remaining
    !! active sides we need to reconnect the two elements by equivalencing
    !! matching pairs of nodes on the element sides.
    !!
    !! The EQUIV array is a linked-list data structure describing the state
    !! of the evolving node equivalence class.  If EQUIV(J) equals 0, node J
    !! is the representative for its class.  Otherwise node J is equivalent
    !! to node EQUIV(J), and the representative for its class can be found
    !! by walking down the linked list to the end.
    !!

    allocate(equiv(ubound(parent,1)))
    equiv = 0

    outmesh%eblk(inmesh%num_eblk+1:)%num_elem = 0   ! use these as index temporaries
    do j = 1, active%num_elem
      do k = 1, MAX_SIDE

        if (active%gside(k,j) < 0) cycle ! not an active side
        if (j > active%jnbr(k,j))  cycle ! handle this one from the other side

        !! The co-oriented side node lists for the adjacent elements.
        side1 => outmesh%side_node_list(active%enum(j), k)
        side2 => outmesh%side_node_list(active%enum(active%jnbr(k,j)), active%knbr(k,j), reverse=.true.)
        ASSERT( associated(side1) .and. associated(side2) )
        ASSERT( size(side1) == size(side2) )

        !! Rotate the SIDE2 list so that it is in 1-1 correspondence with SIDE1.
        shift = 0
        do while (parent(side2(1+shift)) /= parent(side1(1)))
          shift = 1 + shift
          ASSERT( shift < size(side1) )
        end do
        side2 = cshift(side2, shift)
        ASSERT( all(parent(side1) == parent(side2)) )

        select case (active%gside(k,j))
        case (1:)

          !! This is a gap side.  Create a gap element connecting the adjacent
          !! elements for this side; this works for HEX8 and WEDGE6 elements.

          b = gap_count(size(side1),active%gside(k,j)) ! destination element block
          n = 1 + outmesh%eblk(b)%num_elem
          outmesh%eblk(b)%connect(:size(side1),  n) = side1
          outmesh%eblk(b)%connect(1+size(side1):,n) = side2
          outmesh%eblk(b)%num_elem = n

        case (0) ! open side

          !! This is an open side.  Reconnect the adjacent elements by making
          !! the matching nodes of the two sides equivalent.

          do i = 1, size(side1)
            l1 = rep_node(side1(i)) ! representative of SIDE1(I)'s class
            l2 = rep_node(side2(i)) ! representative of SIDE2(I)'s class
            if (l1 /= l2) equiv(l1) = l2
          end do

        end select

        deallocate(side1, side2)

      end do
    end do

    ASSERT( all(outmesh%eblk%defined()) )
    ASSERT( all(equiv(:inmesh%num_node) == 0) )
    ASSERT( minval(equiv,equiv/=0) > inmesh%num_node )
    ASSERT( maxval(equiv,equiv/=0) <= ubound(equiv,1) )

    !! Count the number of daughter nodes arising from each active node,
    !! using the unused initial part of the EQUIV array as a temporary.
    do j = 1+inmesh%num_node, ubound(parent,1)
      if (equiv(j) == 0) equiv(parent(j)) = 1 + equiv(parent(j))
    end do

    !! Warn if any gap nodes were replaced by too many daughters.
    call daughter_check (equiv(:inmesh%num_node))

    !! Reuse each (now defunct) gap node by making it the representative
    !! for one of its daughter nodes; restore the EQUIV array as we go.
    do j = 1+inmesh%num_node, ubound(parent,1)
      if (equiv(j) /= 0) cycle
      if (equiv(parent(j)) > 0) then
        equiv(j) = parent(j)
        equiv(parent(j)) = 0
      end if
    end do

    ASSERT( all(equiv(:inmesh%num_node) == 0) )

    !!
    !! STEP 4: Number the nodes
    !!
    !! We are now ready to generate the (consecutive) node numbering for the
    !! output mesh.  All the original nodes retain their number, including the
    !! gap nodes which were replaced by one of their daughters.  The remaining
    !! new daughter nodes are numbered consecutively following the others.
    !!

    outmesh%num_node = count(equiv == 0)
    allocate(map(ubound(equiv,1)))
    do j = 1, inmesh%num_node
      map(j) = j
    end do

    !! Consecutive numbering of representative nodes.
    n = inmesh%num_node
    do j = 1+inmesh%num_node, ubound(map,1)
      if (equiv(j) /= 0) cycle
      n = n + 1
      map(j) = n
    end do
    ASSERT( n == outmesh%num_node )

    !! Number non-representative nodes as their representative.
    do j = 1+inmesh%num_node, ubound(map,1)
      if (equiv(j) /= 0) map(j) = map(rep_node(j))
    end do

    !! Renumber the active nodes of the original active elements.
    do j = 1, active%num_elem
      associate(list => outmesh%eblk(active%bnum(j))%connect(:,active%lnum(j)))
        do k = 1, size(list)
          if (btest(active%nmask(j),k-1)) list(k) = map(list(k))
        end do
      end associate
    end do

    !! Renumber the nodes of the new gap elements.
    do b = inmesh%num_eblk+1, outmesh%num_eblk
      n = 0
      do l = 1, outmesh%eblk(b)%num_elem
        associate(list => outmesh%eblk(b)%connect(:,l))
          list = map(list)
          if (.not.unique(list)) n = n + 1  ! degenerate gap element
        end associate
      end do
      if (n > 0) call warn ('element block ID ' // i_to_c(outmesh%eblk(b)%ID) // &
                            ' contains ' // i_to_c(n) // ' degenerate gap elements')
    end do

    !! Generate the gap node daughter structure.
    allocate(daughter(outmesh%num_node))
    daughter = 0
    do j = ubound(parent,1), 1+inmesh%num_node, -1
      if (equiv(j) /= 0) cycle
      daughter(map(j)) = daughter(parent(j))
      daughter(parent(j)) = map(j)
    end do
    ASSERT( minval(daughter,daughter/=0) > inmesh%num_node )
    ASSERT( maxval(daughter) <= size(daughter) )

    deallocate(equiv, parent, map)

  contains

    !!
    !! This auxillary function returns the representative node for the
    !! argument node's equivalence class.
    !!

    integer function rep_node (n) result (r)
      integer, intent(in) :: n
      r = n
      do while (equiv(r) /= 0)
        r = equiv(r)
      end do
    end function rep_node

    !!
    !! This auxillary subroutine looks for gap nodes that have been replaced
    !! by more than 3 daughter nodes and prints out an appropriate warning.
    !!

    subroutine daughter_check (dcount)

      integer, intent(in) :: dcount(:)

      integer :: j, n

      n = count(dcount > 3)
      if (n == 0) return

      do j = 1, size(dcount)
        if (dcount(j) > 3) exit
      end do

      call warn ('node ' // i_to_c(j) // ' was replaced by ' // i_to_c(dcount(j)) // ' daughter nodes.')
      call warn ( i_to_c(n) // ' nodes were replace by >3 daughter nodes')
      call warn ('output mesh may have a gap element topology that Truchas cannot handle.')

    end subroutine daughter_check

   !!
   !! This auxillary function returns the value true if the integer values in
   !! the vector argument are unique; otherwise it returns the value false.
   !!

    pure logical function unique (list)

      integer, intent(in) :: list(:)

      integer :: i, j

      unique = .false.
      do i = 1, size(list)-1
        do j = i+1, size(list)
          if (list(i) == list(j)) return
        end do
      end do
      unique = .true.

    end function unique

  end subroutine define_eblk

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! DEFINE_COORD
 !!
 !! This subroutine allocates and defines the COORD array component of OUTMESH.
 !! The output mesh contains all the nodes of the input mesh plus additional
 !! nodes that were created when an original gap node was split into 2 or more
 !! coincident daughter nodes.  The relationship between a gap node and its
 !! daughters is described by the DAUGHTER array; see the comments for the
 !! subroutine DEFINE_EBLK for a discussion of this array.
 !!

  subroutine define_coord (inmesh, daughter, outmesh)

    type(exodus_mesh), intent(in)    :: inmesh
    integer,           intent(in)    :: daughter(:)
    type(exodus_mesh), intent(inout) :: outmesh

    integer :: j, d
    integer :: map(outmesh%num_node)

    !! Generate the new-to-old node mapping array (onto, many-to-one).
    map = 0
    do j = 1, inmesh%num_node
      d = daughter(j)
      do while (d /= 0)
        ASSERT( map(d) == 0 )
        map(d) = j
        d = daughter(d)
      end do
      map(j) = j
    end do

    ASSERT( minval(map) > 0 )
    ASSERT( maxval(map) <= inmesh%num_node )

    !! Define COORD; daughter nodes inherit the coordinates of their parent.
    outmesh%num_dim = inmesh%num_dim
    allocate(outmesh%coord(outmesh%num_dim,outmesh%num_node))
    do j = 1, outmesh%num_node
      outmesh%coord(:,j) = inmesh%coord(:,map(j))
    end do


  end subroutine define_coord

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! DEFINE_NSET
 !!
 !! This subroutine allocates and defines the NSET structure-array component of
 !! OUTMESH, copying the data from the input mesh.  Each gap node in a node set
 !! is replaced by its set of daughter nodes.  The relationship between a gap
 !! node and its daughters is described by the DAUGHTER array; see the comments
 !! for the subroutine DEFINE_EBLK for a discussion of this array.
 !!

  subroutine define_nset (inmesh, daughter, outmesh)

    type(exodus_mesh), intent(in)    :: inmesh
    integer,           intent(in)    :: daughter(:)
    type(exodus_mesh), intent(inout) :: outmesh

    integer :: i, j, n, d, extra

    outmesh%num_nset = inmesh%num_nset
    if (outmesh%num_nset == 0) return
    allocate(outmesh%nset(outmesh%num_nset))

    do n = 1, inmesh%num_nset

      !! Count the number of extra nodes in the transformed node set.
      extra = 0
      do j = 1, inmesh%nset(n)%num_node
        d = daughter(inmesh%nset(n)%node(j))
        do while (d /= 0)
          extra = 1 + extra
          d = daughter(d)
        end do
      end do

      !! Copy the node set, replacing each gap node by its daughters.
      outmesh%nset(n)%ID = inmesh%nset(n)%ID
      outmesh%nset(n)%num_node = inmesh%nset(n)%num_node + extra
      allocate(outmesh%nset(n)%node(outmesh%nset(n)%num_node))
      i = 0
      do j = 1, inmesh%nset(n)%num_node
        i = i + 1
        outmesh%nset(n)%node(i) = inmesh%nset(n)%node(j)
        d = daughter(inmesh%nset(n)%node(j))
        do while (d /= 0)
          i = i + 1
          outmesh%nset(n)%node(i) = d
          d = daughter(d)
        end do
      end do

    end do

    ASSERT( all(outmesh%nset%defined()) )

  end subroutine define_nset

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! DEFINE_SSET
 !!
 !! This subroutine allocates and defines the SSET structure-array component
 !! of OUTMESH, copying the data from the input mesh.  The XSIDE/SIDE arrays
 !! describe additional sides that should be included in the side sets of the
 !! output mesh: SIDE(I,XSIDE(N):XSIDE(N+1)-1) are the elements (I=1) and
 !! their corresponding local sides (I=2) that should be added to side set N.
 !! If STRICT_SS_TRANSF has the value true, then these additional sides are
 !! ignored, and the side sets are copied from the input mesh unaltered.
 !!

  subroutine define_sset (inmesh, xside, side, outmesh)

    type(exodus_mesh), intent(in)    :: inmesh
    integer,           intent(in)    :: xside(:), side(:,:)
    type(exodus_mesh), intent(inout) :: outmesh

    integer :: n, extra

    outmesh%num_sset = inmesh%num_sset
    allocate(outmesh%sset(outmesh%num_sset))

    do n = 1, outmesh%num_sset

      extra = xside(n+1) - xside(n)
      if (extra > 0 .and. strict_ss_transf) then
        extra = 0
        call warn ('side set ' // i_to_c(inmesh%sset(n)%ID) // &
                   ' does not fully describe both sides of the gap surface')
      end if

      outmesh%sset(n)%ID = inmesh%sset(n)%ID
      outmesh%sset(n)%num_side = inmesh%sset(n)%num_side + extra
      allocate(outmesh%sset(n)%elem(outmesh%sset(n)%num_side))
      allocate(outmesh%sset(n)%face(outmesh%sset(n)%num_side))
      outmesh%sset(n)%elem(:inmesh%sset(n)%num_side) = inmesh%sset(n)%elem
      outmesh%sset(n)%face(:inmesh%sset(n)%num_side) = inmesh%sset(n)%face

      if (extra > 0) then
        outmesh%sset(n)%elem(1+inmesh%sset(n)%num_side:) = side(1,xside(n):xside(n+1)-1)
        outmesh%sset(n)%face(1+inmesh%sset(n)%num_side:) = side(2,xside(n):xside(n+1)-1)
      end if

    end do

    ASSERT( all(outmesh%sset%defined()) )

  end subroutine define_sset

  subroutine warn (mesg)
    character(len=*), intent(in) :: mesg
    write(*,'(a)') 'WARNING: ' // mesg
    write(LUN,'(a)') 'WARNING: ' // mesg
  end subroutine warn

  subroutine halt (mesg)
    character(len=*), intent(in) :: mesg
    write(0,'(a)') 'ERROR: ' // mesg
    stop
  end subroutine halt

  !!
  !!  This auxillary function returns the number of 1 bits in its argument.
  !!

  pure integer function bit_count (i)
    integer, intent(in) :: i
    integer :: n
    n = i
    bit_count = 0
    do while (n /= 0)
      if (btest(n,0)) bit_count = 1 + bit_count
      n = ishft(n,-1)
    end do
  end function bit_count

end module addgaps_proc
