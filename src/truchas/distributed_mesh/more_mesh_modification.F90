!!
!! MESH_MODIFICATION
!!
!! This module provides several routines that apply specialized schemes that
!! modify a primitive mesh (in its intermediate serial EXTERNAL_MESH-type
!! form) prior to generating the final distributed mesh.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Modified for general Exodus meshes, June 2015.
!!
!! PROGRAMMING INTERFACE
!!
!! The following routines operate on EXTERNAL_MESH-type meshes
!!
!!  CALL CONVERT_CELLS_TO_LINKS (MESH, BLOCK_ID, STAT, ERRMSG) removes the
!!    specified set of cells from the mesh and converts them into a special
!!    link data structure from which the original mesh connectivity may be
!!    recovered.  The cells must form a 1-cell thick layer, which when removed
!!    create gap interface in the mesh.  Each of these cells must have exactly
!!    one pair of opposite faces that match faces in the remaining mesh.
!!    MESH is the EXTERNAL_MESH-type mesh and BLOCK_ID is a rank-1 integer
!!    array of cell block IDs identifying the blocks of cells that are to be
!!    converted.  The modified mesh is returned in MESH.  STAT returns a
!!    nonzero integer value if and error is encountered and an explanatory
!!    error message is assigned to the string ERRMSG.  This is a serial
!!    procedure and is only applicable to hexahedral meshes (non-degenerate).
!!
!! CALL CREATE_INTERNAL_INTERFACES (MESH, SSID, STAT, ERRMSG) slices the
!!    input mesh open along the mesh surfaces specified by the list SSID
!!    of side set IDs.  This involves duplicating nodes along the surface
!!    and disconnecting neighboring cells across the surface to create an
!!    internal interface.  This interface appears as a boundary of the mesh,
!!    but a special link data structure is added the MESH object that encodes
!!    the original cell connectivity across the interface.  The specified
!!    surface must not lie on the boundary of the input mesh.  STAT returns
!!    a nonzero value if an error is encountered and an explanatory error
!!    message is assigned to the string ERRMSG.  This is a serial procedure
!!    and is applicable to both tet and hex meshes (non-degenerate).
!!

#include "f90_assert.fpp"

module more_mesh_modification

  use kinds, only: r8
  use exodus_mesh_type
  use string_utilities, only: i_to_c
  implicit none
  private

  public :: convert_cells_to_links
  public :: create_internal_interfaces
  
  type, public :: link_mesh
    integer :: nlink=0, nlblock=0
    integer, allocatable :: xlnode(:), lnode(:)
    integer, allocatable :: link_block(:)
    integer, allocatable :: link_block_id(:)
  end type link_mesh

contains

  subroutine convert_cells_to_links (mesh, block_id, stat, errmsg)

    use mesh_importer, only: external_mesh, side_set
    use cell_topology
    use facet_hash_type
    use string_utilities, only: i_to_c

    type(external_mesh), intent(inout), target :: mesh
    integer, intent(in) :: block_id(:)
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    integer :: i, j, k, n, side_bitmask, node_bitmask
    logical, allocatable :: link_cell(:), link_node(:), keep(:)
    integer, allocatable :: link_side(:), xbin(:), map(:), old_cnode(:,:), old_cell_block(:), tmp(:)
    integer, pointer :: cnode(:), fnode(:)
    type(side_set), pointer :: old_sset(:)
    type(facet_hash) :: hpar

    type :: table_entry
      integer, pointer :: fnode(:) => null()
    end type table_entry
    type(table_entry), pointer :: table(:), bin(:)

    !! Hex cell parameters
    integer, parameter :: NUM_FACE = 6
    integer, parameter :: FACE_SIZE = 4
    integer, parameter :: LINK_12 = int(b'000011')
    integer, parameter :: LINK_34 = int(b'001100')
    integer, parameter :: LINK_56 = int(b'110000')
    integer, parameter :: ROTATE_12(8) = (/ 1, 5, 6, 2, 4, 8, 7, 3 /)
    integer, parameter :: ROTATE_34(8) = (/ 1, 4, 8, 5, 2, 3, 7, 6 /)

    stat = 0
    errmsg = ''

    if (size(block_id) == 0) return ! nothing to do

    !! This subroutine can only be applied to hex meshes.
    if (mesh%mesh_type /= 'HEX') then
      stat = -1
      errmsg = 'CONVERT_CELLS_TO_LINKS: invalid mesh type: ' // trim(mesh%mesh_type)
      return
    end if

    !! Check that the specifed block IDs are valid.
    do j = 1, size(block_id)
      do i = size(mesh%block_id), 1, -1
        if (block_id(j) == mesh%block_id(i)) exit
      end do
      if (i == 0) then
        stat = -1
        errmsg = 'CONVERT_CELLS_TO_LINKS: unknown element block ID: ' // i_to_c(block_id(j))
        return
      end if
    end do

    !! Tag the cells that are to be converted into link cells.
    allocate(link_cell(mesh%ncell))
    do j = 1, mesh%ncell
      link_cell(j) = any(mesh%cell_block(j) == block_id)
    end do

    !! Tag the nodes that lie along the gap: all link cell nodes.
    allocate(link_node(mesh%nnode))
    link_node = .false.
    do j = 1, mesh%ncell
      if (link_cell(j)) link_node(mesh%cnode(:,j)) = .true.
    end do

    !! Identify those remaining cell sides that may match a link cell side;
    !! these are sides whose nodes are all link nodes.  This info is stored
    !! in the bitmask array LINK_SIDE.
    n = 0 ! side count
    allocate(link_side(mesh%ncell))
    do j = 1, mesh%ncell
      link_side(j) = 0
      if (link_cell(j)) cycle  ! not interested in these cells
      !! Tag the link nodes on this cell.
      cnode => mesh%cnode(:,j)
      node_bitmask = 0
      do k = 1, size(cnode)
        if (link_node(cnode(k))) node_bitmask = ibset(node_bitmask, k-1)
      end do
      if (node_bitmask == 0) cycle ! no link nodes on this cell at all
      !! Tag the cell sides whose nodes are all link nodes.
      side_bitmask = 0
      do k = 1, size(HEX8_FACE_SIG)
        if (iand(node_bitmask, HEX8_FACE_SIG(k)) == HEX8_FACE_SIG(k)) then
          n = n + 1
          side_bitmask = ibset(side_bitmask, k-1)
        end if
      end do
      link_side(j) = side_bitmask
    end do
    deallocate(link_node)

   !!
   !! Create a bin table of all the cell sides identified above.  We use a
   !! hash of the side nodes to generate the bin number for each side.  This
   !! data structure enables efficient searching of the sides.
   !!

    allocate(table(n))
    call hpar%init (n, mesh%nnode)  ! adjusts N to a power of 2
    allocate(xbin(0:n))

    !! Count the number of hits to each bin; store the count for bin N in XBIN(N+1).
    xbin = 0
    do j = 1, mesh%ncell
      if (link_side(j) == 0) cycle
      do k = 1, NUM_FACE
        if (btest(link_side(j), k-1)) then
          fnode => hex_face_nodes (mesh%cnode(:,j), k, normalize=.true.)
          call hpar%hash (fnode, n)
          xbin(n+1) = xbin(n+1) + 1
          deallocate(fnode)
        end if
      end do
    end do

    !! Prepare XBIN: bin J will be TABLE(XBIN(J):XBIN(J+1)-1).
    xbin(0) = 1
    do j = 1, ubound(xbin,1)
      xbin(j) = xbin(j) + xbin(j-1)
    end do

    !! Fill the table; use XBIN as a temporary to hold the next free location for each bin.
    do j = 1, mesh%ncell
      do k = 1, NUM_FACE
        if (btest(link_side(j), k-1)) then
          fnode => hex_face_nodes (mesh%cnode(:,j), k, normalize=.true.)
          call hpar%hash (fnode, n)
          i = xbin(n)
          table(i)%fnode => fnode
          xbin(n) = i + 1
        end if
      end do
    end do
    deallocate(link_side)

    !! Restore XBIN: the index of the first element of bin J is now XBIN(J-1)
    !! instead of XBIN(J) as it should be -- fix this.
    do j = ubound(xbin,1), 1, -1
      xbin(j) = xbin(j-1)
    end do
    xbin(0) = 1

   !!
   !! Verify that the specified cells are valid link cells: exactly one pair
   !! of opposing sides must match sides of the remaining cells.  The table
   !! created above contains all the possible candidates.  We also reorient
   !! the link cells, if necessary, so that it is sides 5/6 (bottom/top) that
   !! lie along the gap surfaces.
   !!

    do j = 1, mesh%ncell
      if (link_cell(j)) then
        cnode => mesh%cnode(:,j)
        !! For each side of the link cell search the table for a matching side.
        side_bitmask = 0  ! sides that match
        do k = 1, NUM_FACE
          !! The FNODE the neighbor cell would have for this side.
          fnode => hex_face_nodes (cnode, k, reverse=.true., normalize=.true.)
          call hpar%hash (fnode, n)
          bin => table(xbin(n):xbin(n+1)-1)
          do i = size(bin), 1, -1
            if (all(fnode == bin(i)%fnode)) exit
          end do
          if (i /= 0) side_bitmask = ibset(side_bitmask, k-1)
          deallocate(fnode)
        end do
        !! Check the signature of the matching sides and reorient if necessary.
        select case (side_bitmask)
        case (LINK_12)
          cnode = cnode(ROTATE_12)
        case (LINK_34)
          cnode = cnode(ROTATE_34)
        case (LINK_56)
          ! no reorientation necessary
        case default
          stat = -1
          errmsg = 'CONVERT_CELLS_TO_LINKS: unrecognized link cell ' // i_to_c(j)
          return
        end select
      end if
    end do

    !! Delete the table -- we're finished with it.
    do j = 1, size(table)
      if (associated(table(i)%fnode)) deallocate(table(i)%fnode)
    end do
    deallocate(table, xbin)

   !!
   !! Convert the specified cells into link cells.  The cells that remain must
   !! be repacked into contiguous arrays along with their associated cell block
   !! info.  We also need to update the list of block IDs and the side set data.
   !!

    !! New link cell data array.
    INSIST(mesh%nlink == 0)
    deallocate(mesh%lnode, mesh%link_block, mesh%link_block_id)
    mesh%nlink = count(link_cell)
    allocate(mesh%lnode(2*FACE_SIZE,mesh%nlink), mesh%link_block(mesh%nlink))

    !! Resized mesh cell arrays.
    mesh%ncell = mesh%ncell - mesh%nlink
    call move_alloc (mesh%cnode, old_cnode)
    call move_alloc (mesh%cell_block, old_cell_block)
    allocate(mesh%cnode(8,mesh%ncell), mesh%cell_block(mesh%ncell))

    !! Copy the old cell data into the new data arrays.  If cell J is not a link
    !! cell, then MAP(J) is its new cell number; MAP(J) = 0 for link cells.
    n = 0
    i = 0
    allocate(map(size(old_cnode,dim=2)))
    do j = 1, size(old_cnode,dim=2)
      if (link_cell(j)) then
        n = n + 1
        map(j) = 0  ! dummy value -- shouldn't be used
        mesh%lnode(:,n) = old_cnode(:,j)
        mesh%link_block(n) = old_cell_block(j)
      else
        i = i + 1
        map(j) = i
        mesh%cnode(:,i) = old_cnode(:,j)
        mesh%cell_block(i) = old_cell_block(j)
      end if
    end do
    deallocate(old_cnode, old_cell_block)
    ASSERT(n == mesh%nlink)
    ASSERT(i == mesh%ncell)

    !! Fix-up the block ID list; the link cell blocks are gone now.
    allocate(keep(mesh%nblock))
    do i = 1, mesh%nblock
      keep(i) = all(mesh%block_id(i) /= block_id)
    end do
    n = count(keep)
    allocate(mesh%link_block_id(mesh%nblock-n))
    mesh%nlblock = size(mesh%link_block_id)
    mesh%link_block_id = pack(mesh%block_id, mask=.not.keep)
    if (n < mesh%nblock) then
      call move_alloc (mesh%block_id, tmp)
      allocate(mesh%block_id(n))
      mesh%nblock = n
      mesh%block_id = pack(tmp, mask=keep)
      deallocate(tmp)
    end if
    deallocate(keep)

    !! Fix-up the sideset data: some elements are gone and the rest renumbered.
    do i = 1, size(mesh%sset)
      allocate(keep(mesh%sset(i)%num_side))
      keep = .not.link_cell(mesh%sset(i)%elem)
      n = count(keep)
      mesh%sset(i)%num_side = n
      call move_alloc (mesh%sset(i)%elem, tmp)
      allocate(mesh%sset(i)%elem(n))
      mesh%sset(i)%elem = map(pack(tmp, mask=keep))
      deallocate(tmp)
      call move_alloc (mesh%sset(i)%face, tmp)
      allocate(mesh%sset(i)%face(n))
      mesh%sset(i)%face = pack(tmp, mask=keep)
      deallocate(tmp, keep)
    end do
    deallocate(link_cell)

    !! Prune any empty side sets that may have resulted.
    n = count(mesh%sset%num_side > 0)
    if (n < size(mesh%sset)) then
      old_sset => mesh%sset
      allocate(mesh%sset(n))
      n = 0
      do i = 1, size(old_sset)
        if (old_sset(i)%num_side > 0) then
          n = n + 1
          mesh%sset(n) = old_sset(i)
        else  ! deallocate the 0-sized array pointer components
          deallocate(old_sset(i)%elem, old_sset(i)%face)
        end if
      end do
      ASSERT(n == size(mesh%sset))
      deallocate(old_sset)
    end if

  end subroutine convert_cells_to_links

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! CREATE_INTERNAL_INTERFACES
 !!
 !! This subroutine modifies a mesh by creating internal interfaces within
 !! the mesh.  It does this by slicing the mesh open along a mesh-conforming
 !! surface: nodes on the surface are duplicated and neighboring cells across
 !! the surface disconnected from one another by redefining their node
 !! connectivity.  The internal interface that this creates constitutes part
 !! of the actual boundary of the modified mesh, but has the special property
 !! that it is composed of matching pairs of boundary faces that arose from a
 !! a surface face originally shared by neighboring cells on either side of
 !! the surface.  The relationship between these pairs of faces is described
 !! by a special link data structure component of the mesh object.  In this
 !! case a link can be regarded as a special 0-thickness cell bridging the
 !! pair of faces:
 !!                                                              6
 !!                        8-----7                              /:\ 
 !!                       /:    /:                 +           / : \ 
 !!       +-----+        5-----6 :                / \         4-----5
 !!      /     /   ==>   : :   : :      and      /   \   ==>  :  :  :
 !!     +-----+          : 4---:-3              +-----+       :  3  :
 !!                      :/    :/                             : / \ :
 !!                      1-----2                              :/   \:
 !!                                                           1-----2
 !!
 !! The surface along which the interface is to be created is specified by
 !! giving the IDs of the side sets that describe the surface in the rank-1
 !! integer array SSID.  This surface must not lie on any part of the existing
 !! mesh boundary.  Moreover each face on the surface must be associated with
 !! a unique side set from the specified list, and this ID becomes the link
 !! block ID of the link that is created for it.
 !!
 !! Besides creating (or adding to an existing) link data structure, the
 !! mesh is modified in a number of ways.  Duplicate nodes along the
 !! interface are created coincident with the original.  All other nodes
 !! retain their original identity.  No new cells are created and they
 !! all retain their original identity, except that those bordering the
 !! interface will be connected to different (but geometrically identical)
 !! nodes.
 !!
 !! IMPLEMENTATION NOTES
 !!
 !! The basic algorithm is relatively simple.  For each cell we replace each
 !! surface node by a new coincident node.  This slices the mesh open along
 !! the surface as desired, but also slices the mesh open on any face that
 !! contains a surface node.  We then repair this collateral damage by
 !! equivalencing matching pairs of nodes on these non-surface faces.  For the
 !! remaining surface faces we create the links between matching face pairs.
 !! One virtue of this algorithm is that we do not need to determine a priori
 !! how many daughter nodes a surface node needs to be split into.  This
 !! happens automatically during the equivalencing step, and as a result the
 !! algorithm can handle arbitrarily complex intersecting surfaces.
 !!
 !! In order to carry out the algorithm we need cell neighbor data: for each
 !! side of a cell what is the neighboring cell across that side and what is
 !! the corresponding side of the neighbor.  Computing and storing this data
 !! mesh-wide is expensive and generally unnecessary; the portion of the
 !! mesh involved in the algorithm is typically a very small fraction of the
 !! whole mesh, especially for large meshes.  So the first step in this
 !! implementation is to identify the active subset of the mesh.  An active
 !! node is simply a node on the surface.  An active cell is one containing
 !! an active node, and an active side in one containing an active node.
 !! Interface sides are those active sides lying on the surface.  The
 !! collection of all the active cells forms the active submesh; all other
 !! cells and all inactive sides and nodes can be ignored as the are not
 !! involved in the algorithm.  The neighbor data is then computed and stored
 !! for the active submesh only.
 !!

  subroutine create_internal_interfaces (ssid, mesh, link, stat, errmsg)

    use facet_hash_type
    use cell_topology
    use integer_set_type

    integer, intent(in) :: ssid(:)
    type(exodus_mesh), intent(inout), target :: mesh
    type(link_mesh), intent(inout), target :: link
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: i, j, k, n, n1, n2, offset, l1, l2, bit_mask, shift, jj, kk, nn, new_nnode
    integer, pointer :: list(:), side_sig(:), link_side(:)
    integer, allocatable :: side1(:), side2(:), sides(:), cell_side(:), itmp(:)
    logical, allocatable :: active_sset(:), active_node(:), active_cell(:)
    integer, allocatable :: parent(:), equiv(:), map(:), daughter(:)
    real(r8), allocatable :: new_x(:,:)
    type(integer_set), allocatable :: extra_sides(:)

    !! Data arrays for the active submesh.
    integer, allocatable :: link_jnbr(:,:), link_knbr(:,:)

    !! Local variables used to compute the neighbor data.
    type :: table_entry
      integer :: j, k
      integer, allocatable :: side(:)
    end type table_entry
    type(table_entry), pointer :: table(:)
    integer, allocatable :: xbin(:)
    type(facet_hash) :: hpar
    
    !! Data structure to hold working data about active cells
    type :: active_elem
      integer :: cnum, nmask, smask, nface
      integer, pointer :: connect(:) => null()
      integer, allocatable :: jnbr(:), knbr(:), iside(:)
      integer, pointer :: side_size(:) => null()
    end type
    type(active_elem), allocatable :: acell(:)

    stat = 0

    !! Put a 0-sized input link mesh into a known state. 
    if (.not.allocated(link%xlnode)) then
      link%nlink = 0
      link%nlblock = 0
      allocate(link%xlnode(1), link%lnode(0), link%link_block_id(0), link%link_block(0))
      link%xlnode(1) = 1
    end if

    if (size(ssid) == 0) return ! nothing to do

    !! Generate the active side set mask array ACTIVE_SSET.
    allocate(active_sset(size(mesh%sset)))
    active_sset = .false.
    do i = 1, size(ssid)
      do n = size(mesh%sset), 1, -1
        if (mesh%sset(n)%ID == ssid(i)) exit
      end do
      if (n == 0) then
        stat = -1
        errmsg = 'CREATE_INTERNAL_INTERFACES: unknown side set ID: ' // i_to_c(ssid(i))
        return
      end if
      active_sset(n) = .true.
    end do

    !! Generate the active node mask array ACTIVE_NODE.
    allocate(active_node(mesh%num_node))
    active_node = .false.
    do n = 1, size(mesh%sset)
      if (active_sset(n)) then
        list => mesh%side_set_node_list(n)
        do i = 1, size(list)
          active_node(list(i)) = .true.
        end do
        deallocate(list)
      end if
    end do

    !! Generate the active cell mask array ACTIVE_CELL.
    allocate(active_cell(mesh%num_elem))
    n = 0
    do i = 1, mesh%num_eblk
      do j = 1, mesh%eblk(i)%num_elem
        n = n + 1
        active_cell(n) = any(active_node(mesh%eblk(i)%connect(:,j)))
      end do
    end do
    ASSERT(n == size(active_cell))
    
    !! Generate the structure array ACELL of active cell data.  The components
    !! of an element are:
    !!  - the cell index CNUM
    !!  - the connectivity data for the cell; this is a pointer to a section of
    !!    MESH%EBLK%CONNECT array data
    !!  - the size of each side of the cell: SIDE_SIZE(K) is the number of nodes
    !!    on side K of the cell; this is a pointer to static data associated
    !!    with the type of the cell (hex, tet, etc.)
    !!  - the active node bit mask NMASK:
    !!    BTEST(NMASK,K-1) is true when node K of the cell is an active node
    !!  - the active side bit mask SMASK:
    !!    BTEST(SMASK,K-1) is true when side K of the cell is an active side
    !!  The side-based component arrays JNBR, KNBR, and ISIDE are also allocated
    !!  but only initialized with default values.

    allocate(acell(count(active_cell)))
    n = 0; offset = 0
    do i = 1, mesh%num_eblk
      do j = 1, mesh%eblk(i)%num_elem
        if (active_cell(j+offset)) then
          n = n + 1
          acell(n)%cnum = j + offset
          acell(n)%connect => mesh%eblk(i)%connect(:,j)
          acell(n)%side_size => cell_face_sizes(acell(n)%connect)
          acell(n)%nface = size(acell(n)%side_size)
          !! Define NMASK: mark the active nodes on this active cell.
          bit_mask = 0
          do k = 1, size(acell(n)%connect)
            if (active_node(acell(n)%connect(k))) bit_mask = ibset(bit_mask, k-1)
          end do
          acell(n)%nmask = bit_mask
          !! Define SMASK: mark the active sides on this active cell.
          side_sig => cell_face_sig(acell(n)%connect)
          bit_mask = 0
          do k = 1, size(side_sig)
            if (iand(acell(n)%nmask, side_sig(k)) /= 0) bit_mask = ibset(bit_mask, k-1)
          end do
          acell(n)%smask = bit_mask
          !! Allocate space for neighbor info
          allocate(acell(n)%jnbr(acell(n)%nface), acell(n)%knbr(acell(n)%nface), acell(n)%iside(acell(n)%nface))
          acell(n)%jnbr = -1
          acell(n)%knbr = -1
          acell(n)%iside = 0
        end if
      end do
      offset = offset + mesh%eblk(i)%num_elem
    end do
    ASSERT(n == size(acell))

    deallocate(active_cell)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Generate the neighbor component arrays ACELL%JNBR and ACELL%KNBR.  When
   !! positive, ACELL(J)%JNBR(K) is the active cell that neighbors active cell
   !! J across side K, and ACELL(J)%KNBR(K) is the corresponding side of the
   !! neighbor.  Boundary sides are indicated by a zero ACELL%JNBR value,
   !! and inactive sides to be ignored are indicated by a negative value.
   !! These arrays satisfy the reciprocity relations
   !!
   !!      JJ = ACELL(J)%JNBR(K)>0,    J = ACELL(JJ)%JNBR(KK)>0,
   !!      KK = ACELL(J)%KNBR(K),      K = ACELL(JJ)%KNBR(KK).
   !!
   !! We implicitly assume that the input mesh has a valid topology and make
   !! no effort here to verify a good topology.
   !!
   !! The identification of neighbors is accomplished using a bin table of all
   !! the active sides.  We use an order-invariant hash of the side nodes to
   !! generate a bin number in the range [0,N) for each side.  Once this data
   !! structure has been created, neighbors can be found by a quick search of
   !! a (very) small bin.  With an ideal hash function, the number of bins N
   !! is between 1/2 and 1 times the number of active sides, but generally much
   !! closer to the lower bound.  Setting N equal to the number of active sides
   !! is generous.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !! Allocate space for the bin table, and initialize the hash function.
    n = sum(popcnt(acell%smask))  ! number of active sides
    allocate(table(n))
    call hpar%init (n, mesh%num_node)  ! adjusts N up to a power of 2
    allocate(xbin(0:n))
    
    !! Count the number of hits to each bin; store the count for bin N in XBIN(N+1).
    xbin = 0
    do j = 1, size(acell)
      do k = 1, acell(j)%nface
        if (btest(acell(j)%smask,k-1)) then
          call get_face_nodes (acell(j)%connect, k, side1)
          call hpar%hash (side1, n)
          xbin(n+1) = xbin(n+1) + 1
        end if
      end do
    end do

    !! Prepare XBIN: bin J will be TABLE(XBIN(J):XBIN(J+1)-1)
    xbin(0) = 1
    do j = 1, ubound(xbin,1)
      xbin(j) = xbin(j) + xbin(j-1)
    end do

    !! Fill the table; use XBIN as a temporary to hold the next free location for each bin.
    do j = 1, size(acell)
      do k = 1, acell(j)%nface
        if (btest(acell(j)%smask,k-1)) then
          call get_face_nodes (acell(j)%connect, k, side1, normalize=.true.)
          call hpar%hash (side1, n)
          i = xbin(n)
          table(i)%j = j
          table(i)%k = k
          call move_alloc (side1, table(i)%side)
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

    !! The bin table is complete.  Now use it to find the neighbor across each
    !! active side, filling in ACELL%JNBR and ACELL%KNBR.
    do j = 1, size(acell)
      do k = 1, acell(j)%nface
        if (btest(acell(j)%smask,k-1)) then
          if (acell(j)%jnbr(k) > 0) cycle  ! info already assigned
          !! The SIDE1 my neighbor would have for this side.
          call get_face_nodes (acell(j)%connect, k, side1, normalize=.true., reverse=.true.)
          call hpar%hash (side1, n)
          associate (bin => table(xbin(n):xbin(n+1)-1))
            !! Search the bin for a matching SIDE1.
            do i = size(bin), 1, -1
              if (all(side1 == bin(i)%side)) exit  ! found my neighbor
            end do
            if (i > 0) then
              !! My neighbor.
              jj = bin(i)%j
              kk = bin(i)%k
              acell(j)%jnbr(k) = jj
              acell(j)%knbr(k) = kk
              !! My neighbor's neighbor (me).
              acell(jj)%jnbr(kk) = j
              acell(jj)%knbr(kk) = k
            else  ! no neighbor -- must be a boundary face
              acell(j)%jnbr(k) = 0
              acell(j)%knbr(k) = 0
            end if
          end associate
        end if
      end do
    end do

    !! Find the neighboring cell to each active link side, filling in
    !! LINK_JNBR and LINK_KNBR; 0 values assigned for inactive sides.
    allocate(link_jnbr(2,link%nlink), link_knbr(2,link%nlink))
    do j = 1, link%nlink
      associate (link_nodes => link%lnode(link%xlnode(j):link%xlnode(j+1)-1))
        n = size(link_nodes)/2  ! number of nodes on each link face
        do k = 1, 2
          !! The SIDE1 my neighbor cell would have for this side.
          select case (k)
          case (1)
            side1 = link_nodes(:n)
          case (2)
            side1 = link_nodes(n+1:)
            call reverse_facet (side1)
          end select
          call normalize_facet (side1)
          if (any(active_node(side1))) then
            call hpar%hash (side1, n)
            associate (bin => table(xbin(n):xbin(n+1)-1))
              !! Search the bin for the matching SIDE1.
              do i = size(bin), 1, -1
                if (all(side1 == bin(i)%side)) exit  ! found my neighbor
              end do
              INSIST(i /= 0)
              link_jnbr(k,j) = bin(i)%j
              link_knbr(k,j) = bin(i)%k
            end associate
          else
            link_jnbr(k,j) = 0
            link_knbr(k,j) = 0
          end if
        end do
      end associate
    end do

    deallocate(table, xbin, active_node)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Generate the component array ACELL%ISIDE which identifies the interface
   !! sides and other active sides of active cells.  When positive, ISIDE(K)
   !! is the index of the side set that specified side K of the cell as an
   !! interface side.  ISIDE = 0 for remaining active, non-interface sides,
   !! and ISIDE < 0 for inactive sides to be ignored.  The array data satisfies
   !! the reciprocity relation ACELL(J)%ISIDE(K) = ACELL(JJ)%ISIDE(KK) whenever
   !! JJ = ACELL(J)%JNBR(K) > 0, KK = ACELL(J)%KNBR(K,J).  This means that
   !! neighboring cells have the same ISIDE value for their shared side.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! Generate the cell-to-active-cell number map.
    allocate(map(mesh%num_elem))
    map = 0
    do j = 1, size(acell)
      map(acell(j)%cnum) = j
    end do

    !! Tag the side-set-specified interface sides with the side set index.
    do n = 1, size(mesh%sset)
      if (.not.active_sset(n)) cycle
      do i = 1, mesh%sset(n)%num_side
        k = mesh%sset(n)%face(i)
        j = map(mesh%sset(n)%elem(i))
        ASSERT(j > 0)
        if (acell(j)%iside(k) == 0) then
          acell(j)%iside(k) = n
        else if (acell(j)%iside(k) /= n) then
          stat = -1
          errmsg = 'CREATE_INTERNAL_INTERFACES: overlapping side sets: ID=' // &
                   i_to_c(mesh%sset(n)%ID) // ',' // i_to_c(mesh%sset(acell(j)%iside(k))%ID)
          return
        end if
      end do
    end do
    deallocate(map)

    !! Examine the ACELL%ISIDE array for consistency.  Interface sides cannot
    !! lie on the mesh boundary.  If tagged, the neighboring side to an
    !! interface side must be tagged with the same side set index.  If the
    !! neighboring side is not tagged, we tag it with the side set index and add
    !! the side to collection of sides to be added (optionally) to the side set.

    allocate(extra_sides(size(mesh%sset)))
    do j = 1, size(acell)
      do k = 1, acell(j)%nface
        n = acell(j)%iside(k)
        if (n > 0) then
          kk = acell(j)%knbr(k)
          jj = acell(j)%jnbr(k)
          if (jj <= 0) then
            stat = -1
            errmsg = 'CREATE_INTERNAL_INTERFACES: side set lies on boundary: ID=' // &
                     i_to_c(mesh%sset(n)%ID)
            return
          end if
          nn = acell(jj)%iside(kk)
          if (nn == 0) then
            acell(jj)%iside(kk) = n
            call extra_sides(n)%add (6*acell(jj)%cnum + kk)
          else if (nn /= n) then
            stat = -1
            errmsg = 'CREATE_INTERNAL_INTERFACES: overlapping side sets: ID=' // &
                     i_to_c(mesh%sset(n)%ID) // ',' // i_to_c(mesh%sset(nn)%ID)
            return
          end if
        end if
      end do
    end do

    !! Tag all inactive sides with a negative value in ACELL%ISIDE.
    do j = 1, size(acell)
      do k = 1, acell(j)%nface
        if (.not.btest(acell(j)%smask,k-1)) acell(j)%iside(k) = -1
      end do
    end do

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Each 'side' in a side set is a cell/side index pair.  The side sets that
   !! specify the interface surface may be one-sided (cells from one side of
   !! the surface), two-sided (cells from both sides of the surface), or
   !! something in between.  It may be desirable in the output mesh, however,
   !! for a side set to include 'sides' from both sides of the interface it
   !! specified.  For this purpose, the preceding stage also generated the
   !! sets of sides EXTRA_SIDES that are not present in side sets: if N is a
   !! value in the integer set EXTRA_SIDES(I), N/6 is the element index and
   !! MODULO(N,6) is the side index of a side to be added to side set I. Only
   !! when an interface side set is two-sided will there be no missing sides.
   !!
   !! This section of code is optional.  I have left it in for the time being
   !! even though I have currently by-passed the code that actually modifies
   !! the mesh side set info, thinking it might be desired at some later time.
   !! To eliminate it, drop anything having to do with EXTRA_SIDES in the code
   !! above and delete the following section of code.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! For now we'll leave the side sets as originally defined.
    if (.false.) then
      !! Modify the mesh side set data to include the 'missing' sides.
      if (allocated(itmp)) deallocate(itmp)
      do i = 1, size(mesh%sset)
        if (extra_sides(i)%size() > 0) then
          sides = extra_sides(n)
          n = mesh%sset(i)%num_side
          mesh%sset(i)%num_side = mesh%sset(i)%num_side + size(sides)
          !! Extend the ELEM array component with the extra side data.
          call move_alloc(mesh%sset(i)%elem, itmp)
          allocate(mesh%sset(i)%elem(mesh%sset(i)%num_side))
          mesh%sset(i)%elem(:n) = itmp
          mesh%sset(i)%elem(n+1:) = sides/6
          deallocate(itmp)
          !! Extend the FACE array component with the extra side data.
          call move_alloc(mesh%sset(i)%face, itmp)
          allocate(mesh%sset(i)%face(mesh%sset(i)%num_side))
          mesh%sset(i)%face(:n) = itmp
          mesh%sset(i)%face(n+1:) = modulo(sides,6)
          deallocate(itmp)
        end if
      end do
    end if
    deallocate(extra_sides)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Slice open the mesh by replacing each active node of each active cell by
   !! a new daughter node.  This disconnects a cell from its neighbors across
   !! active sides.  The daughter nodes are numbered consecutively following
   !! the nodes of the input mesh.  PARENT(J) records the original node number
   !! for daughter node J.  For simplicity we extend the array to all nodes by
   !! setting PARENT(J) = J for original nodes.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    n = sum(popcnt(acell%nmask))
    allocate(parent(n+mesh%num_node))
    do j = 1, mesh%num_node
      parent(j) = j
    end do
    parent(mesh%num_node+1:) = 0

    !! Replace each active node of each active cell with a new daughter node.
    n = mesh%num_node
    do j = 1, size(acell)
      do k = 1, size(acell(j)%connect)
        if (btest(acell(j)%nmask,k-1)) then
          n = n + 1
          parent(n) = acell(j)%connect(k)
          acell(j)%connect(k) = n
        end if
      end do
    end do

    ASSERT(n == size(parent))
    ASSERT(minval(parent) > 0 .and. maxval(parent) <= mesh%num_node)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Reconnect existing links to their neighboring cells.  Some of the
   !! collateral damage done by the preceding step may have been to
   !! disconnect them.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do j = 1, link%nlink
      associate (link_nodes => link%lnode(link%xlnode(j):link%xlnode(j+1)-1))
        n = size(link_nodes)/2  ! number of nodes on each link face
        do k = 1, 2 ! loop over the two link sides
          if (link_jnbr(k,j) > 0) then
            jj = link_jnbr(k,j)
            kk = link_knbr(k,j)
            !! The co-oriented side node lists for the link and neighbor cell.
            select case (k)
            case (1)
              link_side => link_nodes(:n)
              call get_face_nodes (acell(jj)%connect, kk, cell_side)
            case (2)
              link_side => link_nodes(n+1:)
              call get_face_nodes (acell(jj)%connect, kk, cell_side, reverse=.true.)
            end select
            !! Rotate the CELL_SIDE list so that it is in 1-1 correspondence with LINK_SIDE.
            shift = 0
            do while (parent(cell_side(1+shift)) /= link_side(1))
              shift = 1 + shift
              ASSERT(shift < size(link_side))
            end do
            cell_side = cshift(cell_side, shift)
            ASSERT(all(parent(cell_side) == link_side))
            !! Reconnect the link to its neighbor cell.
            link_side = cell_side
          end if
        end do
      end associate
    end do
    
    deallocate(link_jnbr, link_knbr)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! A preceding step sliced open the mesh along all active sides.  On
   !! interface sides we now need to create a link between the two cells that
   !! were disconnected from each other.  Remaining active sides, which
   !! shouldn't have been sliced open, need to have the two cells reconnected
   !! by equivalencing matching pairs of nodes on the shared side.
   !!
   !! The EQUIV array is a linked-list data structure describing the state
   !! of the evolving node equivalence classes.  If EQUIV(J) equals 0, node J
   !! is the representative for its class.  Otherwise node J is equivalent
   !! to node EQUIV(J), and the representative for its class can be found
   !! by walking down the linked list to the end.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! Define the list of link block IDs.  The side set ID used to define the
    !! interface is inherited as the link block ID of the corresponding links.
    link%nlblock = link%nlblock + count(active_sset)
    call move_alloc (link%link_block_id, itmp)
    allocate(link%link_block_id(link%nlblock))
    link%link_block_id(:size(itmp)) = itmp
    link%link_block_id(size(itmp)+1:) = pack(mesh%sset(:)%ID, mask=active_sset)
    deallocate(itmp, active_sset)

    !! Extend the link data structure storage to accomodate the new links.
    n1 = 0  ! number of links
    n2 = 0  ! number of link nodes
    do j = 1, size(acell)
      do k = 1, acell(j)%nface
        if (acell(j)%iside(k) > 0) then
          n1 = n1 + 1
          n2 = n2 + acell(j)%side_size(k)
        end if
      end do
    end do
    n1 = n1 / 2 ! double counted
    !! Extend link_block
    call move_alloc (link%link_block, itmp)
    allocate(link%link_block(size(itmp)+n1))
    link%link_block(:size(itmp)) = itmp
    deallocate(itmp)
    !! Extend xlnode
    call move_alloc (link%xlnode, itmp)
    allocate(link%xlnode(size(itmp)+n1))
    link%xlnode(:size(itmp)) = itmp
    deallocate(itmp)
    !! Extend lnode
    call move_alloc (link%lnode, itmp)
    allocate(link%lnode(size(itmp)+n2))
    link%lnode(:size(itmp)) = itmp
    deallocate(itmp)

    allocate(equiv(size(parent)))
    equiv = 0

    n = link%nlink  ! still the old value
    link%nlink = size(link%xlnode) - 1
    do j = 1, size(acell)
      do k = 1, acell(j)%nface

        if (acell(j)%iside(k) < 0) cycle ! not an active side
        if (j > acell(j)%jnbr(k))  cycle ! handle this one from the other side

        jj = acell(j)%jnbr(k)
        kk = acell(j)%knbr(k)

        !! The co-oriented side node lists for the neighboring elements.
        call get_face_nodes (acell(j)%connect, k, side1)
        call get_face_nodes (acell(jj)%connect, kk, side2, reverse=.true.)
        ASSERT(allocated(side1) .and. allocated(side2))
        ASSERT(size(side1) == size(side2))

        !! Rotate the SIDE2 list so that it is in 1-1 correspondence with SIDE1.
        shift = 0
        do while (parent(side2(1+shift)) /= parent(side1(1)))
          shift = 1 + shift
          ASSERT( shift < size(side1) )
        end do
        side2 = cshift(side2, shift)
        ASSERT(all(parent(side1) == parent(side2)))

        select case (acell(j)%iside(k))
        case (1:) ! interface side

          !! Create a link between the adjacent cells across this interface side.
          !! The link inherits the defining side set ID as its block ID.
          n = n + 1
          link%xlnode(n+1) = link%xlnode(n) + size(side1) + size(side2)
          link%lnode(link%xlnode(n):link%xlnode(n)+size(side1)-1) = side1
          link%lnode(link%xlnode(n)+size(side1):link%xlnode(n+1)-1) = side2
          link%link_block(n) = mesh%sset(acell(j)%iside(k))%id

        case (0) ! open side

          !! Reconnect the neighboring cells: equivalence matching nodes of the two sides.
          do i = 1, size(side1)
            l1 = rep_node(side1(i)) ! representative of SIDE1(I)'s class
            l2 = rep_node(side2(i)) ! representative of SIDE2(I)'s class
            if (l1 /= l2) equiv(l1) = l2
          end do

        end select

      end do
    end do
    ASSERT(n == link%nlink)
    ASSERT(size(link%lnode) == link%xlnode(link%nlink+1)-1)

    ASSERT(all(equiv(:mesh%num_node) == 0))
    ASSERT(minval(equiv,equiv/=0) > mesh%num_node)
    ASSERT(maxval(equiv,equiv/=0) <= ubound(equiv,1))

    !! Count the number of daughter nodes arising from each active node,
    !! using the unused initial part of the EQUIV array as a temporary.
    !! This effectively tags active nodes with a positive value in EQUIV.
    !! The actual count may be useful at some point.
    do j = mesh%num_node+1, ubound(parent,1)
      if (equiv(j) == 0) equiv(parent(j)) = 1 + equiv(parent(j))
    end do

    !! Reuse each (now defunct) interface node by making it the representative
    !! for one of its daughter nodes; restore the EQUIV array as we go.
    do j = mesh%num_node+1, ubound(parent,1)
      if (equiv(j) /= 0) cycle
      if (equiv(parent(j)) > 0) then  ! the parent hasn't been reused yet
        equiv(j) = parent(j)
        equiv(parent(j)) = 0  ! mark the parent as used
      end if
    end do

    ASSERT(all(equiv(:mesh%num_node) == 0))

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Generate a sequential node numbering for the modified mesh.  All original
   !! nodes retain their number, including the interface nodes that took the
   !! place of one of their daughters.  The remaining new daughter nodes are
   !! numbered sequentially following the others.  The cell node connectivity
   !! and link data structure must be updated to reflect the new numbering.
   !!
   !! The array MAP defines the mapping from the current temporary numbering
   !! to final node numbers; it is a many-to-one mapping because of the
   !! equivalencing of nodes done in the preceding step.
   !!
   !! The final step generates the interface node daughter structure.  When
   !! J <= mesh%num_node, DAUGHTER(J) > 0 identifies node J as an interface node.
   !! Its first daughter replaces it as node J, and its next daughter is node
   !! DAUGHTER(J).  The complete set of daughter nodes is obtained by walking
   !! down the linked-list until a zero DAUGHTER value is encountered.  All
   !! the daughters after the first are new nodes, and are numbered after the
   !! nodes in the input mesh.  DAUGHTER(J) = 0 for all original inactive
   !! nodes.  This essentially the one-to-many mapping of input mesh node
   !! numbers to output node numbers.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    new_nnode = count(equiv == 0)
    allocate(map(ubound(equiv,1)))
    do j = 1, mesh%num_node
      map(j) = j
    end do

    !! Consecutive numbering of representative nodes.
    n = mesh%num_node
    do j = mesh%num_node+1, ubound(map,1)
      if (equiv(j) /= 0) cycle
      n = n + 1
      map(j) = n
    end do
    ASSERT(n == new_nnode)

    !! Number non-representative nodes as their representative.
    do j = mesh%num_node+1, ubound(map,1)
      if (equiv(j) /= 0) map(j) = map(rep_node(j))
    end do

    !! Renumber the active cell nodes -- the only ones that will have changed.
    !! Note that this modifies MESH%EBLK(:)%CONNECT(:,:), as desired.
    do j = 1, size(acell)
      associate (connect => acell(j)%connect)
        do k = 1, size(connect)
          if (btest(acell(j)%nmask,k-1)) connect(k) = map(connect(k))
        end do
      end associate
    end do
    
    !! Renumber the link nodes.
    do i = 1, size(link%lnode)
      link%lnode(i) = map(link%lnode(i))
    end do

    !! Define the DAUGHTER data structure.
    allocate(daughter(new_nnode))
    daughter = 0
    do j = ubound(parent,1), 1+mesh%num_node, -1
      if (equiv(j) /= 0) cycle
      daughter(map(j)) = daughter(parent(j))
      daughter(parent(j)) = map(j)
    end do
    ASSERT(minval(daughter,daughter/=0) > mesh%num_node)
    ASSERT(maxval(daughter) <= size(daughter))

    deallocate(equiv, parent, map)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Generate the node coordinate array for the modified mesh.  Inactive nodes
   !! keep their original location and all daughter nodes are coincident with
   !! their parent node.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! Generate the new-to-old node mapping array (onto, many-to-one).
    allocate(map(new_nnode))
    map = 0
    do j = 1, mesh%num_node
      n = daughter(j)
      do while (n /= 0)
        ASSERT(map(n) == 0)
        map(n) = j
        n = daughter(n)
      end do
      map(j) = j
    end do
    ASSERT(minval(map) > 0)
    ASSERT(maxval(map) <= mesh%num_node)

    !! Define the node location array of the modified mesh;
    !! daughter nodes inherit the location of their parent.
    allocate(new_x(size(mesh%coord,dim=1),new_nnode))
    do j = 1, new_nnode
      new_x(:,j) = mesh%coord(:,map(j))
    end do
    deallocate(mesh%coord)
    call move_alloc (new_x, mesh%coord)
    mesh%num_node = new_nnode

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

  end subroutine create_internal_interfaces

end module more_mesh_modification
