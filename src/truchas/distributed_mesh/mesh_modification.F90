!!
!! MESH_MODIFICATION
!!
!! This module provides several routines that apply specialized schemes that
!! modify a primitive mesh (in its intermediate serial EXTERNAL_MESH-type
!! form) prior to generating the final distributed mesh.
!!
!! Neil N. Carlson <nnc@lanl.gov>
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

module mesh_modification

  use kinds, only: r8
  use mesh_importer, only: external_mesh
  use string_utilities, only: i_to_c
  implicit none
  private

  public :: convert_cells_to_links
  public :: create_internal_interfaces

contains

  subroutine convert_cells_to_links (mesh, block_id, stat, errmsg)

    use mesh_importer, only: external_mesh, side_set
    use cell_topology
    use hashing
    use string_utilities, only: i_to_c

    type(external_mesh), intent(inout) :: mesh
    integer, intent(in) :: block_id(:)
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    integer :: i, j, k, n, side_bitmask, node_bitmask
    logical, allocatable :: link_cell(:), link_node(:), keep(:)
    integer, allocatable :: link_side(:), xbin(:), map(:)
    integer, pointer :: cnode(:), fnode(:), old_cnode(:,:), old_cell_block(:), tmp(:)
    type(side_set), pointer :: old_sset(:)
    type(hash_param) :: hpar

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
    call initialize_hash_param (hpar, n, mesh%nnode)  ! adjusts N to a power of 2
    allocate(xbin(0:n))

    !! Count the number of hits to each bin; store the count for bin N in XBIN(N+1).
    xbin = 0
    do j = 1, mesh%ncell
      if (link_side(j) == 0) cycle
      do k = 1, NUM_FACE
        if (btest(link_side(j), k-1)) then
          fnode => hex_face_nodes (mesh%cnode(:,j), k, normalize=.true.)
          call hash (hpar, fnode, n)
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
          call hash (hpar, fnode, n)
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
          call hash (hpar, fnode, n)
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
    old_cnode => mesh%cnode
    old_cell_block => mesh%cell_block
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
      tmp => mesh%block_id
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
      tmp => mesh%sset(i)%elem
      allocate(mesh%sset(i)%elem(n))
      mesh%sset(i)%elem = map(pack(tmp, mask=keep))
      deallocate(tmp)
      tmp => mesh%sset(i)%face
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

  subroutine create_internal_interfaces (mesh, ssid, stat, errmsg)

    use hashing
    use cell_topology
    use mesh_importer, only: side_set_node_list

    type(external_mesh), intent(inout) :: mesh
    integer, intent(in) :: ssid(:)
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    integer :: num_cell_side, num_side_node
    integer :: i, j, k, n, l1, l2, bit_mask, shift, jj, kk, nn, new_nnode
    integer, pointer :: list(:), side_sig(:), side1(:), side2(:)
    logical, allocatable :: active_sset(:), active_node(:), active_cell(:)
    integer, allocatable :: parent(:), equiv(:), map(:), daughter(:)
    integer, allocatable :: xside(:), side(:,:)
    integer, pointer :: iptr1(:), iptr2(:,:)
    real(r8), pointer :: new_x(:,:)

    !! Data arrays for the active submesh.
    integer :: num_active_cell
    integer, allocatable :: active_cnum(:), active_nmask(:), active_smask(:)
    integer, allocatable :: active_jnbr(:,:), active_knbr(:,:), active_iside(:,:)
    integer, allocatable :: link_jnbr(:,:), link_knbr(:,:)

    !! Local variables used to compute the neighbor data.
    type :: table_entry
      integer :: j, k
      integer, pointer :: side(:) => null()
    end type table_entry
    type(table_entry), pointer :: table(:), bin(:)
    integer, allocatable :: xbin(:)
    type(hash_param) :: hpar

    stat = 0
    errmsg = ''

    if (size(ssid) == 0) return ! nothing to do
    
    !! Define some cell-type dependent constants.
    select case (mesh%mesh_type)
    case ('TET')
      num_cell_side = 4
      num_side_node = 3
      side_sig => TETRA4_FACE_SIG
    case ('HEX')
      num_cell_side = 6
      num_side_node = 4
      side_sig => HEX8_FACE_SIG
    case default
      stat = -1
      errmsg = 'CREATE_INTERNAL_INTERFACES: unable to handle mesh type: ' // trim(mesh%mesh_type)
      return
    end select

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
    allocate(active_node(mesh%nnode))
    active_node = .false.
    do n = 1, size(mesh%sset)
      if (.not.active_sset(n)) cycle
      list => side_set_node_list(mesh, n)
      do i = 1, size(list)
        active_node(list(i)) = .true.
      end do
      deallocate(list)
    end do

    !! Generate the active cell mask array ACTIVE_CELL.
    allocate(active_cell(mesh%ncell))
    do j = 1, mesh%ncell
      active_cell(j) = any(active_node(mesh%cnode(:,j)))
    end do

    !! Generate the list of active cells ACTIVE_CNUM and the corresponding
    !! bitmask arrays ACTIVE_NMASK and ACTIVE_SMASK:
    !! BTEST(ACTIVE_NMASK(J),K-1) is true when node K of active cell J is an active node
    !! BTEST(ACTIVE_SMASK(j),K-1) is true when side K of active cell J is an active side.
    n = count(active_cell)
    allocate(active_cnum(n), active_nmask(n), active_smask(n))
    num_active_cell = n
    n = 0
    do j = 1, mesh%ncell
      if (.not.active_cell(j)) cycle
      n = n + 1
      active_cnum(n) = j
      !! Define NMASK: mark the active nodes on this active cell.
      list => mesh%cnode(:,j)
      bit_mask = 0
      do k = 1, size(list)
        if (active_node(list(k))) bit_mask = ibset(bit_mask, k-1)
      end do
      active_nmask(n) = bit_mask
      !! Define SMASK: mark the active sides on this active cell.
      bit_mask = 0
      do k = 1, size(side_sig)
        if (iand(active_nmask(n), side_sig(k)) /= 0) bit_mask = ibset(bit_mask, k-1)
      end do
      active_smask(n) = bit_mask
    end do

    deallocate(active_cell)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Generate the neighbor data arrays ACTIVE_JNBR and ACTIVE_KNBR.  When
   !! positive, ACTIVE_JNBR(K,J) is the active cell that neighbors active
   !! cell J across side K, and ACTIVE_KNBR(K,J) is the corresponding side of
   !! the neighbor.  Boundary sides are indicated by a zero ACTIVE_JNBR value,
   !! and inactive sides to be ignored are indicated by a negative value.
   !! These arrays satisfy the reciprocity relations
   !!
   !!      JJ = ACTIVE_JNBR(K,J)>0,    J = ACTIVE_JNBR(KK,JJ)>0,
   !!      KK = ACTIVE_KNBR(K,J),      K = ACTIVE_KNBR(KK,JJ).
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

    !! Count the number of active sides N.
    n = 0
    do j = 1, num_active_cell
      n = n + bit_count(active_smask(j))
    end do

    !! Allocate space for the bin table, and initialize the hash function.
    allocate(table(n))
    call initialize_hash_param (hpar, n, mesh%nnode)  ! adjusts N up to a power of 2
    allocate(xbin(0:n))

    !! Count the number of hits to each bin; store the count for bin N in XBIN(N+1).
    xbin = 0
    do j = 1, num_active_cell
      do k = 1, num_cell_side
        if (btest(active_smask(j),k-1)) then
          side1 => side_node_list(mesh%cnode(:,active_cnum(j)), k)
          call hash (hpar, side1, n)
          xbin(n+1) = xbin(n+1) + 1
          deallocate(side1)
        end if
      end do
    end do

    !! Prepare XBIN: bin J will be TABLE(XBIN(J):XBIN(J+1)-1)
    xbin(0) = 1
    do j = 1, ubound(xbin,1)
      xbin(j) = xbin(j) + xbin(j-1)
    end do

    !! Fill the table; use XBIN as a temporary to hold the next free location for each bin.
    do j = 1, num_active_cell
      do k = 1, num_cell_side
        if (btest(active_smask(j),k-1)) then
          side1 => side_node_list(mesh%cnode(:,active_cnum(j)), k, normalize=.true.)
          call hash (hpar, side1, n)
          i = xbin(n)
          table(i)%side => side1
          table(i)%j = j
          table(i)%k = k
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
    !! active side, filling in ACTIVE_JNBR and ACTIVE_KNBR.
    allocate(active_jnbr(num_cell_side,num_active_cell), active_knbr(num_cell_side,num_active_cell))
    active_jnbr = -1
    active_knbr = -1
    do j = 1, num_active_cell
      do k = 1, num_cell_side
        if (btest(active_smask(j),k-1)) then
          if (active_jnbr(k,j) > 0) cycle  ! info already assigned
          !! The SIDE1 my neighbor would have for this side.
          side1 => side_node_list(mesh%cnode(:,active_cnum(j)), k, normalize=.true., reverse=.true.)
          call hash (hpar, side1, n)
          bin => table(xbin(n):xbin(n+1)-1)
          !! Search the bin for a matching SIDE1.
          do i = size(bin), 1, -1
            if (all(side1 == bin(i)%side)) exit  ! found my neighbor
          end do
          deallocate(side1)
          if (i > 0) then
            !! My neighbor.
            jj = bin(i)%j
            kk = bin(i)%k
            active_jnbr(k,j) = jj
            active_knbr(k,j) = kk
            !! My neighbor's neighbor (me).
            active_jnbr(kk,jj) = j
            active_knbr(kk,jj) = k
          else  ! no neighbor -- must be a boundary face
            active_jnbr(k,j) = 0
            active_knbr(k,j) = 0
          end if
        end if
      end do
    end do

    !! Find the neighboring cell to each active link side, filling in
    !! LINK_JNBR and LINK_KNBR; 0 values assigned for inactive sides.
    if (associated(mesh%lnode)) then
      allocate(side1(num_side_node))
      allocate(link_jnbr(2,mesh%nlink), link_knbr(2,mesh%nlink))
      do j = 1, mesh%nlink
        do k = 1, 2
          !! The SIDE1 my neighbor cell would have for this side.
          select case (k)
          case (1)
            side1 = mesh%lnode(:num_side_node,j)
          case (2)
            side1 = mesh%lnode(num_side_node+1:,j)
            call reverse_facet (side1)
          end select
          call normalize_facet (side1)
          if (any(active_node(side1))) then
            call hash (hpar, side1, n)
            bin => table(xbin(n):xbin(n+1)-1)
            !! Search the bin for the matching SIDE1.
            do i = size(bin), 1, -1
              if (all(side1 == bin(i)%side)) exit  ! found my neighbor
            end do
            INSIST(i /= 0)
            link_jnbr(k,j) = bin(i)%j
            link_knbr(k,j) = bin(i)%k
          else
            link_jnbr(k,j) = 0
            link_knbr(k,j) = 0
          end if
        end do
      end do
      deallocate(side1)
    else
      allocate(link_jnbr(2,0), link_knbr(2,0))
    end if

    !! We're finished with the bin table now.
    do j = 1, size(table)
      deallocate(table(j)%side)
    end do
    deallocate(table, xbin)

    deallocate(active_node)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Generate the ACTIVE_ISIDE array which identifies the interface sides and
   !! other active sides of active cells.  When positive, ACTIVE_ISIDE(K,J) is
   !! the index of the side set that specified side K of active cell J as an
   !! interface side.  ACTIVE_ISIDE = 0 for remaining active, non-interface
   !! sides, and ISIDE is negative for inactive sides to be ignored.  The array
   !! satisfies the reciprocity relation
   !!
   !!   ACTIVE_ISIDE(K,J) = ACTIVE_ISIDE(KK,JJ) when
   !!   JJ = ACTIVE_JNBR(K,J) > 0, KK = ACTIVE_KNBR(K,J).
   !!
   !! This means that neighboring cells have the same ACTIVE_ISIDE value for
   !! their shared side.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(active_iside(num_cell_side,num_active_cell))
    active_iside = 0

    !! Generate the cell-to-active-cell number map.
    allocate(map(mesh%ncell))
    map = 0
    do j = 1, num_active_cell
      map(active_cnum(j)) = j
    end do

    !! Tag the side-set-specified interface sides with the side set index.
    do n = 1, size(mesh%sset)
      if (.not.active_sset(n)) cycle
      do i = 1, mesh%sset(n)%num_side
        k = mesh%sset(n)%face(i)
        j = map(mesh%sset(n)%elem(i))
        ASSERT(j>0)
        if (active_iside(k,j) == 0) then
          active_iside(k,j) = n
        else if (active_iside(k,j) /= n) then
          stat = -1
          errmsg = 'CREATE_INTERNAL_INTERFACES: overlapping side sets: ID=' // &
                   i_to_c(mesh%sset(n)%ID) // ',' // i_to_c(mesh%sset(active_iside(k,j))%ID)
          return
        end if
      end do
    end do
    deallocate(map)

    !! Examine the ACTIVE_ISIDE array for consistency.  Interface sides cannot
    !! lie on the mesh boundary.  If tagged, the neighboring side to an
    !! interface side must be tagged with the same side set index.  If the
    !! neighboring side is not tagged, we tag it temporarily with the negative
    !! of the side set index (see below) but will unnegate the value later to
    !! symmetrize ACTIVE_ISIDE.

    allocate(xside(1+size(mesh%sset)))
    xside = 0 ! per side set count of 'missing' sides
    do j = 1, num_active_cell
      do k = 1, num_cell_side
        n = active_iside(k,j)
        if (n > 0) then
          kk = active_knbr(k,j)
          jj = active_jnbr(k,j)
          if (jj <= 0) then
            stat = -1
            errmsg = 'CREATE_INTERNAL_INTERFACES: side set lies on boundary: ID=' // &
                     i_to_c(mesh%sset(n)%ID)
            return
          end if
          nn = active_iside(kk,jj)
          if (nn == 0) then
            active_iside(kk,jj) = -n    ! negative as a special tag for later use
            xside(n+1) = 1 + xside(n+1) ! count of missing sides for side set N
          else if (nn /= n) then
            stat = -1
            errmsg = 'CREATE_INTERNAL_INTERFACES: overlapping side sets: ID=' // &
                     i_to_c(mesh%sset(n)%ID) // ',' // i_to_c(mesh%sset(nn)%ID)
            return
          end if
        end if
      end do
    end do

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Each 'side' in a side set is a cell/side index pair.  The side sets that
   !! specify the interface surface may be one-sided (cells from one side of
   !! the surface), two-sided (cells from both sides of the surface), or
   !! something in between.  It may be desirable in the output mesh, however,
   !! for a side set to include 'sides' from both sides of the interface it
   !! specified.  For this purpose we also generate the XSIDE/SIDE structure
   !! that contain any of these sides that were missing from the input mesh
   !! side set: SIDE(I,XSIDE(N):XSIDE(N+1)-1) are the missing cells (I=1) and
   !! corresponding sides (I=2).  Only when an interface side set is two-sided
   !! will there be no missing sides.
   !!
   !! This section of code is optional.  I have left it in for the time being
   !! even though I have currently by-passed the code that actually modifies
   !! the mesh side set info, thinking it might be desired at some later time.
   !! To eliminate it, drop anything having to do with XSIDE in the loop above
   !! and assign "N" instead of "-N" in the assignment to ACTIVE_ISIDE.  Then
   !! the following section of code can be entirely deleted except for the very
   !! last stanza.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! Setup XSIDE: SIDE(:,XSIDE(N):XSIDE(N+1)-1) are the missing sides for side set N.
    !! XSIDE(N+1) currently contains the count of missing sides for side set N.
    xside(1) = 1
    do j = 2, size(xside)
      xside(j) = xside(j) + xside(j-1)
    end do
    ASSERT(all(xside(:size(xside)-1) <= xside(2:)))

    !! Fill-in the missing sides, correcting ACTIVE_ISIDE as we go.  We use XSIDE
    !! as a temporary to hold the next free location in the list for each side set.
    allocate(side(2,xside(size(xside))-1))
    do j = 1, num_active_cell
      do k = 1, num_cell_side
        if (active_iside(k,j) < 0) then
          n = -active_iside(k,j)
          active_iside(k,j) = n
          side(1,xside(n)) = active_cnum(j)
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

    !! For now we'll leave the side sets as originally defined.
    if (.false.) then
      !! Modify the mesh side set data to include the 'missing' sides.
      do i = 1, size(mesh%sset)
        if (xside(i+1) > xside(i)) then
          n = mesh%sset(i)%num_side
          mesh%sset(i)%num_side = n + xside(i+1) - xside(i)
          !! Extend the ELEM array component with the extra side data.
          iptr1 => mesh%sset(i)%elem
          allocate(mesh%sset(i)%elem(mesh%sset(i)%num_side))
          mesh%sset(i)%elem(:n) = iptr1
          mesh%sset(i)%elem(n+1:) = side(1,xside(i):xside(i+1)-1)
          deallocate(iptr1)
          !! Extend the FACE array component with the extra side data.
          iptr1 => mesh%sset(i)%face
          allocate(mesh%sset(i)%face(mesh%sset(i)%num_side))
          mesh%sset(i)%face(:n) = iptr1
          mesh%sset(i)%face(n+1:) = side(2,xside(i):xside(i+1)-1)
          deallocate(iptr1)
        end if
      end do
    end if
    deallocate(xside, side)


    !! We've finished with the above temporary use of ACTIVE_ISIDE;
    !! tag all inactive sides with a negative value in ACTIVE_ISIDE.
    ASSERT(all(active_iside >= 0))
    do j = 1, num_active_cell
      do k = 1, num_cell_side
        if (.not.btest(active_smask(j),k-1)) active_iside(k,j) = -1
      end do
    end do

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Slice open the mesh by replacing each active node of each active cell by
   !! a new daughter node.  This disconnects a cell from its neighbors across
   !! active sides.  The daughter nodes are numbered consecutively following
   !! the nodes of the input mesh.  PARENT(J) records the original node number
   !! for daughter node J.  For simplicity we extend the array to all nodes by
   !! setting PARENT(J) = J for original nodes.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! Count the number of daughter nodes that will be created.
    n = 0
    do j = 1, num_active_cell
      n = n + bit_count(active_nmask(j))
    end do

    allocate(parent(n+mesh%nnode))
    do j = 1, mesh%nnode
      parent(j) = j
    end do
    parent(mesh%nnode+1:) = 0

    !! Replace each active node of each active cell with a new daughter node.
    n = mesh%nnode
    do j = 1, num_active_cell
      list => mesh%cnode(:,active_cnum(j))
      do k = 1, size(list)
        if (btest(active_nmask(j),k-1)) then
          n = n + 1
          parent(n) = list(k)
          list(k) = n
        end if
      end do
    end do

    ASSERT(n == size(parent))
    ASSERT(minval(parent) > 0 .and. maxval(parent) <= mesh%nnode)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Reconnect existing links to their neighboring cells.  Some of the
   !! collateral damage done by the preceding step may have been to
   !! disconnect them.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do j = 1, mesh%nlink
      do k = 1, 2 ! loop over the two link sides
        if (link_jnbr(k,j) > 0) then
          jj = link_jnbr(k,j)
          kk = link_knbr(k,j)
          !! The co-oriented side node lists for the link and neighbor cell.
          select case (k)
          case (1)
            side1 => mesh%lnode(:num_side_node,j)
            side2 => side_node_list(mesh%cnode(:,active_cnum(jj)), kk)
          case (2)
            side1 => mesh%lnode(num_side_node+1:,j)
            side2 => side_node_list(mesh%cnode(:,active_cnum(jj)), kk, reverse=.true.)
          end select
          !! Rotate the SIDE2 list so that it is in 1-1 correspondence with SIDE1.
          shift = 0
          do while (parent(side2(1+shift)) /= side1(1))
            shift = 1 + shift
            ASSERT(shift < size(side1))
          end do
          side2 = cshift(side2, shift)
          ASSERT(all(parent(side2) == side1))
          !! Reconnect the link to its neighbor cell.
          side1 = side2
          deallocate(side2)
        end if
      end do
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

    allocate(equiv(size(parent)))
    equiv = 0

    !! Define the list of link block IDs.  The side set ID used to define the
    !! interface is inherited as the link block ID of the corresponding links.
    list => mesh%link_block_id
    n = mesh%nlblock  !size(list)
    mesh%nlblock = n + count(active_sset)
    allocate(mesh%link_block_id(mesh%nlblock))
    if (associated(list)) then
      mesh%link_block_id(:n) = list
      deallocate(list)
    end if
    mesh%link_block_id(n+1:) = pack(mesh%sset(:)%ID, mask=active_sset)
    deallocate(active_sset)

    !! Extend the link data structure storage to accomodate the new links.
    iptr2 => mesh%lnode
    iptr1 => mesh%link_block
    n = mesh%nlink + count(active_iside > 0) / 2  ! double counted
    allocate(mesh%lnode(2*num_side_node,n), mesh%link_block(n))
    if (associated(iptr2)) then
      mesh%lnode(:,:mesh%nlink) = iptr2
      mesh%link_block(:mesh%nlink) = iptr1
      deallocate(iptr2, iptr1)
    end if

    n = mesh%nlink  ! still the old value
    mesh%nlink = size(mesh%lnode,dim=2)
    do j = 1, num_active_cell
      do k = 1, num_cell_side

        if (active_iside(k,j) < 0) cycle ! not an active side
        if (j > active_jnbr(k,j))  cycle ! handle this one from the other side

        jj = active_jnbr(k,j)
        kk = active_knbr(k,j)

        !! The co-oriented side node lists for the neighboring elements.
        side1 => side_node_list(mesh%cnode(:,active_cnum(j)), k)
        side2 => side_node_list(mesh%cnode(:,active_cnum(jj)), kk, reverse=.true.)
        ASSERT(associated(side1) .and. associated(side2))
        ASSERT(size(side1) == size(side2))

        !! Rotate the SIDE2 list so that it is in 1-1 correspondence with SIDE1.
        shift = 0
        do while (parent(side2(1+shift)) /= parent(side1(1)))
          shift = 1 + shift
          ASSERT( shift < size(side1) )
        end do
        side2 = cshift(side2, shift)
        ASSERT(all(parent(side1) == parent(side2)))

        select case (active_iside(k,j))
        case (1:) ! interface side

          !! Create a link between the adjacent cells across this interface side.
          !! The link inherits the defining side set ID as its block ID.
          n = n + 1
          mesh%lnode(:num_side_node,n) = side1
          mesh%lnode(num_side_node+1:,n) = side2
          mesh%link_block(n) = mesh%sset(active_iside(k,j))%ID

        case (0) ! open side

          !! Reconnect the neighboring cells: equivalence matching nodes of the two sides.
          do i = 1, size(side1)
            l1 = rep_node(side1(i)) ! representative of SIDE1(I)'s class
            l2 = rep_node(side2(i)) ! representative of SIDE2(I)'s class
            if (l1 /= l2) equiv(l1) = l2
          end do

        end select
        deallocate(side1, side2)

      end do
    end do
    ASSERT(n == mesh%nlink)

    ASSERT(all(equiv(:mesh%nnode) == 0))
    ASSERT(minval(equiv,equiv/=0) > mesh%nnode)
    ASSERT(maxval(equiv,equiv/=0) <= ubound(equiv,1))

    !! Count the number of daughter nodes arising from each active node,
    !! using the unused initial part of the EQUIV array as a temporary.
    !! This effectively tags active nodes with a positive value in EQUIV.
    !! The actual count may be useful at some point.
    do j = mesh%nnode+1, ubound(parent,1)
      if (equiv(j) == 0) equiv(parent(j)) = 1 + equiv(parent(j))
    end do

    !! Reuse each (now defunct) interface node by making it the representative
    !! for one of its daughter nodes; restore the EQUIV array as we go.
    do j = mesh%nnode+1, ubound(parent,1)
      if (equiv(j) /= 0) cycle
      if (equiv(parent(j)) > 0) then  ! the parent hasn't been reused yet
        equiv(j) = parent(j)
        equiv(parent(j)) = 0  ! mark the parent as used
      end if
    end do

    ASSERT(all(equiv(:mesh%nnode) == 0))

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
   !! J <= MESH%NNODE, DAUGHTER(J) > 0 identifies node J as an interface node.
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
    do j = 1, mesh%nnode
      map(j) = j
    end do

    !! Consecutive numbering of representative nodes.
    n = mesh%nnode
    do j = mesh%nnode+1, ubound(map,1)
      if (equiv(j) /= 0) cycle
      n = n + 1
      map(j) = n
    end do
    ASSERT(n == new_nnode)

    !! Number non-representative nodes as their representative.
    do j = mesh%nnode+1, ubound(map,1)
      if (equiv(j) /= 0) map(j) = map(rep_node(j))
    end do

    !! Renumber the active cell nodes -- the only ones that will have changed.
    do j = 1, num_active_cell
      list => mesh%cnode(:,active_cnum(j))
      do k = 1, size(list)
        if (btest(active_nmask(j),k-1)) list(k) = map(list(k))
      end do
    end do

    !! Renumber the link nodes.
    do j = 1, mesh%nlink
      mesh%lnode(:,j) = map(mesh%lnode(:,j))
    end do

    !! Define the DAUGHTER data structure.
    allocate(daughter(new_nnode))
    daughter = 0
    do j = ubound(parent,1), 1+mesh%nnode, -1
      if (equiv(j) /= 0) cycle
      daughter(map(j)) = daughter(parent(j))
      daughter(parent(j)) = map(j)
    end do
    ASSERT(minval(daughter,daughter/=0) > mesh%nnode)
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
    do j = 1, mesh%nnode
      n = daughter(j)
      do while (n /= 0)
        ASSERT(map(n) == 0)
        map(n) = j
        n = daughter(n)
      end do
      map(j) = j
    end do
    ASSERT(minval(map) > 0)
    ASSERT(maxval(map) <= mesh%nnode)

    !! Define the node location array of the modified mesh;
    !! daughter nodes inherit the location of their parent.
    allocate(new_x(size(mesh%x,dim=1),new_nnode))
    do j = 1, new_nnode
      new_x(:,j) = mesh%x(:,map(j))
    end do
    deallocate(mesh%x)
    mesh%x => new_x
    mesh%nnode = new_nnode

    deallocate(map, daughter)

    deallocate(active_cnum, active_nmask, active_smask)
    deallocate(active_jnbr, active_knbr, active_iside)

  contains

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

   !!
   !! This auxillary function returns a pointer to the list of the nodes
   !! on a side of a cell.  The cell nodes are specified by CNODES and the
   !! side index by K.  If the size of CNODES is 4 we assume it is a tet
   !! cell and if the size is 8, a hex cell.  The side nodes are ordered
   !! CCW with respect to the outward orientation unless REVERSE is specified
   !! with the value true, in which case they are ordered CCW with respect
   !! to the inward orientation.  If NORMALIZE is specified with value true,
   !! the lowest numbered node is listed first.
   !!

    function side_node_list (cnodes, k, normalize, reverse) result (list)

      use cell_topology

      integer, intent(in) :: cnodes(:)
      integer, intent(in) :: k
      logical, intent(in), optional :: normalize
      logical, intent(in), optional :: reverse
      integer, pointer :: list(:)

      select case (size(cnodes))
      case (4)  ! 4-node cell ==> a tet
        list => tet_face_nodes(cnodes, k, normalize, reverse)
      case (8)  ! 8-node cell ==> a hex
        list => hex_face_nodes(cnodes, k, normalize, reverse)
      case default
        INSIST(.false.)
      end select

    end function side_node_list

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

end module mesh_modification
