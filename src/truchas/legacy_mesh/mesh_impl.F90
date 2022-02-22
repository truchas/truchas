!!
!! MESH_IMPL
!!
!! Implements the legacy MESH structure array and related data and procedures.
!!

#include "f90_assert.fpp"

module mesh_impl

  use common_impl, only: ncells, nnodes, new_mesh
  use var_vector_types, only: int_var_vector
  use index_map_type, only: index_map
  implicit none
  private

  public :: init_mesh_impl
  public :: is_face_ngbr, mesh_collate_vertex

  integer, parameter, public :: DEGENERATE_FACE = -huge(1)

  type(index_map), public :: legacy_cell_ip

  !! Copy of the MESH_CONNECTIVITY type from MESH_MODULE (used components only)
  type, public :: mesh_connectivity
    integer :: ngbr_cell(6) = -1            ! cell number across each face
    integer :: ngbr_cell_orig(6) = -1       ! original cell numbers across each face
    integer :: ngbr_face(6) = -1            ! face number of cell across each face
    integer :: ngbr_vrtx(8) = -1            ! vertex numbers
    integer :: ngbr_vrtx_orig(8) = -1       ! original vertex numbers
    type(int_var_vector) :: ngbr_cells_all  ! all neighbor cells
    type(int_var_vector) :: ngbr_cells_face ! bits to identify face ngbrs
    integer :: cell_shape                   ! indicator for hex, tet, ...
    integer :: cblockid = 0                 ! cell block id
  end type

  !! Replacement for the MESH structure array from MESH_MODULE
  type(mesh_connectivity), pointer, public :: mesh(:) => null()

  !! Originally from MESH_MODULE
  integer, parameter, public :: CELL_TET       = 1
  integer, parameter, public :: CELL_PYRAMID   = 2
  integer, parameter, public :: CELL_PRISM     = 3
  integer, parameter, public :: CELL_HEX       = 4
  integer, parameter, public :: GAP_ELEMENT_1  = 5
  integer, parameter, public :: GAP_ELEMENT_3  = 6
  integer, parameter, public :: GAP_ELEMENT_5  = 7

contains

  pure logical function is_face_ngbr (face_bits, face)
    integer, intent(in) :: face_bits, face
    is_face_ngbr = btest(face_bits, pos=face)
  end function is_face_ngbr

  !! Copy of function from MESH_MODULE
  function mesh_collate_vertex (mesh)
    use common_impl, only: ncells_tot
    use parallel_communication, only: is_IOP, collate
    type(mesh_connectivity), intent(in) :: mesh(:)
    integer, pointer :: mesh_collate_vertex(:,:)
    integer :: k
    allocate(mesh_collate_vertex(8,merge(ncells_tot,0,is_IOP)))
    do k = 1, 8
      call collate (mesh%ngbr_vrtx_orig(k), mesh_collate_vertex(k,:))
    end do
  end function mesh_collate_vertex

  subroutine init_mesh_impl
    call init_mesh
  end subroutine init_mesh_impl

  subroutine init_mesh

    integer :: j
    integer, allocatable :: ngbr_cell_orig(:,:), ngbr_cell(:,:)
    integer, allocatable :: ngbr_vrtx_orig(:,:), ngbr_vrtx(:,:)
    integer, allocatable :: ngbr_face(:,:)

    allocate(mesh(ncells))

    !! MESH%CELL_SHAPE
    call init_cell_shape (mesh%cell_shape)

    !! MESH%NGBR_CELL, MESH%NGBR_CELL_ORIG
    allocate(ngbr_cell_orig(6,ncells), ngbr_cell(6,ncells))
    call init_ngbr_cell (ngbr_cell_orig, ngbr_cell, legacy_cell_ip)
    do j = 1, ncells
      mesh(j)%ngbr_cell      = ngbr_cell(:,j)
      mesh(j)%ngbr_cell_orig = ngbr_cell_orig(:,j)
    end do
    deallocate(ngbr_cell_orig, ngbr_cell)

    !! MESH%NGBR_FACE
    allocate(ngbr_face(6,ncells))
    call init_ngbr_face (ngbr_face)
    do j = 1, ncells
      mesh(j)%ngbr_face = ngbr_face(:,j)
    end do
    deallocate(ngbr_face)

    !! MESH%NGBR_VRTX, MESH%NGBR_VRTX_ORIG
    allocate(ngbr_vrtx_orig(8,ncells), ngbr_vrtx(8,ncells))
    call init_ngbr_vrtx (ngbr_vrtx_orig, ngbr_vrtx)
    do j = 1, ncells
      mesh(j)%ngbr_vrtx      = ngbr_vrtx(:,j)
      mesh(j)%ngbr_vrtx_orig = ngbr_vrtx_orig(:,j)
    end do
    deallocate(ngbr_vrtx_orig, ngbr_vrtx)

    !! MESH%NGBR_CELLS_ALL, MESH%NGBR_CELLS_FACE
    call init_ngbr_cells_all (mesh%ngbr_cells_all, mesh%ngbr_cells_face)

    !! MESH%CBLOCKID
    call init_cblockid (mesh%cblockid)

  end subroutine init_mesh

  subroutine init_cell_shape (cell_shape)

    integer, intent(out) :: cell_shape(:)

    integer :: j

    ASSERT(size(cell_shape) == ncells)

    !! Real cells
    do j = 1, new_mesh%ncell_onP
      select case (new_mesh%xcnode(j+1)-new_mesh%xcnode(j))
      case (4)
        cell_shape(j) = CELL_TET
      case (5)
        cell_shape(j) = CELL_PYRAMID
      case (6)
        cell_shape(j) = CELL_PRISM
      case (8)
        cell_shape(j) = CELL_HEX
      case default
        INSIST(.false.)
      end select
    end do

    !! Gap cells.  Correct  for hex, but incorrect for prism.
    cell_shape(new_mesh%ncell_onP+1:) = GAP_ELEMENT_5

  end subroutine init_cell_shape

  subroutine init_cblockid (cblockid)

    use bitfield_type, only: popcnt, trailz

    integer, intent(out) :: cblockid(:)

    integer :: j, n

    ASSERT(size(cblockid) == ncells)

    !! Real cells
    do j = 1, new_mesh%ncell_onP
      associate (bitmask => new_mesh%cell_set_mask(j))
        INSIST(popcnt(bitmask) == 1)
        cblockid(j) = new_mesh%cell_set_id(trailz(bitmask))
      end associate
    end do

    !! Gap cells
    n = new_mesh%ncell_onP
    do j = 1, new_mesh%nlink_onP
      if (new_mesh%link_cell_id(j) == 0) cycle ! not from gap cell
      associate (bitmask => new_mesh%link_set_mask(j))
        INSIST(popcnt(bitmask) == 1)
        n = n + 1
        cblockid(n) = new_mesh%link_set_id(trailz(bitmask))
      end associate
    end do

  end subroutine init_cblockid

  !! Initialize the passed NGBR_CELL and NGBR_CELL_ORIG arrays, and return the
  !! associated cell index partition object CELL_IP which must be used by most
  !! of the EE_GATHER_* procedures that operate with the returned arguments.
  !! The legacy mesh API handles certain mesh links as cells, namely those that
  !! originally came from special gap cells in the Exodus mesh. Because of this
  !! a separate legacy cell index partition object is required.
  !!
  !! This implementation does not supply neighbor information for gap cells;
  !! every face of a gap cell will appear to be a degenerate face.
  !!
  !! This implementation does, however, provide info for neighbors that are gap
  !! cells (so the neighbor indexing info is not symmetric).  It appears that
  !! this is not necessary (regression tests pass without it), but it has been
  !! retained for the time being.
  !!
  !! NB1: When the neighbor is a link, but not a gap cell link, it is correct
  !! to say that the face is on the boundary and there is no neighbor.  But if
  !! we do this the htvoid3 regression test fails with a convergence failure in
  !! Ubik. It is highly likely due to a disconnected piece of the mesh with no
  !! flow.  For the time being, until this is resolved, we retain the approach
  !! of saying the neighbor is the cell on the other side of the link.  In
  !! essence giving the neighbor information for the stitched up mesh where
  !! the internal interface has been removed.

  subroutine init_ngbr_cell (ngbr_cell_orig, ngbr_cell, cell_ip)

    use common_impl, only: OLD_TET_SIDE_MAP, OLD_PYR_SIDE_MAP, OLD_PRI_SIDE_MAP, OLD_HEX_SIDE_MAP

    integer, intent(out) :: ngbr_cell_orig(:,:), ngbr_cell(:,:)
    type(index_map), intent(out) :: cell_ip

    integer :: i, j, jj, k, n
    integer :: gid(new_mesh%ncell), link_gid(new_mesh%nlink)

    !! Legacy global cell IDs on the mesh.
    call cell_ip%init (ncells)
    do j = 1, new_mesh%ncell_onP
      gid(j) = cell_ip%global_index(j)
    end do
    call new_mesh%cell_ip%gather_offp(gid)

    !! Generate the cell neighbor array (legacy global IDs)
    do j = 1, new_mesh%ncell_onP
      associate (list => new_mesh%cnhbr(new_mesh%xcnhbr(j):new_mesh%xcnhbr(j+1)-1))
        ngbr_cell_orig(:,j) = DEGENERATE_FACE
        do k = 1, size(list)
          if (list(k) > 0) then
            ngbr_cell_orig(k,j) = gid(list(k))
          else
            ngbr_cell_orig(k,j) = 0
          end if
        end do
      end associate
    end do

    !! Legacy global cell IDs on the mesh links.
    n = new_mesh%ncell_onP
    do j = 1, new_mesh%nlink_onP
      if (new_mesh%link_cell_id(j) > 0) then
        n = n + 1
        link_gid(j) = cell_ip%global_index(n)
      else
        link_gid(j) = 0
      end if
    end do
    call new_mesh%link_ip%gather_offp(link_gid)

    !! Fill in neighbor data from gap cells.
    do n = 1, new_mesh%nlink
      do i = 1, 2
        j = new_mesh%lnhbr(i,n)
        if (j > new_mesh%ncell_onP) cycle
        associate (cface => new_mesh%cface(new_mesh%xcface(j):new_mesh%xcface(j+1)-1))
          do k = size(cface), 1, -1
            if (new_mesh%lface(i,n) == cface(k)) exit
          end do
          INSIST(k > 0)
        end associate
        if (new_mesh%link_cell_id(n) > 0) then ! from a gap cell
          ngbr_cell_orig(k,j) = link_gid(n)
        else ! Not wanted, but see NB1 above.
          jj = new_mesh%lnhbr(1+modulo(i,2),n)
          ngbr_cell_orig(k,j) = gid(jj)
        end if
      end do
    end do

    !! Convert to legacy degenerate hex cells.
    do j = 1, new_mesh%ncell_onP
      select case (new_mesh%xcnode(j+1)-new_mesh%xcnode(j))
      case (4)
        ngbr_cell_orig(:,j) = ngbr_cell_orig(OLD_TET_SIDE_MAP,j)
      case (5)
        ngbr_cell_orig(:,j) = ngbr_cell_orig(OLD_PYR_SIDE_MAP,j)
      case (6)
        ngbr_cell_orig(:,j) = ngbr_cell_orig(OLD_PRI_SIDE_MAP,j)
      case (8)
        ngbr_cell_orig(:,j) = ngbr_cell_orig(OLD_HEX_SIDE_MAP,j)
      end select
    end do

    !! We can get away with giving no neighbor info for gap cells.
    ngbr_cell_orig(:,new_mesh%ncell_onP+1:) = DEGENERATE_FACE

    !! Localize the ngbr_cell_orig array.
    ngbr_cell = ngbr_cell_orig
    where (ngbr_cell < 0) ngbr_cell = 0 ! temp convert degen side tags to ignored bndry sides
    call cell_ip%localize_index_array (ngbr_cell)
    where (ngbr_cell_orig < 0) ngbr_cell = ngbr_cell_orig ! restore degen side tags

    !! Convert off-process references to look-aside boundary array references.
    where (ngbr_cell > ncells) ngbr_cell = ncells - ngbr_cell

  end subroutine init_ngbr_cell

  !! Initialize the passed NGBR_FACE array.  This procedure is the companion
  !! to INIT_NGBR_CELL; see the comments preceding it for relevant details.

  subroutine init_ngbr_face (ngbr_face)

    use common_impl, only: NEW_TET_SIDE_MAP, NEW_PYR_SIDE_MAP, NEW_PRI_SIDE_MAP, NEW_HEX_SIDE_MAP

    integer, intent(out) :: ngbr_face(:,:)

    integer :: i1, i2, j1, j2, k1, k2, n, n1, n2, map1(6), map2(6)

    ASSERT(size(ngbr_face,1) == 6)
    ASSERT(size(ngbr_face,2) == ncells)

    ngbr_face = DEGENERATE_FACE
    do j1 = 1, new_mesh%ncell_onP
      select case (new_mesh%xcnode(j1+1)-new_mesh%xcnode(j1))
      case (4)
        map1 = NEW_TET_SIDE_MAP
      case (5)
        map1 = NEW_PYR_SIDE_MAP
      case (6)
        map1 = NEW_PRI_SIDE_MAP
      case (8)
        map1 = NEW_HEX_SIDE_MAP
      case default
        INSIST(.false.)
      end select
      associate (nhbr1 => new_mesh%cnhbr(new_mesh%xcnhbr(j1):new_mesh%xcnhbr(j1+1)-1))
        do k1 = 1, size(nhbr1)
          i1 = map1(k1) ! old mesh side index
          if (ngbr_face(i1,j1) > 0) cycle ! already defined
          j2 = nhbr1(k1)
          if (j2 > 0) then
            select case (new_mesh%xcnode(j2+1)-new_mesh%xcnode(j2))
            case (4)
              map2 = NEW_TET_SIDE_MAP
            case (5)
              map2 = NEW_PYR_SIDE_MAP
            case (6)
              map2 = NEW_PRI_SIDE_MAP
            case (8)
              map2 = NEW_HEX_SIDE_MAP
            case default
              INSIST(.false.)
            end select
            associate (nhbr2 => new_mesh%cnhbr(new_mesh%xcnhbr(j2):new_mesh%xcnhbr(j2+1)-1))
              do k2 = size(nhbr2), 1, -1
                if (nhbr2(k2) == j1) exit
              end do
              INSIST(k2 > 0)
              i2 = map2(k2) ! old mesh side index
              ngbr_face(i1,j1) = i2
              if (j2 <= new_mesh%ncell_onP) ngbr_face(i2,j2) = i1
            end associate
          else  ! boundary face
            ngbr_face(i1,j1) = 0
          end if
        end do
      end associate
    end do

    !! Fill in data from gap cell neighbors
    do n = 1, new_mesh%nlink
      j1 = new_mesh%lnhbr(1,n)
      select case (new_mesh%xcnode(j1+1)-new_mesh%xcnode(j1))
      case (4)
        map1 = NEW_TET_SIDE_MAP
      case (5)
        map1 = NEW_PYR_SIDE_MAP
      case (6)
        map1 = NEW_PRI_SIDE_MAP
      case (8)
        map1 = NEW_HEX_SIDE_MAP
      case default
        INSIST(.false.)
      end select
      associate (nhbr1 => new_mesh%cface(new_mesh%xcface(j1):new_mesh%xcface(j1+1)-1))
        do k1 = size(nhbr1), 1, -1
          if (new_mesh%lface(1,n) == nhbr1(k1)) exit
        end do
        INSIST(k1 > 0)
        i1 = map1(k1)
      end associate
      j2 = new_mesh%lnhbr(2,n)
      select case (new_mesh%xcnode(j2+1)-new_mesh%xcnode(j2))
      case (4)
        map2 = NEW_TET_SIDE_MAP
      case (5)
        map2 = NEW_PYR_SIDE_MAP
      case (6)
        map2 = NEW_PRI_SIDE_MAP
      case (8)
        map2 = NEW_HEX_SIDE_MAP
      case default
        INSIST(.false.)
      end select
      associate (nhbr2 => new_mesh%cface(new_mesh%xcface(j2):new_mesh%xcface(j2+1)-1))
        do k2 = size(nhbr2), 1, -1
          if (new_mesh%lface(2,n) == nhbr2(k2)) exit
        end do
        INSIST(k2 > 0)
        i2 = map2(k2)
      end associate
      select case (new_mesh%xlnode(n+1)-new_mesh%xlnode(n))
      case (6) ! prism gap cell
        n1 = 4; n2 = 5
      case (8) ! hex gap cell
        n1 = 5; n2 = 6
      case default
        INSIST(.false.)
      end select
      if (j1 <= new_mesh%ncell_onP) then
        ASSERT(ngbr_face(i1,j1) == 0)
        if (new_mesh%link_cell_id(n) > 0) then ! from a gap cell
          ngbr_face(i1,j1) = n1
        else ! do not want to do this; see NB1 above
          ngbr_face(i1,j1) = i2
        end if
      end if
      if (j2 <= new_mesh%ncell_onP) then
        ASSERT(ngbr_face(i2,j2) == 0)
        if (new_mesh%link_cell_id(n) > 0) then ! from a gap cell
          ngbr_face(i2,j2) = n2
        else ! do not want to do this; see NB1 above
          ngbr_face(i2,j2) = i1
        end if
      end if
    end do

    !! We can get away with giving no neighbor info for gap cells.
    ngbr_face(:,new_mesh%ncell_onP+1:) = DEGENERATE_FACE

  end subroutine init_ngbr_face

  !! Initialize the passed  NGBR_VRTX_ORIG and NGBR_VRTX indexing arrays.
  !! Note that prism gap cells are not mapped to degenerate hex cells expected
  !! by legacy mesh API clients.

  subroutine init_ngbr_vrtx (ngbr_vrtx_orig, ngbr_vrtx)

    use common_impl, only: OLD_TET_NODE_MAP, OLD_PYR_NODE_MAP, OLD_PRI_NODE_MAP

    integer, intent(out) :: ngbr_vrtx_orig(:,:), ngbr_vrtx(:,:)

    integer :: j, k, n

    ASSERT(size(ngbr_vrtx_orig,1) == 8)
    ASSERT(size(ngbr_vrtx_orig,2) == ncells)
    ASSERT(size(ngbr_vrtx,1) == 8)
    ASSERT(size(ngbr_vrtx,2) == ncells)

    !! Unpack the mesh cell node structure into NGBR_VRTX (real cells).
    ngbr_vrtx_orig = 0
    do j = 1, new_mesh%ncell_onP
      associate (cnode => new_mesh%cnode(new_mesh%xcnode(j):new_mesh%xcnode(j+1)-1))
        do k = 1, size(cnode)
          ngbr_vrtx(k,j) = cnode(k)
        end do
      end associate
    end do

    !! Unpack the mesh link node structure into NGBR_VRTX (gap cells).
    n = new_mesh%ncell_onP
    do j = 1, new_mesh%nlink_onP
      if (new_mesh%link_cell_id(j) == 0) cycle ! not from a gap cell
      n = n + 1
      associate (lnode => new_mesh%lnode(new_mesh%xlnode(j):new_mesh%xlnode(j+1)-1))
        do k = 1, size(lnode)
          ngbr_vrtx(k,n) = lnode(k)
        end do
      end associate
    end do

    !! Convert to legacy API degenerate hexes.
    do j = 1, ncells
      select case (mesh(j)%cell_shape)
      case (CELL_TET)
        ngbr_vrtx(:,j) = ngbr_vrtx(OLD_TET_NODE_MAP,j)
      case (CELL_PYRAMID)
        ngbr_vrtx(:,j) = ngbr_vrtx(OLD_PYR_NODE_MAP,j)
      case (CELL_PRISM)
        ngbr_vrtx(:,j) = ngbr_vrtx(OLD_PRI_NODE_MAP,j)
      end select
    end do

    !! Map to global node IDs.
    ngbr_vrtx_orig = new_mesh%node_ip%global_index(ngbr_vrtx)

    !! Convert off-process references to look-aside boundary array references.
    where (ngbr_vrtx > nnodes) ngbr_vrtx = nnodes - ngbr_vrtx

  end subroutine init_ngbr_vrtx

  !! Initialize the passed NGBR_CELLS_ALL and NGBR_CELLS_FACE structure arrays.
  !! The EE_GATHER_ALL_V_S_* procedures that operate with the returned arguments
  !! should use the cell index partition from the mesh object and NOT the legacy
  !! cell index partition LEGACY_CELL_IP from this module.  The reason is that
  !! these structure arrays ignore gap cells, and therefore have precisely the
  !  same cells as the mesh.  Neighbor lists will not include any neighbors that
  !! are gap cells, and while the arrays are sized to include neighbor lists for
  !! gap cells, these lists are entirely empty.
  !!
  !! This data is used only by the least squares operators, which I believe are
  !! only used by flow.  There is at least one regression problem that uses this
  !! (natural-conv-tet) but perhaps no more.  Bottom line is that this is very
  !! likely to break on any problem that uses the least squares operators and
  !! includes gap elements.  However, gap elements should not be included in
  !! flow anyway, so once things are converted over to the new mesh with gap
  !! elements separated, we should be okay.

  subroutine init_ngbr_cells_all (ngbr_cells_all, ngbr_cells_face)

    use common_impl, only: NEW_TET_SIDE_MAP,  NEW_PYR_SIDE_MAP, NEW_PRI_SIDE_MAP, NEW_HEX_SIDE_MAP
    use var_vector_module, only: int_var_vector, flatten

    type(int_var_vector), intent(out) :: ngbr_cells_all(:), ngbr_cells_face(:)

    integer :: j, k, n, b, new_side_map(6)
    integer, pointer :: cells(:)

    call init_ngbr_cells_all_aux (new_mesh, ngbr_cells_all)
    call init_ngbr_cells_face_aux (new_mesh, ngbr_cells_all, ngbr_cells_face)

    !! Map the face bitmasks to the legacy face ordering conventions.
    do j = 1, new_mesh%ncell_onP
      select case (new_mesh%xcnode(j+1)-new_mesh%xcnode(j))
      case (4)
        new_side_map = NEW_TET_SIDE_MAP
      case (5)
        new_side_map = NEW_PYR_SIDE_MAP
      case (6)
        new_side_map = NEW_PRI_SIDE_MAP
      case (8)
        new_side_map = NEW_HEX_SIDE_MAP
      case default
        INSIST(.false.)
      end select
      associate(fbits => ngbr_cells_face(j)%v)
        do n = 1, size(fbits)
          b = 0
          do k = 1, 6
            if (btest(fbits(n),pos=k)) b = ibset(b,pos=new_side_map(k))
          end do
          fbits(n) = b
        end do
      end associate
    end do

    !! Translate off-process references to boundary buffer references.
    cells => flatten(ngbr_cells_all)
    where (cells > ncells) cells = new_mesh%ncell_onP - cells
    INSIST(all(cells /= 0))

  end subroutine init_ngbr_cells_all

  !! This auxiliary subroutine initializes the NGBR_CELLS_ALL structure array
  !! that corresponds to the given new MESH object.  The array data is local
  !! cell IDs including off-process cells, and the subroutine requires that
  !! the the off-process cells in the mesh include all cells that share a node
  !! with an on-process cell.  The subroutine uses no module-scope variables.

  subroutine init_ngbr_cells_all_aux (mesh, ngbr_cells_all)

    use unstr_mesh_type
    use var_vector_module, only: int_var_vector, create
    use integer_set_type

    type(unstr_mesh), intent(in) :: mesh
    type(int_var_vector), intent(out) :: ngbr_cells_all(:)

    integer :: j, k, sizes(size(ngbr_cells_all))
    type(integer_set) :: nsets(mesh%nnode), csets(mesh%ncell_onP)

    ASSERT(size(ngbr_cells_all) >= mesh%ncell_onP)

    !! For each node, form the set of cells that contain it.
    do j = 1, mesh%ncell
      associate (cnode => mesh%cnode(mesh%xcnode(j):mesh%xcnode(j+1)-1))
        do k = 1, size(cnode)
          call nsets(cnode(k))%add(j)
        end do
      end associate
    end do

    !! For each cell, form the set of cells with which it shares a node.
    do j = 1, mesh%ncell_onP
      associate (cnode => mesh%cnode(mesh%xcnode(j):mesh%xcnode(j+1)-1))
        do k = 1, size(cnode)
          call csets(j)%add(nsets(cnode(k)))
        end do
      end associate
    end do

    !! Pack the data into an int_var_vector
    sizes(:mesh%ncell_onP) = csets%size()
    sizes(mesh%ncell_onP+1:) = 0
    call create (ngbr_cells_all, sizes)
    do j = 1, size(csets)
      call csets(j)%copy_to_array (ngbr_cells_all(j)%v)
    end do

  end subroutine init_ngbr_cells_all_aux

  !! This auxiliary subroutine initializes the NGBR_CELLS_FACE structure array
  !! that is the companion to the given NGBR_CELLS_ALL structure array.  The
  !! ordering of the face bits in the array data are relative to the new mesh.
  !! This subroutine uses no module scope data.

  subroutine init_ngbr_cells_face_aux (mesh, ngbr_cells_all, ngbr_cells_face)

    use unstr_mesh_type
    use cell_topology, only: cell_face_sig
    use var_vector_module, only: int_var_vector, create, sizes

    type(unstr_mesh), intent(in) :: mesh
    type(int_var_vector), intent(in)  :: ngbr_cells_all(:)
    type(int_var_vector), intent(out) :: ngbr_cells_face(:)

    integer :: i, j, k, ik, jk, n, m, b
    integer, pointer :: face_sig(:)

    ASSERT(size(ngbr_cells_face) == size(ngbr_cells_all))
    ASSERT(size(ngbr_cells_face) >= mesh%ncell_onP)

    call create (ngbr_cells_face, sizes(ngbr_cells_all))
    do j = 1, mesh%ncell_onP  ! for each on-process cell J
      associate (jnode => mesh%cnode(mesh%xcnode(j):mesh%xcnode(j+1)-1))
        face_sig => cell_face_sig(jnode)
        do n = 1, size(ngbr_cells_all(j)%v) ! for each neighbor cell I of J
          i = ngbr_cells_all(j)%v(n)
          m = 0 ! bit mask tagging J nodes contained in I
          associate (inode => mesh%cnode(mesh%xcnode(i):mesh%xcnode(i+1)-1))
            do jk = 1, size(jnode)  ! for each node of cell J
              do ik = 1, size(inode)  ! for
                if (inode(ik) == jnode(jk)) then
                  m = ibset(m, pos=jk-1)
                  cycle
                end if
              end do
            end do
          end associate
          b = 0 ! bit mask tagging J faces touching cell I
          do k = 1, size(face_sig)
            if (iand(m, face_sig(k)) /= 0) b = ibset(b, pos=k)
          end do
          ngbr_cells_face(j)%v(n) = b
        end do
      end associate
    end do

  end subroutine init_ngbr_cells_face_aux

end module mesh_impl
