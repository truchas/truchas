!!
!! MESH_IMPL
!!
!! Implements the legacy MESH structure array and related data and procedures.
!!

#include "f90_assert.fpp"

module mesh_impl

  use common_impl, only: ncells, nnodes, new_mesh
  use var_vector_types, only: int_var_vector
  use index_partitioning, only: ip_desc
  implicit none
  private

  public :: init_mesh_impl
  public :: is_face_ngbr, mesh_collate_vertex

  integer, parameter, public :: DEGENERATE_FACE = -huge(1)

  type(ip_desc), public :: cell_ip

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
    call init_ngbr_cell (ngbr_cell_orig, ngbr_cell, cell_ip)
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
    call init_ngbr_vrtx_mapped (ngbr_vrtx_orig, ngbr_vrtx)
    do j = 1, ncells
      mesh(j)%ngbr_vrtx      = ngbr_vrtx(:,j)
      mesh(j)%ngbr_vrtx_orig = ngbr_vrtx_orig(:,j)
    end do
    deallocate(ngbr_vrtx_orig, ngbr_vrtx)

    !! MESH%NGBR_CELLS_ALL, MESH%NGBR_CELLS_FACE
    call init_ngbr_cells_all_mapped (mesh%ngbr_cells_all, mesh%ngbr_cells_face)

    !! MESH%CBLOCKID
    call init_cblockid (mesh%cblockid)

  end subroutine init_mesh

  subroutine init_cell_shape (cell_shape)

    use common_impl, only: pcell_old_to_new
    use parallel_permutations, only: rearrange

    integer, intent(out) :: cell_shape(:)

    integer :: j, src(new_mesh%ncell_onP)

    ASSERT(size(cell_shape) == ncells)

    do j = 1, size(src)
      select case (new_mesh%xcnode(j+1)-new_mesh%xcnode(j))
      case (4)
        src(j) = CELL_TET
      case (5)
        src(j) = CELL_PYRAMID
      case (6)
        src(j) = CELL_PRISM
      case (8)
        src(j) = CELL_HEX
      case default
        INSIST(.false.)
      end select
    end do

    !! Map the cell shape array to the old mesh.  The default is assigned to
    !! gap cells, and appears to be correct for the hex gap cells that Truchas
    !! generates.  It is not correct for prism gap cells, but no regression
    !! test currently includes these.
    call rearrange (pcell_old_to_new, cell_shape, src, default=GAP_ELEMENT_5)

  end subroutine init_cell_shape

  subroutine init_cblockid (cblockid)

    use common_impl, only: new_mesh, pcell_old_to_new, pgap_old_to_new, gap_link_mask, gap_cells
    use parallel_permutations, only: rearrange
    use bitfield_type, only: popcnt, trailz

    integer, intent(out) :: cblockid(:)

    integer :: j
    integer, allocatable :: src(:), dest(:)

    ASSERT(size(cblockid) == ncells)

    !! Cell set IDs
    allocate(src(new_mesh%ncell_onP))
    do j = 1, new_mesh%ncell_onP
      associate (bitmask => new_mesh%cell_set_mask(j))
        INSIST(popcnt(bitmask) == 1)
        src(j) = new_mesh%cell_set_id(trailz(bitmask))
      end associate
    end do
    call rearrange (pcell_old_to_new, cblockid, src, default = 0)
    deallocate(src)
    ASSERT(all(cblockid(gap_cells) == 0))

    !! Link set IDs for links coming from gap cells
    allocate(src(new_mesh%nlink_onP),dest(size(gap_cells)))
    do j = 1, new_mesh%nlink_onP
      associate (bitmask => new_mesh%link_set_mask(j))
        INSIST(popcnt(bitmask) == 1)
        src(j) = new_mesh%link_set_id(trailz(bitmask))
      end associate
    end do
    src = pack(src, mask=gap_link_mask)
    call rearrange (pgap_old_to_new, dest, src)
    cblockid(gap_cells) = dest

    INSIST(all(cblockid > 0))

  end subroutine init_cblockid

  subroutine init_ngbr_cell (ngbr_cell_orig, ngbr_cell, cell_ip)

    use common_impl, only: OLD_TET_SIDE_MAP, OLD_PYR_SIDE_MAP, OLD_PRI_SIDE_MAP, OLD_HEX_SIDE_MAP
    use common_impl, only: pcell_old_to_new, pcell_new_to_old
    use common_impl, only: gap_cells, gap_link_mask, pgap_new_to_old
    use parallel_permutations, only: rearrange
    use index_partitioning, only: localize_index_array, gather_boundary

    integer, intent(out) :: ngbr_cell_orig(:,:), ngbr_cell(:,:)
    type(ip_desc), intent(out) :: cell_ip

    integer :: i, j, jj, k, n
    integer :: old_gid(ncells), new_gid(new_mesh%ncell), link_gid(new_mesh%nlink)
    integer :: src(6,new_mesh%ncell_onP)
    integer, allocatable :: offP_index(:), src_buf(:), dest_buf(:)
    logical :: mask(new_mesh%nlink)

    !! Push the array of old-mesh global cell IDs to the new mesh.
    call cell_ip%init (ncells)
    do j = 1, ncells
      old_gid(j) = cell_ip%global_index(j)
    end do
    call rearrange (pcell_new_to_old, new_gid(:new_mesh%ncell_onP), old_gid) ! drops gap cells
    call gather_boundary (new_mesh%cell_ip, new_gid)

    !! Push the gap cell global IDs to the new mesh.
    src_buf = old_gid(gap_cells)
    allocate(dest_buf(count(gap_link_mask)))
    call rearrange (pgap_new_to_old, dest_buf, src_buf)
    link_gid(:new_mesh%nlink_onP) = unpack(dest_buf, mask=gap_link_mask, field=0)
    call gather_boundary (new_mesh%link_ip, link_gid)

    !! Generate the cell neighbor array (old mesh global IDs)
    do j = 1, new_mesh%ncell_onP
      associate (list => new_mesh%cnhbr(new_mesh%xcnhbr(j):new_mesh%xcnhbr(j+1)-1))
        src(:,j) = DEGENERATE_FACE
        do k = 1, size(list)
          if (list(k) > 0) then
            src(k,j) = new_gid(list(k))
          else
            src(k,j) = 0
          end if
        end do
      end associate
    end do

    !! Fill in neighbor data from gap cells.
    mask(:new_mesh%nlink_onP) = gap_link_mask
    call gather_boundary (new_mesh%link_ip, mask)
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
        INSIST(src(k,j) == 0)
        if (mask(n)) then ! neighbor value is gap cell value
          src(k,j) = link_gid(n)
        else ! neighbor value is other link neighbor value
          jj = new_mesh%lnhbr(1+modulo(i,2),n)
          src(k,j) = new_gid(jj)
        end if
      end do
    end do

    !! Convert to the old mesh convention for ordering neighbors.
    do j = 1, new_mesh%ncell_onP
      select case (new_mesh%xcnode(j+1)-new_mesh%xcnode(j))
      case (4)
        src(:,j) = src(OLD_TET_SIDE_MAP,j)
      case (5)
        src(:,j) = src(OLD_PYR_SIDE_MAP,j)
      case (6)
        src(:,j) = src(OLD_PRI_SIDE_MAP,j)
      case (8)
        src(:,j) = src(OLD_HEX_SIDE_MAP,j)
      end select
    end do

    !! Map the neighbor array back to the old mesh.  We have no data
    !! for gap cells; set them to 0 and see if we get away with it.
    call rearrange (pcell_old_to_new, ngbr_cell_orig, src, default=0)

    !! Localize the ngbr_cell_orig array.
    ngbr_cell = ngbr_cell_orig
    where (ngbr_cell < 0) ngbr_cell = 0 ! temp convert degen side tags to ignored bndry sides
    call localize_index_array (cell_ip, ngbr_cell, offP_index)
    where (ngbr_cell_orig < 0) ngbr_cell = ngbr_cell_orig ! restore degen side tags

    call cell_ip%add_offP_index (offP_index)
    deallocate(offP_index)

    !! Convert off-process references to look-aside boundary array references.
    where (ngbr_cell > ncells) ngbr_cell = ncells - ngbr_cell

  end subroutine init_ngbr_cell

  !! Initialize the passed NGBR_FACE array.  The handling of gap cells is
  !! incomplete and fragile.  We do not return any valid neighbor face data
  !! for gap cells; we use DEGENERATE_FACE for all.  Gap cells recovered
  !! from links may not have the same orientation as the original gap cells,
  !! leading to incorrect face values for neighbors that are gap cells. This
  !! is definitely true for prism gap cells, but for hex gap cells we seem
  !! to have gotten lucky; those face values are good. This issue becomes
  !! moot once we no longer try to compare against the original data.

  subroutine init_ngbr_face (ngbr_face)

    use common_impl, only: NEW_TET_SIDE_MAP, NEW_PYR_SIDE_MAP, NEW_PRI_SIDE_MAP, NEW_HEX_SIDE_MAP
    use common_impl, only: pcell_old_to_new
    use parallel_permutations, only: rearrange

    integer, intent(out) :: ngbr_face(:,:)

    integer :: i1, i2, j1, j2, k1, k2, n, n1, n2, map1(6), map2(6)
    integer, allocatable :: src(:,:)

    ASSERT(size(ngbr_face,1) == 6)
    ASSERT(size(ngbr_face,2) == ncells)

    allocate(src(6,new_mesh%ncell_onP))
    src = DEGENERATE_FACE
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
          if (src(i1,j1) > 0) cycle ! already defined
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
              src(i1,j1) = i2
              if (j2 <= size(src,2)) src(i2,j2) = i1
            end associate
          else  ! boundary face
            src(i1,j1) = 0
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
      if (j1 <= size(src,2)) then
        ASSERT(src(i1,j1) == 0)
        if (new_mesh%link_cell_id(n) > 0) then ! from a gap cell
          src(i1,j1) = n1
        else
          src(i1,j1) = i2
        end if
      end if
      if (j2 <= size(src,2)) then
        ASSERT(src(i2,j2) == 0)
        if (new_mesh%link_cell_id(n) > 0) then ! from a gap cell
          src(i2,j2) = n2
        else
          src(i2,j2) = i1
        end if
      end if
    end do

    call rearrange (pcell_old_to_new, ngbr_face, src, default=DEGENERATE_FACE)

  end subroutine init_ngbr_face

  !! Initialize the passed  NGBR_VRTX_ORIG and NGBR_VRTX indexing arrays.
  !! Note that prism gap cells are not mapped to degenerate hex cells expected
  !! by legacy mesh API clients.

  subroutine init_ngbr_vrtx_mapped (ngbr_vrtx_orig, ngbr_vrtx)

    use common_impl, only: OLD_TET_NODE_MAP, OLD_PYR_NODE_MAP, OLD_PRI_NODE_MAP
    use common_impl, only: pcell_old_to_new
    use common_impl, only: pgap_old_to_new, gap_cells, gap_link_mask
    use parallel_permutations, only: rearrange
    use index_partitioning, only: localize_index_array
    use truchas_logging_services

    integer, intent(out) :: ngbr_vrtx_orig(:,:), ngbr_vrtx(:,:)

    integer :: j, k
    integer, allocatable :: src(:,:), src_buf(:), dest_buf(:), offP_index(:)

    ASSERT(size(ngbr_vrtx_orig,1) == 8)
    ASSERT(size(ngbr_vrtx_orig,2) == ncells)
    ASSERT(size(ngbr_vrtx,1) == 8)
    ASSERT(size(ngbr_vrtx,2) == ncells)

    !! Generate the node neighbor array (old mesh global IDs)
    allocate(src(8,new_mesh%ncell_onP))
    src = 0
    do j = 1, new_mesh%ncell_onP
      associate (cnode => new_mesh%cnode(new_mesh%xcnode(j):new_mesh%xcnode(j+1)-1))
        do k = 1, size(cnode)
          src(k,j) = new_mesh%node_ip%global_index(cnode(k))
        end do
      end associate
    end do

    !! Move the array back to the old mesh; no data for gap cells.
    call rearrange (pcell_old_to_new, ngbr_vrtx_orig, src, default=0)
    deallocate(src)

    !! Generate the node neighbor array for links (old mesh global IDs)
    allocate(src(8,new_mesh%nlink_onP))
    src = 0
    do j = 1, new_mesh%nlink_onP
      associate (lnode => new_mesh%lnode(new_mesh%xlnode(j):new_mesh%xlnode(j+1)-1))
        do k = 1, size(lnode)
          src(k,j) = new_mesh%node_ip%global_index(lnode(k))
        end do
      end associate
    end do

    !! Move the array back to the old mesh, filling in the data for gap cells.
    allocate(dest_buf(size(gap_cells)))
    do k = 1, size(src,1)
      src_buf = pack(src(k,:), mask=gap_link_mask)
      call rearrange (pgap_old_to_new, dest_buf, src_buf)
      ngbr_vrtx_orig(k,gap_cells) = dest_buf
    end do

    !! Convert to the legacy mesh node ordering convention.
    do j = 1, ncells
      select case (mesh(j)%cell_shape)
      case (CELL_TET)
        ngbr_vrtx_orig(:,j) = ngbr_vrtx_orig(OLD_TET_NODE_MAP,j)
      case (CELL_PYRAMID)
        ngbr_vrtx_orig(:,j) = ngbr_vrtx_orig(OLD_PYR_NODE_MAP,j)
      case (CELL_PRISM)
        ngbr_vrtx_orig(:,j) = ngbr_vrtx_orig(OLD_PRI_NODE_MAP,j)
      end select
    end do

    !! Copy to the NGBR_VRTX array and localize it.
    ngbr_vrtx = ngbr_vrtx_orig
    call localize_index_array (new_mesh%node_ip, ngbr_vrtx, offP_index)
    INSIST(size(offP_index) == 0)

    !! Convert off-process references to look-aside boundary array references.
    where (ngbr_vrtx > nnodes) ngbr_vrtx = nnodes - ngbr_vrtx

  end subroutine init_ngbr_vrtx_mapped

  !! Initialize the passed NGBR_CELLS_ALL and NGBR_CELLS_FACE structure arrays,
  !! and return the associated cell index partition object CELL_IP which must be
  !! used by the EE_GATHER_ALL_V_S_* procedures that operate with the returned
  !! arguments. CELL_IP is the analog of EE_ALL_NGBR_TRACE from gs_info_module.
  !! This data is all with respect to the legacy mesh ordering and partitioning
  !! of cells and its convention for face numbering.
  !!
  !! This implementation does not fill in data for gap elements.  That means
  !! the neighbor lists will not include any neighbors that are gap elements,
  !! and the neighbor lists for gap elements will be entirely empty, not even
  !! including themselves.
  !!
  !! This data is used only by the least squares operators, which I believe are
  !! only used by flow.  There is at least one regression problem that uses this
  !! (natural-conv-tet) but perhaps no more.  Bottom line is that this is very
  !! likely to break on any problem that uses the least squares operators and
  !! includes gap elements.  However, gap elements should not be included in
  !! flow anyway, so once things are converted over to the new mesh with gap
  !! elements separated, we should be okay.

  subroutine init_ngbr_cells_all_mapped (ngbr_cells_all, ngbr_cells_face)

    use common_impl, only: pcell_old_to_new, pcell_new_to_old
    use common_impl, only: NEW_TET_SIDE_MAP,  NEW_PYR_SIDE_MAP, NEW_PRI_SIDE_MAP, NEW_HEX_SIDE_MAP
    use parallel_permutations, only: rearrange
    use permutations, only: reorder
    use sort_module, only: heapsort
    use index_partitioning, only: ip_desc, gather_boundary, localize_index_struct
    use var_vector_module, only: int_var_vector, create, sizes, flatten
    use parallel_communication, only: global_maxval

    type(int_var_vector), intent(out) :: ngbr_cells_all(:), ngbr_cells_face(:)

    integer :: j, k, n, b
    integer :: old_gid(ncells), new_gid(new_mesh%ncell), new_side_map(6)
    integer, allocatable :: perm(:), new_sizes(:), old_sizes(:), src(:,:), dest(:,:), offP_index(:)
    type(int_var_vector), allocatable :: new_ngbr_cells_all(:), new_ngbr_cells_face(:)
    integer, pointer :: cells(:)

    ASSERT(size(ngbr_cells_all) == ncells)
    ASSERT(size(ngbr_cells_face) == ncells)

    !! Push the array of old-mesh global cell IDs to the new mesh.
    do j = 1, ncells
      old_gid(j) = new_mesh%cell_ip%global_index(j)
    end do
    call rearrange (pcell_new_to_old, new_gid(:new_mesh%ncell_onP), old_gid) ! drops gap cells
    call gather_boundary (new_mesh%cell_ip, new_gid)

    !! Create the NGBR_CELLS_ALL data relative to the new mesh ordering.
    allocate(new_ngbr_cells_all(new_mesh%ncell_onP))
    call init_ngbr_cells_all_aux (new_mesh, new_ngbr_cells_all)

    !! Create the companion NGBR_CELLS_FACE data relative to new mesh conventions.
    allocate(new_ngbr_cells_face(new_mesh%ncell_onP))
    call init_ngbr_cells_face_aux (new_mesh, new_ngbr_cells_all, new_ngbr_cells_face)

    !! Map the face bitmasks to the legacy face ordering conventions.
    do j = 1, size(new_ngbr_cells_face)
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
      associate(fbits => new_ngbr_cells_face(j)%v)
        do n = 1, size(fbits)
          b = 0
          do k = 1, 6
            if (btest(fbits(n),pos=k)) b = ibset(b,pos=new_side_map(k))
          end do
          fbits(n) = b
        end do
      end associate
    end do

    !! Map the data to the cell ordering used by the legacy mesh.
    !! We also re-sort the mapped lists to recover that behavior.
    allocate(perm(maxval(sizes(new_ngbr_cells_all))))
    do j = 1, size(new_ngbr_cells_all)
      associate (jngbr => new_ngbr_cells_all(j)%v, jface => new_ngbr_cells_face(j)%v)
        jngbr = new_gid(jngbr)
        call heapsort (jngbr, perm(:size(jngbr)))
        call reorder  (jngbr, perm(:size(jngbr)))
        call reorder  (jface, perm(:size(jface)))
      end associate
    end do
    deallocate(perm)

    !! Move this data to the old partition.  This is awkward because there is
    !! a variable amount of data per cell and no easy way to do it with PGSLib.
    !! So we unpack into a sufficiently large rank-2 array and move it instead.
    new_sizes = sizes(new_ngbr_cells_all)
    allocate(old_sizes(ncells))
    call rearrange (pcell_old_to_new, old_sizes, new_sizes, default=0)
    n = global_maxval(new_sizes)
    allocate(src(n,size(new_sizes)), dest(n,ncells))

    !! Push the NGBR_CELL_ALL data to the old mesh.
    src = 0
    do j = 1, size(new_ngbr_cells_all)
      n = size(new_ngbr_cells_all(j)%v)
      src(1:n,j) = new_ngbr_cells_all(j)%v
    end do
    call rearrange (pcell_old_to_new, dest, src, default=0)

    !! Pack the data into an int_var_vector.
    call create (ngbr_cells_all, old_sizes)
    do j = 1, ncells
      ngbr_cells_all(j)%v = dest(:old_sizes(j),j)
    end do

    !! Push the NGBR_CELL_FACE data to the old mesh.
    src = 0
    do j = 1, size(new_ngbr_cells_face)
      n = size(new_ngbr_cells_face(j)%v)
      src(1:n,j) = new_ngbr_cells_face(j)%v
    end do
    call rearrange (pcell_old_to_new, dest, src, default=0)

    !! Pack the data into an int_var_vector.
    call create (ngbr_cells_face, old_sizes)
    do j = 1, ncells
      ngbr_cells_face(j)%v = dest(:old_sizes(j),j)
    end do
    deallocate(src, dest)

    !! Create an index partition
    cells => flatten(ngbr_cells_all)
    call localize_index_struct (new_mesh%cell_ip, old_sizes, cells, offP_index)
    INSIST(size(offP_index)==0)

    !! Translate off-process references to boundary buffer references.
    where (cells > ncells) cells = ncells - cells
    INSIST(all(cells /= 0))

  end subroutine init_ngbr_cells_all_mapped

  !! Initialize the passed NGBR_CELLS_ALL and NGBR_CELLS_FACE structure arrays,
  !! that corresponds to the given new MESH object.  These arrays work with the
  !! cell index partition component MESH%CELL_IP.  This is what will be used
  !! once we switch to using the new mesh partitioning and conventions.
  !! This is untested and hence commented out.

  !subroutine init_ngbr_cells_all (mesh, ngbr_cells_all, ngbr_cells_face)
  !
  !  type(unstr_mesh), intent(in) :: mesh
  !  type(int_var_vector), intent(out) :: ngbr_cells_all(:), ngbr_cells_face(:)
  !
  !  ASSERT(size(ngbr_cells_all) == mesh%ncell_onP)
  !  ASSERT(size(ngbr_cells_face) == mesh%ncell_onP)
  !
  !  call init_ngbr_cells_all_aux (mesh, ngbr_cells_all)
  !  call init_ngbr_cells_face_aux (mesh, ngbr_cells_all, ngbr_cells_face)
  !
  !end subroutine init_ngbr_cells_all

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

    integer :: j, k
    type(integer_set) :: nsets(mesh%nnode), csets(mesh%ncell_onP)

    ASSERT(size(ngbr_cells_all) == mesh%ncell_onP)

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
#ifdef INTEL_INTEGER_SET_ICE
          call csets(j)%add_set(nsets(cnode(k)))
#else
          call csets(j)%add(nsets(cnode(k)))
#endif
        end do
      end associate
    end do

    !! Pack the data into an int_var_vector
    call create (ngbr_cells_all, csets%size())
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
    ASSERT(size(ngbr_cells_face) <= mesh%ncell)

    call create (ngbr_cells_face, sizes(ngbr_cells_all))
    do j = 1, size(ngbr_cells_all)  ! for each on-process cell J
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

  !! Copy of function from MESH_MODULE
  function mesh_collate_vertex (mesh)
    use common_impl, only: ncells_tot
    use parallel_communication, only: is_IOP, collate
    type(mesh_connectivity), intent(in) :: mesh(:)
    integer, pointer :: mesh_collate_vertex(:,:)
    integer :: k
    allocate(mesh_collate_vertex(8,merge(ncells_tot,0,is_IOP)))
    do k = 1, 8
      call collate (mesh_collate_vertex(k,:), mesh%ngbr_vrtx_orig(k))
    end do
  end function mesh_collate_vertex

end module mesh_impl
