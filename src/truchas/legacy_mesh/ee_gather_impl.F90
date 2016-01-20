!!
!! EE_GATHER_IMPL
!!
!! Implements the legacy element-element gather procedures.
!!

#include "f90_assert.fpp"

module ee_gather_impl

  use kinds, only: r8
  use common_impl, only: ncells, new_mesh
  implicit none
  private

  public :: ee_gather, gather_boundarydata, init_ee_gather_impl

  interface ee_gather
    !! Non-validating versions
    !module procedure ee_gather_int32
    !module procedure ee_gather_real64
    !module procedure ee_gather_log
    !module procedure ss_gather_real64
    !module procedure ss_gather_log
    !module procedure ee_gather_all_v_s_real64

    !! Validating versions
    module procedure ee_gather_validate_int32
    module procedure ee_gather_validate_real64
    module procedure ee_gather_validate_log
    module procedure ss_gather_validate_real64
    module procedure ss_gather_validate_log
    module procedure ee_gather_all_v_s_validating
  end interface

  interface gather_boundarydata
    module procedure gather_boundarydata_real64
  end interface

  !! Masks for validating EE_GATHER results; .TRUE. marks elements to validate.
  logical, allocatable :: ee_test_mask(:,:), ee_test_mask_log(:,:)

  !! Face categories for defining EE_GATHER test masks.
  integer, parameter :: DEGEN_SIDE  = 1 ! degenerate cell sides
  integer, parameter :: GAP_SIDE_A  = 2 ! sides adjacent to interface from gap elements
  integer, parameter :: GAP_SIDE_B  = 4 ! sides adjacent to interface from side sets
  integer, parameter :: GAP_CELL    = 8 ! sides of gap elements
  integer, parameter :: GAP_SIDE    = ior(GAP_SIDE_A, GAP_SIDE_B)

contains

  subroutine init_ee_gather_impl

    call make_ee_test_mask (new_mesh, ior(GAP_SIDE,GAP_CELL), ee_test_mask)
    call make_ee_test_mask (new_mesh, ior(GAP_SIDE_A,GAP_CELL), ee_test_mask_log)

  contains

    subroutine make_ee_test_mask (mesh, ignore, ee_test_mask)

      use common_impl, only: pcell_old_to_new
      use common_impl, only: NEW_TET_SIDE_MAP, NEW_PRI_SIDE_MAP, NEW_HEX_SIDE_MAP, NEW_PYR_SIDE_MAP
      use parallel_permutations, only: rearrange
      use unstr_mesh_type

      type(unstr_mesh), intent(in) :: mesh
      integer, intent(in) :: ignore
      logical, allocatable, intent(out) :: ee_test_mask(:,:)

      integer :: i, j, k, n
      integer, allocatable :: bitmask1(:), bitmask2(:)
      logical :: ignore_degen_side, ignore_gap_side_a, ignore_gap_side_b, ignore_gap_cell

      ignore_degen_side = (iand(ignore,DEGEN_SIDE) /= 0)
      ignore_gap_side_a = (iand(ignore,GAP_SIDE_A) /= 0)
      ignore_gap_side_b = (iand(ignore,GAP_SIDE_B) /= 0)
      ignore_gap_cell   = (iand(ignore,GAP_CELL)   /= 0)

      !! Start with ignoring nothing and then start clearing bits.
      allocate(bitmask1(mesh%ncell_onP))
      bitmask1 = 126 ! set bits 1 through 6

      !! Degenerate faces.
      if (ignore_degen_side) then
        do j = 1, mesh%ncell_onP
          select case (mesh%xcnode(j+1)-mesh%xcnode(j))
          case (4)    ! tet
            bitmask1(j) = ibclr(bitmask1(j),pos=2)
            bitmask1(j) = ibclr(bitmask1(j),pos=6)
          case (5:6)  ! pyramid, prism
            bitmask1(j) = ibclr(bitmask1(j),pos=6)
          end select
        end do
      end if

      !! Gap faces.
      do n = 1, mesh%nlink
        if ((mesh%link_cell_id(n) > 0  .and. ignore_gap_side_a) .or. &
            (mesh%link_cell_id(n) == 0 .and. ignore_gap_side_b)) then
          do i = 1, 2
            j = mesh%lnhbr(i,n)
            if (j > mesh%ncell_onP) cycle
            associate (cface => mesh%cface(mesh%xcface(j):mesh%xcface(j+1)-1))
              do k = size(cface), 1, -1
                if (mesh%lface(i,n) == cface(k)) exit
              end do
              INSIST(k > 0)
            end associate
            select case (mesh%xcnode(j+1)-mesh%xcnode(j))
            case (4)
              k = NEW_TET_SIDE_MAP(k)
            case (5)
              k = NEW_PYR_SIDE_MAP(k)
            case (6)
              k = NEW_PRI_SIDE_MAP(k)
            case (8)
              k = NEW_HEX_SIDE_MAP(k)
            end select
            bitmask1(j) = ibclr(bitmask1(j),pos=k)
          end do
        end if
      end do

      !! Map to the old mesh, filling in data for gap cells.
      allocate(bitmask2(ncells))
      if (ignore_gap_cell) then
        call rearrange (pcell_old_to_new, bitmask2, bitmask1, default=0)
      else
        call rearrange (pcell_old_to_new, bitmask2, bitmask1, default=126)
      end if

      !! Convert bit mask to the final logical array.
      allocate(ee_test_mask(6,ncells))
      do j = 1, ncells
        do k = 1, 6
          ee_test_mask(k,j) = btest(bitmask2(j),pos=k)
        end do
      end do

    end subroutine make_ee_test_mask

  end subroutine init_ee_gather_impl


  subroutine ee_gather_real64 (dest, src, boundary)

    use mesh_impl, only: old_mesh => mesh, cell_ip, DEGENERATE_FACE
    use index_partitioning, only: gather_boundary

    real(r8), intent(inout) :: dest(:,:)
    real(r8), intent(in) :: src(:)
    real(r8), pointer, optional :: boundary(:)

    integer :: j, k, n
    real(r8), pointer :: buffer(:)

    if (present(boundary)) then
      if (.not.associated(boundary)) then
        allocate(boundary(cell_ip%offP_size()))
        call gather_boundary (cell_ip, src, boundary)
      end if
      buffer => boundary
    else
      allocate(buffer(cell_ip%offP_size()))
      call gather_boundary (cell_ip, src, buffer)
    end if

    do j = 1, ncells
      do k = 1, 6
        n = old_mesh(j)%ngbr_cell(k)
        select case (n)
        case (1:)
          dest(k,j) = src(n)
        case (DEGENERATE_FACE+1:-1)
          dest(k,j) = buffer(-n)
        case default
          dest(k,j) = 0.0_r8
        end select
      end do
    end do

    if (.not.present(boundary)) deallocate(buffer)

  end subroutine ee_gather_real64


  subroutine ee_gather_int32 (dest, src)

    use common_impl, only: ncells
    use mesh_impl, only: mesh, cell_ip, DEGENERATE_FACE
    use index_partitioning, only: gather_boundary

    integer, intent(inout) :: dest(:,:)
    integer, intent(in) :: src(:)

    integer :: j, k, n
    integer :: buffer(cell_ip%offP_size())

    call gather_boundary (cell_ip, src, buffer)

    do j = 1, ncells
      do k = 1, 6
        n = mesh(j)%ngbr_cell(k)
        select case (n)
        case (1:)
          dest(k,j) = src(n)
        case (DEGENERATE_FACE+1:-1)
          dest(k,j) = buffer(-n)
        case default
          dest(k,j) = 0
        end select
      end do
    end do

  end subroutine ee_gather_int32


  subroutine ee_gather_log (dest, src)

    use common_impl, only: ncells
    use mesh_impl, only: mesh, cell_ip, DEGENERATE_FACE
    use index_partitioning, only: gather_boundary

    logical, intent(inout) :: dest(:,:)
    logical, intent(in) :: src(:)

    integer :: j, k, n
    logical :: buffer(cell_ip%offP_size())

    call gather_boundary (cell_ip, src, buffer)

    do j = 1, ncells
      do k = 1, 6
        n = mesh(j)%ngbr_cell(k)
        select case (n)
        case (1:)
          dest(k,j) = src(n)
        case (DEGENERATE_FACE+1:-1)
          dest(k,j) = buffer(-n)
        case default
          dest(k,j) = .false.
        end select
      end do
    end do

  end subroutine ee_gather_log


  subroutine ss_gather_real64 (dest, src)

    use common_impl, only: ncells
    use mesh_impl, only: mesh, cell_ip, DEGENERATE_FACE
    use index_partitioning, only: gather_boundary

    real(r8), intent(inout) :: dest(:,:)
    real(r8), intent(in) :: src(:,:)

    integer :: j, jj, k, kk
    real(r8) :: buffer(6,cell_ip%offP_size())

    ASSERT(size(dest,1) == 6)
    ASSERT(size(dest,2) == ncells)
    ASSERT(size(src,1) == 6)
    ASSERT(size(src,2) == ncells)

    call gather_boundary (cell_ip, src, buffer)

    do j = 1, ncells
      do k = 1, 6
        jj = mesh(j)%ngbr_cell(k)
        kk = mesh(j)%ngbr_face(k)
        select case (jj)
        case (1:) ! on-process neighbor
          dest(k,j) = src(kk,jj)
        case (DEGENERATE_FACE+1:-1) ! off-process neighbor
          dest(k,j) = buffer(kk,-jj)
        case default  ! boundary face or degenerate face
          dest(k,j) = 0.0_r8
        end select
      end do
    end do

  end subroutine ss_gather_real64


  subroutine ss_gather_log (dest, src)

    use common_impl, only: ncells
    use mesh_impl, only: mesh, cell_ip, DEGENERATE_FACE
    use index_partitioning, only: gather_boundary

    logical, intent(inout) :: dest(:,:)
    logical, intent(in) :: src(:,:)

    integer :: j, jj, k, kk
    logical :: buffer(6,cell_ip%offP_size())

    ASSERT(size(dest,1) == 6)
    ASSERT(size(dest,2) == ncells)
    ASSERT(size(src,1) == 6)
    ASSERT(size(src,2) == ncells)

    call gather_boundary (cell_ip, src, buffer)

    do j = 1, ncells
      do k = 1, 6
        jj = mesh(j)%ngbr_cell(k)
        kk = mesh(j)%ngbr_face(k)
        select case (jj)
        case (1:) ! on-process neighbor
          dest(k,j) = src(kk,jj)
        case (DEGENERATE_FACE+1:-1) ! off-process neighbor
          dest(k,j) = buffer(kk,-jj)
        case default  ! boundary face or degenerate face
          dest(k,j) = .false.
        end select
      end do
    end do

  end subroutine ss_gather_log


  subroutine ee_gather_all_v_s_real64 (dest, source, boundary, range)

    use var_vector_types, only: real_var_vector
    use mesh_impl, only: mesh, cell_ip_all
    use index_partitioning, only: gather_boundary

    real(r8) :: source(:)
    type(real_var_vector) :: dest(:)
    real(r8), pointer, optional :: boundary(:)
    integer, optional :: range(:)

    integer :: j, k
    real(r8), pointer :: buffer(:) => null()

    if (present(boundary)) buffer => boundary
    if (.not.associated(buffer)) then
      allocate(buffer(cell_ip_all%offP_size()))
      call gather_boundary (cell_ip_all, source, buffer)
      if (present(boundary)) boundary => buffer
    end if

    ASSERT(size(source) == cell_ip_all%onP_size())
    ASSERT(size(buffer) == cell_ip_all%offP_size())
    ASSERT(size(dest) <= size(mesh))

    !! Will assume this actual usage when RANGE is specified.
    if (present(range)) then
      INSIST(size(range) == 2)
      INSIST(range(1) == 1)
      INSIST(range(2) == size(dest))
    end if

    do j = 1, size(dest)
      associate (jngbr => mesh(j)%ngbr_cells_all%v, jdest => dest(j)%v)
        ASSERT(size(jdest) == size(jngbr))
        do k = 1, size(jngbr)
          if (jngbr(k) > 0) then
            jdest(k) = source(jngbr(k))
          else
            jdest(k) = buffer(-jngbr(k))
          end if
        end do
      end associate
    end do

    if (.not.present(boundary)) deallocate(buffer)

  end subroutine ee_gather_all_v_s_real64


  subroutine gather_boundarydata_real64 (boundary, source)
    use mesh_impl, only: cell_ip
    use index_partitioning, only: gather_boundary
    real(r8), pointer :: boundary(:)
    real(r8), intent(in) :: source(:)
    ASSERT(size(source) == cell_ip%onP_size())
    if (associated(boundary)) return  ! signal to do nothing
    allocate(boundary(cell_ip%offP_size()))
    call gather_boundary (cell_ip, source, boundary)
  end subroutine gather_boundarydata_real64

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! TESTING CODE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ee_gather_validate_real64 (dest, src, boundary)

    use truchas_logging_services
    use parallel_communication, only: global_any, global_all
    use gs_module, only: ee_gather_old => ee_gather
    use mesh_module, only: old_mesh => mesh

    real(r8), intent(inout) :: dest(:,:)
    real(r8), intent(in) :: src(:)
    real(r8), pointer, optional :: boundary(:)

    integer :: j, unit
    real(r8), allocatable :: dest_ref(:,:)
    real(r8), pointer :: boundary_ref(:)

    dest_ref = dest ! dest is intent(inout); capture it before it is modified

    if (present(boundary)) then
      boundary_ref => boundary
      call ee_gather_old (dest_ref, src, boundary=boundary_ref)
      call ee_gather_real64 (dest, src, boundary)
      if (.not.associated(boundary_ref,boundary)) deallocate(boundary_ref)
    else
      call ee_gather_old (dest_ref, src)
      call ee_gather_real64 (dest, src)
    end if

!#define DEBUG
#ifdef DEBUG
#undef DEBUG
    if (global_any(dest /= dest_ref .and. ee_test_mask)) then
      unit = TLS_debug_unit()
      call TLS_debug ('')
      call TLS_debug ('EE_GATHER_REAL64 ***********************************************')
      do j = 1, size(dest,2)
        if (any(dest(:,j) /= dest_ref(:,j) .and. ee_test_mask(:,j))) then
          write(unit,'(a,i0,a,1x,i0)')  'shape[', j, ']=', old_mesh(j)%cell_shape
          write(unit,'(a,i0,a,6(1x,g0))') 'new[', j, ']=', dest(:,j)
          write(unit,'(a,i0,a,6(1x,g0))') 'old[', j, ']=', dest_ref(:,j)
        end if
      end do
      call TLS_fatal ('error validating ee_gather_real64')
    end if
#else
    INSIST(global_all(dest == dest_ref .or. .not.ee_test_mask))
#endif

  end subroutine ee_gather_validate_real64

  subroutine ee_gather_validate_int32 (dest, src)

    use truchas_logging_services
    use parallel_communication, only: global_any, global_all
    use gs_module, only: ee_gather_old => ee_gather
    use mesh_module, only: old_mesh => mesh

    integer, intent(inout) :: dest(:,:)
    integer, intent(in) :: src(:)

    integer :: j, unit
    integer, allocatable :: dest_ref(:,:)

    dest_ref = dest ! dest is intent(inout); capture it before it is modified
    call ee_gather_old (dest_ref, src)

    call ee_gather_int32 (dest, src)

!#define DEBUG
#ifdef DEBUG
#undef DEBUG
    if (global_any(dest /= dest_ref .and. ee_test_mask)) then
      unit = TLS_debug_unit()
      call TLS_debug ('')
      call TLS_debug ('EE_GATHER_INT32 ***********************************************')
      do j = 1, size(dest,2)
        if (any(dest(:,j) /= dest_ref(:,j) .and. ee_test_mask(:,j))) then
          write(unit,'(a,i0,a,1x,i0)')  'shape[', j, ']=', old_mesh(j)%cell_shape
          write(unit,'(a,i0,a,6(1x,i0))') 'new[', j, ']=', dest(:,j)
          write(unit,'(a,i0,a,6(1x,i0))') 'old[', j, ']=', dest_ref(:,j)
        end if
      end do
      call TLS_fatal ('error validating ee_gather_int32')
    end if
#else
    INSIST(global_all(dest == dest_ref .or. .not.ee_test_mask))
#endif

  end subroutine ee_gather_validate_int32

  subroutine ee_gather_validate_log (dest, src)

    use truchas_logging_services
    use parallel_communication, only: global_any, global_all
    use gs_module, only: ee_gather_old => ee_gather
    use mesh_module, only: old_mesh => mesh

    logical, intent(inout) :: dest(:,:)
    logical, intent(in) :: src(:)

    integer :: j, unit
    logical, allocatable :: dest_ref(:,:)

    dest_ref = dest ! dest is intent(inout); capture it before it is modified
    call ee_gather_old (dest_ref, src)

    call ee_gather_log (dest, src)

!#define DEBUG
#ifdef DEBUG
#undef DEBUG
    if (global_any((dest .neqv. dest_ref) .and. ee_test_mask_log)) then
      unit = TLS_debug_unit()
      call TLS_debug ('')
      call TLS_debug ('EE_GATHER_LOG ***********************************************')
      do j = 1, size(dest,2)
        if (any((dest(:,j) .neqv. dest_ref(:,j)) .and. ee_test_mask(:,j))) then
          write(unit,'(a,i0,a,1x,i0)')  'shape[', j, ']=', old_mesh(j)%cell_shape
          write(unit,'(a,i0,a,6(1x,l1))') 'new[', j, ']=', dest(:,j)
          write(unit,'(a,i0,a,6(1x,l1))') 'old[', j, ']=', dest_ref(:,j)
        end if
      end do
      call TLS_fatal ('error validating ee_gather_log')
    end if
#else
    INSIST(global_all((dest .eqv. dest_ref) .or. .not.ee_test_mask_log))
#endif

  end subroutine ee_gather_validate_log

  subroutine ss_gather_validate_real64 (dest, src)

    use truchas_logging_services
    use parallel_communication, only: global_any, global_all
    use gs_module, only: ee_gather_old => ee_gather
    use mesh_module, only: old_mesh => mesh

    real(r8), intent(inout) :: dest(:,:)
    real(r8), intent(in) :: src(:,:)

    integer :: j, unit
    real(r8), allocatable :: dest_ref(:,:)

    dest_ref = dest ! dest is intent(inout); capture it before it is modified
    call ee_gather_old (dest_ref, src)

    call ss_gather_real64 (dest, src)

!#define DEBUG
#ifdef DEBUG
#undef DEBUG
    !if (global_any(dest /= dest_ref .and. ee_test_mask)) then
    if (global_any(dest /= dest_ref)) then
      unit = TLS_debug_unit()
      call TLS_debug ('')
      call TLS_debug ('SS_GATHER_REAL64 ***********************************************')
      do j = 1, size(dest,2)
        !if (any(dest(:,j) /= dest_ref(:,j) .and. ee_test_mask(:,j))) then
        if (any(dest(:,j) /= dest_ref(:,j))) then
          write(unit,'(a,i0,a,1x,i0)')  'shape[', j, ']=', old_mesh(j)%cell_shape
          write(unit,'(a,i0,a,6(1x,g0))') 'new[', j, ']=', dest(:,j)
          write(unit,'(a,i0,a,6(1x,g0))') 'old[', j, ']=', dest_ref(:,j)
        end if
      end do
      call TLS_fatal ('error validating ss_gather_real64')
    end if
#else
    INSIST(global_all(dest == dest_ref))
#endif

  end subroutine ss_gather_validate_real64

  subroutine ss_gather_validate_log (dest, src)

    use truchas_logging_services
    use parallel_communication, only: global_any, global_all
    use gs_module, only: ee_gather_old => ee_gather
    use mesh_module, only: old_mesh => mesh

    logical, intent(inout) :: dest(:,:)
    logical, intent(in) :: src(:,:)

    integer :: j, unit
    logical, allocatable :: dest_ref(:,:)

    dest_ref = dest ! dest is intent(inout); capture it before it is modified
    call ee_gather_old (dest_ref, src)

    call ss_gather_log (dest, src)

    !! Regression tests pass without masking, which suggests that the
    !! differences excluded by EE_TEST_MASK are not significant.
!#define DEBUG
#ifdef DEBUG
#undef DEBUG
    if (global_any((dest .neqv. dest_ref) .and. ee_test_mask)) then
      unit = TLS_debug_unit()
      call TLS_debug ('')
      call TLS_debug ('SS_GATHER_LOG ***********************************************')
      do j = 1, size(dest,2)
        if (any((dest(:,j) .neqv. dest_ref(:,j)) .and. ee_test_mask(:,j))) then
          write(unit,'(a,i0,a,1x,i0)')  'shape[', j, ']=', old_mesh(j)%cell_shape
          write(unit,'(a,i0,a,6(1x,g0))') 'new[', j, ']=', dest(:,j)
          write(unit,'(a,i0,a,6(1x,g0))') 'old[', j, ']=', dest_ref(:,j)
        end if
      end do
      call TLS_fatal ('error validating ss_gather_log')
    end if
#else
    INSIST(global_all((dest .eqv. dest_ref) .or. .not.ee_test_mask))
#endif

  end subroutine ss_gather_validate_log

  subroutine ee_gather_all_v_s_validating (dest, source, boundary, range)

    use var_vector_module, only: real_var_vector, sizes, create, flatten
    use gs_module, only: ee_gather_old => ee_gather
    use parallel_communication, only: global_any, global_all
    use truchas_logging_services

    integer :: j, k, unit
    real(r8) :: source(:)
    type(real_var_vector) :: dest(:)
    real(r8), pointer, optional :: boundary(:)
    integer, optional :: range(:)
    real(r8), pointer :: cont(:), ref_cont(:)

    type(real_var_vector) :: ref_dest(size(dest))
    real(r8), pointer :: buffer(:)

    call create (ref_dest, sizes(dest))

    if (present(boundary)) then
      buffer => boundary
      call ee_gather_old (ref_dest, source, buffer, range)
      if (.not.associated(buffer, boundary)) deallocate(buffer)
    else
      call ee_gather_old (ref_dest, source, boundary, range)
    end if

    call ee_gather_all_v_s_real64 (dest, source, boundary, range)

!#define DEBUG
#ifdef DEBUG
#undef DEBUG
    cont => flatten(dest); ref_cont => flatten(ref_dest)
    if (global_any(cont /= ref_cont)) then
      unit = TLS_debug_unit()
      call TLS_debug ('')
      call TLS_debug ('EE_GATHER_ALL_V_S ***************************')
      do j = 1, size(dest)
        associate (old => ref_dest(j)%v, new => dest(j)%v)
          do k = 1, size(old)
            if (old(k) /= new(k)) then
              write(unit,'(2(a,i0),a,2(g0,1x),a)') '[', k, ',', j, ']=', old(k), new(k), ' (old, new)'
            end if
          end do
        end associate
      end do
      call TLS_fatal ('error validating ee_gather_all_v_s')
    end if
#else
    cont => flatten(dest); ref_cont => flatten(ref_dest)
    INSIST(global_all(cont == ref_cont))
#endif

    deallocate(ref_cont)

  end subroutine ee_gather_all_v_s_validating

end module
