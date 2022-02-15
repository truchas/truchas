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

  public :: ee_gather, gather_boundarydata

  interface ee_gather
    module procedure ee_gather_int32
    module procedure ee_gather_real64
    module procedure ee_gather_log
    module procedure ss_gather_real64
    module procedure ss_gather_log
    module procedure ee_gather_all_v_s_real64
  end interface

  interface gather_boundarydata
    module procedure gather_boundarydata_real64
  end interface

contains

  subroutine ee_gather_real64 (dest, src, boundary)

    use mesh_impl, only: mesh, legacy_cell_ip, DEGENERATE_FACE

    real(r8), intent(inout) :: dest(:,:)
    real(r8), intent(in) :: src(:)
    real(r8), pointer, optional :: boundary(:)

    integer :: j, k, n
    real(r8), pointer :: buffer(:)

    if (present(boundary)) then
      if (.not.associated(boundary)) then
        allocate(boundary(legacy_cell_ip%offp_size))
        call legacy_cell_ip%gather_offp(src, boundary)
      end if
      buffer => boundary
    else
      allocate(buffer(legacy_cell_ip%offp_size))
      call legacy_cell_ip%gather_offp(src, buffer)
    end if

    do j = 1, ncells
      do k = 1, 6
        n = mesh(j)%ngbr_cell(k)
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

    use mesh_impl, only: mesh, legacy_cell_ip, DEGENERATE_FACE

    integer, intent(inout) :: dest(:,:)
    integer, intent(in) :: src(:)

    integer :: j, k, n
    integer :: buffer(legacy_cell_ip%offp_size)

    call legacy_cell_ip%gather_offp(src, buffer)

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

    use mesh_impl, only: mesh, legacy_cell_ip, DEGENERATE_FACE

    logical, intent(inout) :: dest(:,:)
    logical, intent(in) :: src(:)

    integer :: j, k, n
    logical :: buffer(legacy_cell_ip%offp_size)

    call legacy_cell_ip%gather_offp(src, buffer)

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

    use mesh_impl, only: mesh, legacy_cell_ip, DEGENERATE_FACE

    real(r8), intent(inout) :: dest(:,:)
    real(r8), intent(in) :: src(:,:)

    integer :: j, jj, k, kk
    real(r8) :: buffer(6,legacy_cell_ip%offp_size)

    ASSERT(size(dest,1) == 6)
    ASSERT(size(dest,2) == ncells)
    ASSERT(size(src,1) == 6)
    ASSERT(size(src,2) == ncells)

    call legacy_cell_ip%gather_offp(src, buffer)

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

    use mesh_impl, only: mesh, legacy_cell_ip, DEGENERATE_FACE

    logical, intent(inout) :: dest(:,:)
    logical, intent(in) :: src(:,:)

    integer :: j, jj, k, kk
    logical :: buffer(6,legacy_cell_ip%offp_size)

    ASSERT(size(dest,1) == 6)
    ASSERT(size(dest,2) == ncells)
    ASSERT(size(src,1) == 6)
    ASSERT(size(src,2) == ncells)

    call legacy_cell_ip%gather_offp(src, buffer)

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

    use mesh_impl, only: mesh
    use var_vector_types, only: real_var_vector

    real(r8) :: source(:)
    type(real_var_vector) :: dest(:)
    real(r8), pointer, optional :: boundary(:)
    integer, optional :: range(:)

    integer :: j, k
    real(r8), pointer :: buffer(:) => null()

    if (present(boundary)) buffer => boundary
    if (.not.associated(buffer)) then
      allocate(buffer(new_mesh%cell_ip%offp_size))
      call new_mesh%cell_ip%gather_offp(source, buffer)
      if (present(boundary)) boundary => buffer
    end if

    ASSERT(size(source) == new_mesh%cell_ip%onp_size)
    ASSERT(size(buffer) == new_mesh%cell_ip%offp_size)
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
    use mesh_impl, only: legacy_cell_ip
    real(r8), pointer :: boundary(:)
    real(r8), intent(in) :: source(:)
    ASSERT(size(source) == legacy_cell_ip%onp_size)
    if (associated(boundary)) return  ! signal to do nothing
    allocate(boundary(legacy_cell_ip%offp_size))
    call legacy_cell_ip%gather_offp(source, boundary)
  end subroutine gather_boundarydata_real64

end module
