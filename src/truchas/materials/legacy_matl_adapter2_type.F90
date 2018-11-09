#include "f90_assert.fpp"

module legacy_matl_adapter2_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use legacy_matl_adapter_class
  use integer_set_type
  implicit none
  private

  type, extends(legacy_matl_adapter), public :: legacy_matl_adapter2
    integer :: nmat, ncells
    real(r8), allocatable :: vfrac(:)
    integer, allocatable :: matid(:), index(:)
    ! temporaries used during construction
    type(integer_set), allocatable :: mset(:)
  contains
    procedure :: redef
    procedure :: enddef
    procedure :: add_cell_mat
    procedure :: gather_vof
    procedure :: get_vof
    procedure :: set_vof
    procedure :: update_vof
    procedure :: read_data
    procedure :: get_cell_vof
    procedure :: set_cell_vof
    procedure :: normalize_vof
    procedure :: get_solid_mask
  end type

contains

  subroutine redef(this)
    use legacy_mesh_api, only: ncells
    use matl_module, only: nmat
    class(legacy_matl_adapter2), intent(inout) :: this
    this%nmat = nmat
    this%ncells = ncells
    allocate(this%mset(ncells))
  end subroutine

  subroutine add_cell_mat(this, cellid, matid)
    class(legacy_matl_adapter2), intent(inout) :: this
    integer, intent(in) :: cellid, matid
    ASSERT(cellid > 0 .and. cellid <= this%ncells)
    ASSERT(matid > 0 .and. matid <= this%nmat)
    call this%mset(cellid)%add(matid)
  end subroutine

  subroutine enddef(this)
    class(legacy_matl_adapter2), intent(inout) :: this
    integer :: j, n
    integer, allocatable :: list(:)
    allocate(this%index(this%ncells+1))
    this%index(1) = 1
    do j = 1, this%ncells
      this%index(j+1) = this%index(j) + this%mset(j)%size()
    end do
    n = this%index(this%ncells+1) - 1
    allocate(this%matid(n))
    do j = 1, this%ncells
      list = this%mset(j) ! lhs must be allocatable :-(
      this%matid(this%index(j):this%index(j+1)-1) = list
    end do
    deallocate(this%mset)
    allocate(this%vfrac(n))
    this%vfrac = 0
  end subroutine

  subroutine gather_vof(this, m, vof)
    class(legacy_matl_adapter2), intent(in) :: this
    integer, intent(in) :: m
    real(r8), intent(out) :: vof(:)
    integer :: i, j
    vof = 0
    do j = 1, this%ncells
      associate (matid => this%matid(this%index(j):this%index(j+1)-1), &
                 vfrac => this%vfrac(this%index(j):this%index(j+1)-1))
        do i = 1, size(matid)
          if (matid(i) == m) then
            vof(j) = vfrac(i)
            exit
          end if
          if (matid(i) > m) exit
        end do
      end associate
    end do
  end subroutine

  subroutine get_vof(this, vof)
    class(legacy_matl_adapter2), intent(in) :: this
    real(r8), intent(out) :: vof(:,:)
    integer :: i, j
    vof = 0
    do j = 1, this%ncells
      associate (matid => this%matid(this%index(j):this%index(j+1)-1), &
                 vfrac => this%vfrac(this%index(j):this%index(j+1)-1))
        do i = 1, size(matid)
          vof(matid(i),j) = vfrac(i)
        end do
      end associate
    end do
  end subroutine

  ! NB: this ignores any vof element not contained in the compressed structure
  subroutine set_vof(this, vof)
    class(legacy_matl_adapter2), intent(inout) :: this
    real(r8), intent(in) :: vof(:,:)
    integer :: i, j
    do j = 1, this%ncells
      associate (matid => this%matid(this%index(j):this%index(j+1)-1), &
                 vfrac => this%vfrac(this%index(j):this%index(j+1)-1))
        do i = 1, size(matid)
          vfrac(i) = vof(matid(i),j)
        end do
      end associate
    end do
  end subroutine

  subroutine update_vof(this, vof)
    class(legacy_matl_adapter2), intent(inout) :: this
    real(r8), intent(in) :: vof(0:,:)
    call set_vof(this, vof(1:,:))
  end subroutine

  subroutine get_cell_vof(this, n, vof)
    class(legacy_matl_adapter2), intent(in) :: this
    integer, intent(in) :: n
    real(r8), intent(out) :: vof(:)
    integer :: i
    vof = 0
    do i = this%index(n), this%index(n+1) - 1
      vof(this%matid(i)) = this%vfrac(i)
    end do
  end subroutine

  subroutine set_cell_vof(this, n, m, vfrac)
    class(legacy_matl_adapter2), intent(inout) :: this
    integer, intent(in) :: n, m
    real(r8), intent(in) :: vfrac
    integer :: i
    do i = this%index(n), this%index(n+1) - 1
      if (this%matid(i) == m) then
        this%vfrac(i) = vfrac
        return
      end if
      if (this%matid(i) > m) exit
    end do
    ASSERT(.false.)
  end subroutine

  !! Possibly called as a final fixup from init_module::matl_init.
  !! Should not be needed.
  subroutine normalize_vof(this)
    class(legacy_matl_adapter2), intent(inout) :: this
    integer :: j
    do j = 1, this%ncells
      associate (vfrac => this%vfrac(this%index(j):this%index(j+1)-1))
        vfrac = vfrac / sum(vfrac)
      end associate
    end do
  end subroutine

  !! Extracted from solid_mechanics_module::solid_mech_init and thermo_mechanics.
  !! NB: Volume fraction data is never referenced!
  !! NNC: Does the matl data structure always squeeze out a material with 0 volume
  !! fraction (by setting the ID to 0)? I am unconvinced that it does. This seems
  !! suspect to me. Really just impacts solidifying material.

  subroutine get_solid_mask(this, mask)
    use fluid_data_module, only: isImmobile
    class(legacy_matl_adapter2), intent(in) :: this
    logical, intent(out) :: mask(:)
    integer :: i, j
    ASSERT(size(mask) == this%ncells)
    CELL: do j = 1, this%ncells
      mask(j) = .true.
      do i = this%index(j), this%index(j+1) - 1
        if (.not.isImmobile(this%matid(i))) cycle
        if (this%vfrac(i) > 0) cycle CELL
      end do
      mask(j) = .false.
    end do CELL
  end subroutine

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! READ_MATL_DATA
 !!
 !! Neil N. Carlson <nnc@lanl.gov>
 !! 18 Apr 2005
 !!
 !! This subroutine reads the volume fraction data from a restart file opened
 !! (and pre-positioned) on UNIT, and uses this data (properly distributed and
 !! permuted) to define the module structure MATL.  VERSION is the version
 !! number of the restart file format.
 !!
 !! NB: It is implicitly assumed that the material indices in the restart file
 !! directly correspond to the material indices generated by the input file.
 !! We require that the number of materials agree between the input and restart
 !! files; this constraint could probably be relaxed to allow the input to
 !! specify more (new additional) materials.
 !!

  subroutine read_data(this, unit, version)

    use matl_module, only: nmat
    use legacy_mesh_api, only: ncells, pcell => unpermute_mesh_vector
    use restart_utilities, only: read_var, read_dist_array, halt
    use string_utilities, only: i_to_c

    class(legacy_matl_adapter2), intent(inout) :: this
    integer, intent(in) :: unit, version

    integer :: n

    this%nmat = nmat
    this%ncells = ncells

    !! Read the number of materials defined in the restart file.
    call read_var(unit, n, 'READ_MATL_DATA: error reading NMAT record')
    if (n /= nmat) call halt('READ_MATL_DATA: incompatible NMAT value: ' // i_to_c(n))

    !! Read the volume fraction array.
!    allocate(this%vfrac(this%nmat,this%ncells))
!    call read_dist_array(unit, this%vfrac, pcell, 'READ_MATL_DATA: error reading VF records')

  end subroutine

end module legacy_matl_adapter2_type
