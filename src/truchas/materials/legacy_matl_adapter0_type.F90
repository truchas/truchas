#include "f90_assert.fpp"

module legacy_matl_adapter0_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use legacy_matl_adapter_class
  use graph_type
  implicit none
  private

  type, extends(legacy_matl_adapter), public :: legacy_matl_adapter0
    integer :: nmat, ncells
    type(graph), allocatable :: g
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
    class(legacy_matl_adapter0), intent(inout) :: this
    integer :: n
    this%nmat = nmat
    this%ncells = ncells
    if (allocated(this%g)) deallocate(this%g)
    allocate(this%g)
    n = max(ncells, nmat)
    call this%g%init(n, directed=.true., self_edge=.true.)
  end subroutine redef

  subroutine add_cell_mat(this, cellid, matid)
    class(legacy_matl_adapter0), intent(inout) :: this
    integer, intent(in) :: cellid, matid
    ASSERT(allocated(this%g))
    ASSERT(cellid > 0 .and. cellid <= this%ncells)
    ASSERT(matid > 0 .and. matid <= this%nmat)
    call this%g%add_edge(cellid, matid)
  end subroutine add_cell_mat

  subroutine enddef(this)

    use matl_module, only: matl, mat_slot, mat_slot_new, slot_increase, slot_set
    class(legacy_matl_adapter0), intent(inout) :: this

    integer :: i, j
    integer, allocatable :: xadj(:), adjncy(:)

    ASSERT(allocated(this%g))

    call this%g%get_adjacency(xadj, adjncy)
    deallocate(this%g)

    !! Allocate the required number of slots and zero it all out
    mat_slot_new = maxval(xadj(2:this%ncells+1) - xadj(1:this%ncells))
    if (mat_slot_new > mat_slot) call slot_increase(matl, mat_slot, mat_slot_new)
    do i = 1, mat_slot
      call slot_set(matl, i)
    end do

    !! Definte the material ids
    do j = 1, this%ncells
      associate (list => adjncy(xadj(j):xadj(j+1)-1))
        do i = 1, size(list)
          matl(i)%cell(j)%id = list(i)
        end do
      end associate
    end do

  end subroutine enddef

  subroutine gather_vof(this, m, vof)
    use matl_module, only: matl_gather_vof => gather_vof
    class(legacy_matl_adapter0), intent(in) :: this
    integer, intent(in) :: m
    real(r8), intent(out) :: vof(:)
    call matl_gather_vof(m, vof)
  end subroutine

  subroutine get_vof(this, vof)
    use matl_utilities, only: matl_get_vof
    class(legacy_matl_adapter0), intent(in) :: this
    real(r8), intent(out) :: vof(:,:)
    call matl_get_vof(vof)
  end subroutine

  subroutine set_vof(this, vof)
    use matl_utilities, only: matl_set_vof
    class(legacy_matl_adapter0), intent(inout) :: this
    real(r8), intent(in) :: vof(:,:)
    call matl_set_vof(vof)
  end subroutine

  subroutine update_vof(this, vof)
    use matl_utilities, only: update_matl
    class(legacy_matl_adapter0), intent(inout) :: this
    real(r8), intent(in) :: vof(:,:)
    call update_matl(vof)
  end subroutine

  subroutine read_data(this, unit, version)
    use matl_utilities, only: read_matl_data
    class(legacy_matl_adapter0), intent(inout) :: this
    integer, intent(in) :: unit, version
    call read_matl_data(unit, version)
  end subroutine

  subroutine get_cell_vof(this, n, vof)
    use matl_utilities, only: matl_get_cell_vof
    class(legacy_matl_adapter0), intent(in) :: this
    integer, intent(in) :: n
    real(r8), intent(out) :: vof(:)
    call matl_get_cell_vof(n, vof)
  end subroutine

  subroutine set_cell_vof(this, n, m, vfrac)
    use matl_module, only: matl, mat_slot
    class(legacy_matl_adapter0), intent(inout) :: this
    integer, intent(in) :: n, m
    real(r8), intent(in) :: vfrac
    integer :: s
    do s = 1, mat_slot
      if (matl(s)%cell(n)%id == m) then
        matl(s)%cell(n)%vof = vfrac
        return
      end if
    end do
    INSIST(.false.)
  end subroutine

  !! Possibly called as a final fixup from init_module::matl_init.
  !! Should not be needed.
  subroutine normalize_vof(this)
    use matl_module, only: matl, mat_slot
    class(legacy_matl_adapter0), intent(inout) :: this
    integer :: s
    real(r8) :: total(this%ncells)
    total = 0.0_r8
    do s = 1, mat_slot
      total = total + matl(s)%cell%vof
    end do
    do s = 1, mat_slot
      matl(s)%cell%vof = matl(s)%cell%vof / total
    end do
  end subroutine

  !! Extracted from solid_mechanics_module::solid_mech_init and thermo_mechanics.
  !! NB: Volume fraction data is never referenced!
  !! NNC: Does the matl data structure always squeeze out a material with 0 volume
  !! fraction (by setting the ID to 0)? I am unconvinced that it does. This seems
  !! suspect to me. Really just impacts solidifying material.

  subroutine get_solid_mask(this, mask)
    use matl_module, only: matl, mat_slot
    use fluid_data_module, only: isImmobile
    class(legacy_matl_adapter0), intent(in) :: this
    logical, intent(out) :: mask(:)
    integer :: j, m, s
    ASSERT(size(mask) == this%ncells)
    mask = .false.
    do j = 1, this%ncells
      do m = 1, this%nmat
        if (isImmobile(m)) then
          do s = 1, mat_slot
            if (matl(s)%cell(j)%id == m) mask(j) = .true.
          end do
        end if
      end do
    end do
  end subroutine

end module legacy_matl_adapter0_type
