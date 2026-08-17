!!
!! CELL_MATL_PROP_FUNC_TYPE
!!
!! This module defines the CELL_MATL_PROP_FUNC type for evaluating cell
!! material properties by volume-fraction averaging over the simulation-owned
!! MATERIAL_COMPOSITION.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module cell_matl_prop_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use new_mesh_func_class
  use material_composition_type
  use avg_matl_prop_type
  implicit none
  private

  type, extends(new_mesh_func), public :: cell_matl_prop_func
    private
    type(material_composition), pointer :: composition => null()
    type(avg_matl_prop) :: prop
    integer, allocatable :: mid(:)
  contains
    procedure :: init
    procedure :: compute_value
    procedure :: compute_value_cell
    procedure :: compute_deriv
  end type cell_matl_prop_func

contains

  subroutine init(this, matl_model, composition, name, stat, errmsg)
    use material_model_type

    class(cell_matl_prop_func), intent(out) :: this
    type(material_model), intent(in) :: matl_model
    type(material_composition), pointer, intent(in) :: composition
    character(*), intent(in) :: name
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: i

    allocate(this%mid(matl_model%nmatl_real))
    this%mid = [(i, i=1,matl_model%nmatl_real)]
    call this%prop%init(name, this%mid, matl_model, stat, errmsg)
    if (stat /= 0) return
    this%composition => composition
  end subroutine init

  subroutine compute_value(this, state, value)
    class(cell_matl_prop_func), intent(in) :: this
    real(r8), intent(in) :: state(:)
    real(r8), intent(out) :: value(:)
    integer :: j

    ASSERT(size(state,1) >= size(this%composition%vfrac,2))
    ASSERT(size(value) == size(this%composition%vfrac,2))
    do j = 1, size(this%composition%vfrac,2)
      call this%prop%compute_value(this%composition%vfrac(this%mid,j), state(j:j), value(j))
    end do
  end subroutine compute_value

  subroutine compute_value_cell(this, cell, state, value)
    class(cell_matl_prop_func), intent(in) :: this
    integer, intent(in) :: cell
    real(r8), intent(in) :: state
    real(r8), intent(out) :: value

    ASSERT(cell >= 1 .and. cell <= size(this%composition%vfrac,2))
    call this%prop%compute_value(this%composition%vfrac(this%mid,cell), [state], value)
  end subroutine compute_value_cell

  subroutine compute_deriv(this, state, index, value)
    class(cell_matl_prop_func), intent(in) :: this
    real(r8), intent(in) :: state(:)
    integer, intent(in) :: index
    real(r8), intent(out) :: value(:)
    integer :: j

    ASSERT(size(state,1) >= size(this%composition%vfrac,2))
    ASSERT(size(value) == size(this%composition%vfrac,2))
    do j = 1, size(this%composition%vfrac,2)
      call this%prop%compute_deriv(this%composition%vfrac(this%mid,j), state(j:j), index, value(j))
    end do
  end subroutine compute_deriv

end module cell_matl_prop_func_type
