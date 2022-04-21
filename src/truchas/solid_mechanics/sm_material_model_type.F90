!!
!! Zach Jibben <zjibben@lanl.gov>
!! October 2021
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module sm_material_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use scalar_func_containers
  use viscoplastic_material_model_types
  implicit none
  private

  type, public :: sm_material_model
    integer :: nmat
    logical :: viscoplasticity_enabled = .false.
    type(scalar_func_box), allocatable :: lame1f(:), lame2f(:), densityf(:)
    real(r8), allocatable :: reference_density(:)
    type(viscoplastic_material_model_box), allocatable :: vp(:)
  contains
    procedure :: init
  end type sm_material_model

contains

  subroutine init(this, params, nmat, lame1f, lame2f, densityf, reference_density, vp)

    class(sm_material_model), intent(out) :: this
    type(parameter_list), intent(inout) :: params
    integer, intent(in) :: nmat
    type(scalar_func_box), allocatable, intent(inout) :: lame1f(:), lame2f(:), densityf(:)
    real(r8), intent(in) :: reference_density(:)
    type(viscoplastic_material_model_box), allocatable, intent(inout) :: vp(:)

    integer :: m

    ASSERT(nmat > 0)
    ASSERT(size(lame1f) == nmat)
    ASSERT(size(lame2f) == nmat)
    ASSERT(size(densityf) == nmat)
    ASSERT(size(reference_density) == nmat)
    ASSERT(size(vp) == nmat)

    this%nmat = nmat
    this%reference_density = reference_density
    call move_alloc(lame1f, this%lame1f)
    call move_alloc(lame2f, this%lame2f)
    call move_alloc(densityf, this%densityf)
    call move_alloc(vp, this%vp)

    this%viscoplasticity_enabled = .false.
    do m = 1, this%nmat
      !if (allocated(this%vp(m)%m)) this%viscoplasticity_enabled = .true.
      if (.not.this%vp(m)%m%is_elastic) this%viscoplasticity_enabled = .true.
    end do

  end subroutine init

end module sm_material_model_type
