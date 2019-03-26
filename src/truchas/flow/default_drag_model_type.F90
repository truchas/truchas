!!
!! DEFAULT_DRAG_MODEL_TYPE
!!
!! This module provides an implementation of the abstract dragulence_model_class
!!
!! Peter Brady <ptb@lanl.gov>
!! May 2018
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Does nothing

module default_drag_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use porous_drag_model_class
  use parameter_list_type
  implicit none
  private

  public :: alloc_default_drag_model

  type, extends(porous_drag_model), public :: default_drag_model
  end type default_drag_model

contains

  subroutine alloc_default_drag_model(t)
    class(porous_drag_model), allocatable, intent(out) :: t
    type(default_drag_model), allocatable :: m
    allocate(m)
    call move_alloc(m, t)
  end subroutine alloc_default_drag_model

  subroutine read_params(this, params)
    class(default_drag_model), intent(inout) :: this
    type(parameter_list), pointer, intent(in) :: params
    ! nothing
  end subroutine read_params

  subroutine init(this, mesh)
    class(default_drag_model), intent(inout) :: this
    type(unstr_mesh), intent(in), target :: mesh
    ! nothing
  end subroutine init

  subroutine setup(this, vel_cc)
    class(default_drag_model), intent(inout) :: this
    real(r8), intent(in) :: vel_cc(:,:)
    ! nothing
  end subroutine setup

  subroutine apply(this, visc_cc)
    class(default_drag_model), intent(inout) :: this
    real(r8), intent(inout) :: visc_cc(:)
    !nothing
  end subroutine apply

  subroutine accept(this)
    class(default_drag_model), intent(inout) :: this
    !nothing
  end subroutine accept

end module default_drag_model_type
