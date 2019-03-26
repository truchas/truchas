!!
!! DEFAULT_TURB_MODEL_TYPE
!!
!! This module provides an implementation of the abstract turbulence_model_class
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

module default_turb_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use turbulence_model_class
  use parameter_list_type
  use unstr_mesh_type
  use flow_props_type
  implicit none
  private

  public :: alloc_default_turb_model

  type, extends(turbulence_model), public :: default_turb_model
  contains
    procedure :: read_params
    procedure :: init
    procedure :: setup
    procedure :: apply
    procedure :: accept
  end type default_turb_model

contains

  subroutine alloc_default_turb_model(t)
    class(turbulence_model), allocatable, intent(out) :: t
    type(default_turb_model), allocatable :: m
    allocate(m)
    call move_alloc(m, t)
  end subroutine alloc_default_turb_model

  subroutine read_params(this, params)
    class(default_turb_model), intent(inout) :: this
    type(parameter_list), pointer, intent(in) :: params
  end subroutine read_params

  subroutine init(this, mesh)
    class(default_turb_model), intent(inout) :: this
    type(unstr_mesh), intent(in), target :: mesh
  end subroutine init

  subroutine setup(this, vel_cc)
    class(default_turb_model), intent(inout) :: this
    real(r8), intent(in) :: vel_cc(:,:)
  end subroutine setup

  subroutine apply(this, props)
    class(default_turb_model), intent(inout) :: this
    type(flow_props), intent(inout) :: props
  end subroutine apply

  subroutine accept(this)
    class(default_turb_model), intent(inout) :: this
  end subroutine accept

end module default_turb_model_type
