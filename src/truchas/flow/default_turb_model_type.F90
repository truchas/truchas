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
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Does nothing

module default_turb_model_type
  use kinds
  use turbulence_model_class
  use parameter_list_type
  use flow_mesh_type
  use flow_props_type
  implicit none
  private

  public :: default_turb_model, alloc_default_turb_model

  type, extends(turbulence_model) :: default_turb_model
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

    ! nothing
  end subroutine read_params

  subroutine init(this, mesh)
    class(default_turb_model), intent(inout) :: this
    type(flow_mesh), pointer, intent(in) :: mesh

    ! nothing
  end subroutine init

  subroutine setup(this, vel_cc)
    class(default_turb_model), intent(inout) :: this
    real(r8), intent(in) :: vel_cc(:,:)

    ! nothing
  end subroutine setup

  subroutine apply(this, props)
    class(default_turb_model), intent(inout) :: this
    type(flow_props), intent(inout) :: props

    !nothing
  end subroutine apply

  subroutine accept(this)
    class(default_turb_model), intent(inout) :: this

    !nothing
  end subroutine accept

end module default_turb_model_type
