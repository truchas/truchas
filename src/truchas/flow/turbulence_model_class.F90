!!
!! TURBULENCE_MODEL_CLASS
!!
!! This module provides an abstract type that encapsulates truchas turbulence models
!!
!! Peter Brady <ptb@lanl.gov>
!! 2018
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  Subtypes should provide implementations of the following methods
!!
!!  READ_PARAMS(PARAMETER_LIST P) - Initializes the parameter list for this object.
!!    This must be call first for the object to behave sensibly.  Error checking of
!!    parameters should be done here
!!
!!  INIT(FLOW_MESH M) - Allocates the data members for this object and initializes
!!    internal state
!!
!!  SETUP(REAL(R8) CELL_VELOCITY(:))
!!    perform computations in then model for a given cell-centered velocity
!!
!!  APPLY(REAL(R8) CELL_VISCOSITY(:))
!!    updates the cell-centered viscosity based on model
!!
!!  ACCEPT()
!!    updates any internal state after timestep is accepted


module turbulence_model_class
  use kinds, only: r8
  use unstr_mesh_type
  use flow_props_type
  use parameter_list_type
  implicit none
  private

  public :: turbulence_model

  type, abstract :: turbulence_model
    type(unstr_mesh), pointer :: mesh => null()
  contains
    procedure(read_params), deferred :: read_params
    procedure(init), deferred :: init
    procedure(setup), deferred :: setup
    procedure(apply), deferred :: apply
    procedure(accept), deferred :: accept
  end type turbulence_model

  abstract interface
    subroutine read_params(this, params)
      import
      class(turbulence_model), intent(inout) :: this
      type(parameter_list), pointer, intent(in) :: params
    end subroutine read_params

    subroutine init(this, mesh)
      import
      class(turbulence_model), intent(inout) :: this
      type(unstr_mesh), intent(in), target :: mesh
    end subroutine init

    subroutine setup(this, vel_cc)
      import
      class(turbulence_model), intent(inout) :: this
      real(r8), intent(in) :: vel_cc(:,:)
    end subroutine setup

    subroutine apply(this, props)
      import
      class(turbulence_model), intent(inout) :: this
      type(flow_props), intent(inout) :: props
    end subroutine apply

    subroutine accept(this)
      import
      class(turbulence_model), intent(inout) :: this
    end subroutine accept
  end interface
end module turbulence_model_class
