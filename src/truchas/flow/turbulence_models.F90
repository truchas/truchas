!!
!! TURBULENCE_MODELS
!!
!! Include this module to use any of the turbulence models provided by truchas.
!!
!! Peter Brady <ptb@lanl.gov>
!! 2018
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PROGRAMMING INTERFACE
!!
!!  READ_FLOW_TURBULENCE_MODEL_NAMELIST (LUN, PARAMETER_LIST)
!!    reads namelists /flow_turbulence_model/ and stores into parameter list
!!
!!  TURBULENCE_MODELS_READ_PARAMS(model, PARAMETER_LIST)
!!    Allocates a turbulence_model based on parameter_list
!!

module turbulence_models

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use truchas_logging_services
  use turbulence_model_class
  use algebraic_turb_model_type
  use default_turb_model_type
  use parameter_list_type
  implicit none
  private

  public :: turbulence_model, alloc_turbulence_model

contains

  subroutine alloc_turbulence_model(turb, params, off)

    class(turbulence_model), allocatable, intent(out) :: turb
    type(parameter_list), pointer, intent(in) :: params
    logical, optional, intent(in) :: off

    integer :: stat
    character(:), allocatable :: model
    type(parameter_list), pointer :: pp
    logical ::off_

    if (present(off)) then
      off_ = off
    else
      off_ = .false.
    end if

    pp => params%sublist("turbulence model")

    call pp%get('type', model, stat=stat)

    if (stat /= 0 .or. off_) then
      call alloc_default_turb_model(turb)
    else
      select case (model)
      case ("alg")
        call alloc_algebraic_turb_model(turb)
      case default
        call alloc_default_turb_model(turb)
      end select
    end if

    pp => pp%sublist("params")
    call turb%read_params(pp)

  end subroutine alloc_turbulence_model

end module turbulence_models
