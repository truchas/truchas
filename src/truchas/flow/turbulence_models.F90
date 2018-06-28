module turbulence_models
  use turbulence_model_class
  use algebraic_turb_model_type
  use default_turb_model_type
  use parameter_list_type
  implicit none
  private

  public :: turbulence_model, turbulence_models_read_params

contains

  subroutine turbulence_models_read_params(turb, params, off)
    class(turbulence_model), allocatable, intent(out) :: turb
    type(parameter_list), pointer, intent(in) :: params
    logical, optional, intent(in) :: off
    !-
    integer :: stat
    character(:), allocatable :: model
    type(parameter_list), pointer :: turb_params
    logical ::off_

    if (present(off)) then
      off_ = off
    else
      off_ = .false.
    end if

    call params%get('type', model, stat=stat)

    if (stat == 0 .or. off_) then
      call alloc_default_turb_model(turb)
    else
      select case (model)
      case ("alg")
        call alloc_algebraic_turb_model(turb)
      case default
        call alloc_default_turb_model(turb)
      end select
    end if

    turb_params => params%sublist("params")
    call turb%read_params(turb_params)
  end subroutine turbulence_models_read_params


end module turbulence_models
