module turbulence_models
  use turbulence_model_class
  use algebraic_turb_model
  use default_turb_model
  use parameter_list_type
  implicit none
  private

  public :: turbulence_model, turbulence_models_read_params

contains

  subroutine turbulence_models_read_params(turb, params)
    class(turbulence_model), allocatable, intent(out) :: turb
    type(parameter_list), pointer, intent(in) :: params
    !-
    integer :: stat
    character(:), allocatable :: model
    type(parameter_list), pointer :: turb_params

    call params%get('type', model, stat)

    if (stat == 0) then
      call alloc_default_turb_model(turb)
    else
      select case (model)
      case ("alg")
        call alloc_algebraic_turb_moel(turb)
      case default
        call alloc_default_turb_model(turb)
      end select
    end if

    turb_params => params%sublist("params")
    call turb%read_params(turb_params)
  end subroutine turbulence_models_read_params


end module turbulence_models
