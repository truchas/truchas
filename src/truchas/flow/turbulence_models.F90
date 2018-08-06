!!
!! TURBULENCE_MODELS
!!
!! Include this module to use any of the turbulence models provided by truchas.
!!
!! Peter Brady <ptb@lanl.gov>
!! 2018
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PROGRAMMING INTERFACE
!!
!!  READ_FLOW_TURBULENCE_MODEL_NAMELIST (LUN, PARAMETER_LIST)
!!    reads namelists /flow_turbulence_model/ and stores into parameter list
!!
!!  TURBULENCE_MODELS_READ_PARAMS(model, PARAMETER_LIST)
!!    Allocates a turbulence_model based on parameter_list
!!

module turbulence_models
  use kinds
  use truchas_logging_services
  use turbulence_model_class
  use algebraic_turb_model_type
  use default_turb_model_type
  use parameter_list_type
  implicit none
  private

  public :: turbulence_model, turbulence_models_read_params, &
      read_flow_turbulence_model_namelist

contains

  subroutine read_flow_turbulence_model_namelist(lun, p)
    use string_utilities, only: i_to_c
    use parallel_communication, only: is_IOP, broadcast
    use flow_input_utils

    integer, intent(in) :: lun
    type(parameter_list), pointer, intent(inout) :: p
    type(parameter_list), pointer :: pp
    integer :: ios
    logical :: found
    character(128) :: iom
    character(128) :: type
    real(r8) :: length, cmu, ke_fraction

    namelist /flow_turbulence_model/ type, length, cmu, ke_fraction

    type = null_c
    length = null_r
    cmu = null_r
    ke_fraction = null_r

    pp => p%sublist("turbulence model")

    if (is_IOP) then
      rewind lun
      call seek_to_namelist(lun, 'FLOW_TURBULENCE_MODEL', found, iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('Error reading input file: iostat=' // i_to_c(ios))

    call broadcast(found)
    if (found) then
      call TLS_info('')
      call TLS_info('Reading FLOW_TURBULENCE_MODEL namelist ...')
      !! Read the namelist.
      if (is_IOP) then
        read(lun,nml=flow_turbulence_model,iostat=ios,iomsg=iom)
      end if
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading FLOW_TURBULENCE_MODEL namelist: ' // trim(iom))

      call broadcast(type)

      if (type /= null_c) then
        call broadcast(length)
        call broadcast(cmu)
        call broadcast(ke_fraction)

        call pp%set("type", trim(type))
        pp => pp%sublist("params")
        call plist_set_if(pp, "length", length)
        call plist_set_if(pp, "cmu", cmu)
        call plist_set_if(pp, "ke fraction", ke_fraction)
      end if
    end if
  end subroutine read_flow_turbulence_model_namelist

  subroutine turbulence_models_read_params(turb, params, off)
    class(turbulence_model), allocatable, intent(out) :: turb
    type(parameter_list), pointer, intent(in) :: params
    logical, optional, intent(in) :: off
    !-
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
  end subroutine turbulence_models_read_params


end module turbulence_models
