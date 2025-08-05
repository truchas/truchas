module induction_heating_namelist

  use parameter_list_type
  use tdme_joule_solver_namelist
  use fdme_solver_namelist
  implicit none
  private

  public :: read_induction_heating_namelist

contains

  subroutine read_induction_heating_namelist(lun, params)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use string_utilities, only: i_to_c
    use input_utilities, only: seek_to_namelist, NULL_C, NULL_R
    use parallel_communication, only: is_IOP, broadcast
    use truchas_logging_services

    integer, intent(in) :: lun
    type(parameter_list), intent(inout) :: params

    integer :: ios
    logical :: found
    character(128) :: iom

    !! Namelist variables
    real(r8) :: prop_change_threshold
    character(32) :: data_mapper_kind
    logical :: use_fd_solver
    namelist /induction_heating/ prop_change_threshold, data_mapper_kind, use_fd_solver

    call TLS_info('Reading INDUCTION_HEATING namelist ...')

    if (is_IOP) rewind(lun)

    if (is_IOP) call seek_to_namelist(lun, 'induction_heating', found, iostat=ios)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
    call broadcast(found)
    if (.not.found) call TLS_fatal('INDUCTION_HEATING namelist not found')

    prop_change_threshold = NULL_R
    data_mapper_kind = NULL_C
    use_fd_solver = .false.

    if (is_IOP) read(lun,nml=induction_heating,iostat=ios,iomsg=iom)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading INDUCTION_HEATING namelist: ' // trim(iom))

    call broadcast(prop_change_threshold)
    call broadcast(data_mapper_kind)
    call broadcast(use_fd_solver)

    if (prop_change_threshold /= NULL_R) then
      if (prop_change_threshold <= 0.0_r8) call TLS_fatal('PROP_CHANGE_THRESHOLD must be > 0.0')
      call params%set('prop-change-threshold', prop_change_threshold)
    end if

    if (data_mapper_kind /= NULL_C) then
      select case (data_mapper_kind)
      case ('default')
      case ('portage')
#ifndef USE_PORTAGE
        call TLS_fatal('DATA_MAPPER_KIND = "portage" is not supported by this Truchas build')
#endif
      case default
        call TLS_fatal('invalid value for DATA_MAPPER_KIND: ' // trim(data_mapper_kind))
      end select
      call params%set('data-mapper-kind', data_mapper_kind)
    end if

    call params%set('use-fd-solver', use_fd_solver)
    if (use_fd_solver) then
      call read_fdme_solver_namelist(lun, params)
    else
      call read_tdme_joule_solver_namelist(lun, params)
    end if

  end subroutine read_induction_heating_namelist

end module induction_heating_namelist
