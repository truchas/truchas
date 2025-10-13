module microwave_heating_namelist

  use parameter_list_type
  use fdme_solver_namelist
  implicit none
  private

  public :: read_microwave_heating_namelist

contains

  subroutine read_microwave_heating_namelist(lun, params)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use string_utilities, only: i_to_c
    use input_utilities, only: seek_to_namelist, NULL_C, NULL_R
    use parallel_communication, only: is_IOP, broadcast
    use truchas_logging_services

    integer, intent(in) :: lun
    type(parameter_list), intent(inout) :: params

    integer :: ios, n1, n2
    logical :: found
    character(128) :: iom

    !! Namelist variables
    real(r8) :: prop_change_threshold, times(32), powers(8,33), frequency
    character(32) :: data_mapper_kind, wg_port_bc(8)
    namelist /microwave_heating/ prop_change_threshold, data_mapper_kind, &
        wg_port_bc, times, powers, frequency

    call TLS_info('Reading MICROWAVE_HEATING namelist ...')

    if (is_IOP) rewind(lun)

    if (is_IOP) call seek_to_namelist(lun, 'microwave_heating', found, iostat=ios)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
    call broadcast(found)
    if (.not.found) call TLS_fatal('MICROWAVE_HEATING namelist not found')

    prop_change_threshold = NULL_R
    data_mapper_kind = NULL_C
    wg_port_bc = NULL_C
    times = NULL_R
    powers = NULL_R
    frequency = NULL_R

    if (is_IOP) read(lun,nml=microwave_heating,iostat=ios,iomsg=iom)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading MICROWAVE_HEATING namelist: ' // trim(iom))

    call broadcast(prop_change_threshold)
    call broadcast(data_mapper_kind)
    call broadcast(wg_port_bc)
    call broadcast(times)
    call broadcast(powers)
    call broadcast(frequency)

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

    n1 = findloc(wg_port_bc, NULL_C, dim=1)
    n1 = modulo(n1-1, 1+size(wg_port_bc))
    if (n1 > 0) then
      call params%set('wg-port-bc', wg_port_bc(:n1))

      n2 = findloc(times, NULL_R, dim=1)
      n2 = modulo(n2-1, 1+size(times))
      if (n2 > 0) call params%set('times', times(:n2))

      if (any(powers(:n1,:1+n2) == NULL_R)) call TLS_fatal('POWERS incompletely specified')
      call params%set('powers', powers(:n1,:1+n2))
    else
      call TLS_info('  Questionable! WG_PORT_BC not specified')
    end if

    if (frequency == NULL_R) call TLS_fatal('FREQUENCY not specified')
    call params%set('frequency', frequency)

    call read_fdme_solver_namelist(lun, params)

  end subroutine read_microwave_heating_namelist

end module microwave_heating_namelist
