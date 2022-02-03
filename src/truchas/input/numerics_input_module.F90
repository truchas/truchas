!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module numerics_input_module

  use kinds, only: r8
  use truchas_logging_services
  implicit none
  private

  public :: numerics_input

contains

  subroutine numerics_input(lun)

    use parallel_communication, only: is_IOP, broadcast
    use string_utilities,       only: i_to_c
    use input_utilities,        only: seek_to_namelist
    use cutoffs_module,         only: alittle, cutvof
    use time_step_module,       only: t, dt_init, dt_constant, dt_grow, dt_min, dt_max, &
                                      cycle_max, cycle_number

    integer, intent(in) :: lun

    integer :: ios
    logical :: found
    character(80) :: iom

    namelist /numerics/ alittle, cutvof, cycle_number, cycle_max, &
        t, dt_constant, dt_init, dt_grow, dt_min, dt_max

    call TLS_info('Reading NUMERICS namelist ...')

    !! Locate the NUMERICS namelist (required)
    if (is_IOP) then
      rewind(lun)
      call seek_to_namelist(lun, 'numerics', found, iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
    call broadcast(found)
    if (.not.found) call TLS_fatal('NUMERICS namelist not found')

    !! Set default values for namelist variables, and others
    call numerics_default

    !! Read the NUMERICS namelist
    if (is_IOP) read(lun,nml=numerics,iostat=ios,iomsg=iom)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading NUMERICS namelist: ' // trim(iom))

    !! Broadcast namelist variables
    call numerics_input_parallel

    !! Check input for errors, and other intialization
    call numerics_check

  end subroutine numerics_input


  subroutine numerics_check

    use input_utilities,   only: NULL_R
    use string_utilities,  only: lower_case
    use cutoffs_module,    only: cutvof
    use time_step_module,  only: t, dt, constant_dt, dt_constraint, &
                                 dt_init, dt_constant, dt_grow, dt_min, dt_max, &
                                 cycle_max, cycle_number

    if (cycle_max < 0) call TLS_fatal('CYCLE_MAX must be >= 0')
    if (cycle_number < 0) call TLS_fatal('CYCLE_NUMBER must be >= 0')
    if (cutvof < 0 .or. cutvof > 1) call TLS_fatal('CUTVOF must be in [0,1]')
    if (dt_grow < 1) call TLS_fatal('DT_GROW must be >= 1')
    if (dt_constant <= 0 .and. dt_constant /= NULL_R) call TLS_fatal('DT_CONSTANT must be > 0.0')
    if (dt_init <= 0) call TLS_fatal('DT_INIT must be > 0.0')
    if (dt_max <= 0) call TLS_fatal('DT_MAX must be > 0.0')
    if (dt_min <= 0) call TLS_fatal('DT_MIN must be > 0.0')

    ! Set initial time step.
    if (dt_constant /= NULL_R) then
      constant_dt = .true.
      dt = dt_constant
      dt_constraint = 'constant'
    else
      dt_constant = 0
      dt = dt_init
      dt_constraint = 'initial'
    end if

  end subroutine numerics_check


  subroutine numerics_default

    use input_utilities, only: NULL_R
    use cutoffs_module, only: cutvof
    use time_step_module, only: t, constant_dt, dt_constant, dt_grow, dt_init, dt_max, dt_min, &
        cycle_max, cycle_number, cycle_number_restart

    t = 0
    cycle_number = 0
    cycle_number_restart = 0
    cycle_max = 1000000

    dt_init = 1.0d-6
    dt_min  = 1.0d-6
    dt_max  = 10.0
    dt_grow = 1.05

    dt_constant = NULL_R
    constant_dt = .false.
    cutvof  = 1.0d-8 ! volume fraction cutoff

  end subroutine numerics_default


  subroutine numerics_input_parallel

    use cutoffs_module,         only: alittle, cutvof
    use parallel_communication, only: broadcast
    use time_step_module,       only: t, constant_dt, cycle_max, &
                                      cycle_number, cycle_number_restart,     &
                                      dt_constant, dt_constraint, dt_grow,    &
                                      dt_init, dt_max, dt_min

    call broadcast(alittle)
    call broadcast(cutvof)

    call broadcast(t)
    call broadcast(dt_constant)
    call broadcast(dt_init)
    call broadcast(dt_grow)
    call broadcast(dt_min)
    call broadcast(dt_max)
    call broadcast(cycle_number)
    call broadcast(cycle_max)

    ! Non-namelist variables.
    call broadcast(constant_dt)
    call broadcast(cycle_number_restart)
    call broadcast(dt_constraint)

  end subroutine numerics_input_parallel

end module numerics_input_module
