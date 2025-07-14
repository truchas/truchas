module electromagnetics_namelist

  use parameter_list_type
  implicit none
  private

  public :: read_electromagnetics_namelist

contains

  subroutine read_electromagnetics_namelist(lun, params)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use string_utilities, only: i_to_c, raise_case
    use input_utilities, only: seek_to_namelist, NULL_C, NULL_I, NULL_R
    use parallel_communication, only: is_IOP, broadcast
    use truchas_logging_services

    integer, intent(in) :: lun
    type(parameter_list), intent(inout) :: params

    integer, parameter :: string_len = 256

    integer :: n, ios, n1, n2
    logical :: found
    character(128) :: iom
    type(parameter_list), pointer :: plist, sublist

    !! EM heating model namelist variables
    real(r8) :: prop_change_threshold
    character(32) :: data_mapper_kind
    logical :: use_fd_solver, graphics_output, use_mixed_form, use_mw_solver
    namelist /electromagnetics/ prop_change_threshold, data_mapper_kind, use_fd_solver, &
        use_mixed_form, graphics_output, use_mw_solver

    !! Time-domain method namelist variables
    integer :: steps_per_cycle, max_source_cycles
    real(r8) :: steady_state_tol, c_ratio
    character(32) :: td_solver_type
    namelist /electromagnetics/ steps_per_cycle, steady_state_tol, max_source_cycles, &
        c_ratio, td_solver_type

    !! Frequency-domain method namelist variables
    real(r8) :: minres_pc_damping
    character(32) :: fd_solver_type, fd_precon_type
    namelist /electromagnetics/ fd_solver_type, fd_precon_type, minres_pc_damping

    !! Linear solver variables
    integer :: max_iter, print_level, ams_cycle_type, ams_proj_freq
    real(r8) :: abs_tol, rel_tol
    namelist /electromagnetics/ abs_tol, rel_tol, max_iter, print_level, & ! common
        ams_cycle_type, ams_proj_freq   ! Hypre AMS

    !! Preconditioner variables (TD only)
    integer :: relax_type
    namelist /electromagnetics/ relax_type

    !! HYPRE BoomerAMG preconditioner variables (FD only)
    integer  :: boomer_num_cycles
    real(r8) :: boomer_strong_threshold
    integer  :: boomer_coarsen_type, boomer_interp_type
    integer  :: boomer_relax_down_type, boomer_relax_up_type
    integer  :: boomer_print_level, boomer_debug_level, boomer_logging_level
    namelist /electromagnetics/ boomer_num_cycles, boomer_strong_threshold, boomer_coarsen_type, &
        boomer_interp_type, boomer_relax_down_type, boomer_relax_up_type, boomer_print_level, &
        boomer_debug_level, boomer_logging_level

    !! Built-in SSOR preconditioner
    integer  :: ssor_num_cycles
    real(r8) :: ssor_omega
    namelist /electromagnetics/ ssor_num_cycles, ssor_omega

    real(r8) :: times(32), powers(8,33), frequency
    character(32) :: wg_port_bc(8)
    namelist /electromagnetics/ wg_port_bc, times, powers, frequency

    call TLS_info('Reading ELECTROMAGNETICS namelist ...')

    if (is_IOP) rewind(lun)

    if (is_IOP) call seek_to_namelist(lun, 'electromagnetics', found, iostat=ios)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
    call broadcast(found)
    if (.not.found) call TLS_fatal('ELECTROMAGNETICS namelist not found')

    prop_change_threshold = NULL_R
    data_mapper_kind = NULL_C
    use_fd_solver = .false.
    use_mixed_form = .false.
    graphics_output = .false.
    use_mw_solver = .false.

    steps_per_cycle = NULL_I
    steady_state_tol = NULL_R
    max_source_cycles = NULL_I
    c_ratio = NULL_R
    td_solver_type = NULL_C

    fd_solver_type = NULL_C
    fd_precon_type = NULL_C
    minres_pc_damping = NULL_R

    abs_tol = NULL_R
    rel_tol = NULL_R
    max_iter = NULL_I
    print_level = NULL_I
    ams_cycle_type = NULL_I
    ams_proj_freq = NULL_I

    relax_type = NULL_I

    boomer_num_cycles = NULL_I
    boomer_strong_threshold = NULL_R
    boomer_coarsen_type = NULL_I
    boomer_interp_type = NULL_I
    boomer_relax_down_type = NULL_I
    boomer_relax_up_type = NULL_I
    boomer_print_level = NULL_I
    boomer_debug_level = NULL_I
    boomer_logging_level = NULL_I
    ssor_num_cycles = NULL_I
    ssor_omega = NULL_R

    wg_port_bc = NULL_C
    times = NULL_R
    powers = NULL_R
    frequency = NULL_R

    if (is_IOP) read(lun,nml=electromagnetics,iostat=ios,iomsg=iom)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading ELECTROMAGNETICS namelist: ' // trim(iom))

    call broadcast(prop_change_threshold)
    call broadcast(data_mapper_kind)
    call broadcast(use_fd_solver)
    call broadcast(use_mixed_form)
    call broadcast(graphics_output)
    call broadcast(use_mw_solver)

    call broadcast(steps_per_cycle)
    call broadcast(steady_state_tol)
    call broadcast(max_source_cycles)
    call broadcast(c_ratio)
    call broadcast(td_solver_type)

    call broadcast(fd_solver_type)
    call broadcast(fd_precon_type)
    call broadcast(minres_pc_damping)

    call broadcast(abs_tol)
    call broadcast(rel_tol)
    call broadcast(max_iter)
    call broadcast(print_level)
    call broadcast(ams_cycle_type)
    call broadcast(ams_proj_freq)

    call broadcast(relax_type)

    call broadcast(boomer_num_cycles)
    call broadcast(boomer_strong_threshold)
    call broadcast(boomer_coarsen_type)
    call broadcast(boomer_interp_type)
    call broadcast(boomer_relax_down_type)
    call broadcast(boomer_relax_up_type)
    call broadcast(boomer_print_level)
    call broadcast(boomer_debug_level)
    call broadcast(boomer_logging_level)
    call broadcast(ssor_num_cycles)
    call broadcast(ssor_omega)

    call broadcast(wg_port_bc)
    call broadcast(times)
    call broadcast(powers)
    call broadcast(frequency)

    !! Joule heat driver parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    call params%set('use-mw-solver', use_mw_solver)
    call params%set('graphics-output', graphics_output)

    if (use_fd_solver .or. use_mw_solver) then ! Frequency-domain solver parameters

      if (use_mw_solver) then
        n1 = findloc(wg_port_bc, NULL_C, dim=1)
        n1 = modulo(n1-1, 1+size(wg_port_bc))
        if (n1 == 0) call TLS_fatal('WG_PORT_BC not specified')
        call params%set('wg-port-bc', wg_port_bc(:n1))
        
        n2 = findloc(times, NULL_R, dim=1)
        n2 = modulo(n2-1, 1+size(times))
        if (n2 > 0) call params%set('times', times(:n2))

        if (any(powers(:n1,:1+n2) == NULL_R)) call TLS_fatal('POWERS incompletely specified')
        call params%set('powers', powers(:n1,:1+n2))
        
        if (frequency == NULL_R) call TLS_fatal('FREQUENCY not specified')
        call params%set('frequency', frequency)
      end if

      call params%set('use-mixed-form', use_mixed_form)

      if (abs_tol /= NULL_R) call params%set('abs-tol', abs_tol)
      if (rel_tol /= NULL_R) call params%set('rel-tol', rel_tol)
      if (max_iter /= NULL_I) call params%set('max-iter', max_iter)
      if (print_level /= NULL_I) call params%set('print-level', print_level)

      select case (fd_solver_type)
      case ('minres')
      case ('mumps')
#ifndef USE_MUMPS
        call TLS_fatal('FD_SOLVER_TYPE = "mumps" is not supported by this Truchas build')
#endif
      case (NULL_C)
        call TLS_fatal('FD_SOLVER_TYPE not specified')
      case default
        call TLS_fatal('invalid FD_SOLVER_TYPE: ' // fd_solver_type)
      end select
      call params%set('solver-type', fd_solver_type)

      if (fd_solver_type /= 'mumps') then
        plist => params%sublist('precon')
        select case (fd_precon_type)
        case ('none')
        case ('boomer')
          sublist => plist%sublist('params')
          if (boomer_num_cycles /= NULL_I)  call sublist%set('num-cycles', boomer_num_cycles)
          if (boomer_strong_threshold /= NULL_R)  call sublist%set('strong-threshold', boomer_strong_threshold)
          if (boomer_coarsen_type /= NULL_I)  call sublist%set('coarsen-type', boomer_coarsen_type)
          if (boomer_interp_type /= NULL_I)  call sublist%set('interp-type', boomer_interp_type)
          if (boomer_relax_down_type /= NULL_I)  call sublist%set('relax-down-type', boomer_relax_down_type)
          if (boomer_relax_up_type /= NULL_I)  call sublist%set('relax-up-type', boomer_relax_up_type)
          if (boomer_print_level /= NULL_I)  call sublist%set('print-level', boomer_print_level)
          if (boomer_debug_level /= NULL_I)  call sublist%set('debug-level', boomer_debug_level)
          if (boomer_logging_level /= NULL_I)  call sublist%set('logging-level', boomer_logging_level)
        case ('ssor')
          sublist => plist%sublist('params')
          if (ssor_omega /= NULL_R)  call sublist%set('omega', ssor_omega)
          if (ssor_num_cycles /= NULL_I) call sublist%set('num-cycles', ssor_num_cycles)
        case (NULL_C)
          call TLS_fatal('FD_PRECON_TYPE not specified')
        case default
          call TLS_fatal('invalid FD_PRECON_TYPE: ' // fd_precon_type)
        end select
        call plist%set('type', fd_precon_type)
        call plist%set('method', fd_precon_type) ! for pcsr_precon_factory (FIXME)
        select case (fd_solver_type)
        case ('minres')
          if (minres_pc_damping /= NULL_R) call plist%set('beta', minres_pc_damping)
        end select
      end if

    else ! Time domain solver parameters

      !! All these parameters have default values which are applied, if needed,
      !! when the parameter list is consumed.
      if (steps_per_cycle /= NULL_I) then
        if (steps_per_cycle < 1) call TLS_fatal('STEPS_PER_CYCLE must be > 0')
        call params%set('steps-per-cycle', steps_per_cycle)
      end if

      if (steady_state_tol /= NULL_R) then
        if (steady_state_tol <= 0.0_r8) call TLS_fatal('STEADY_STATE_TOL must be > 0.0')
        call params%set('steady-state-tol', steady_state_tol)
      end if

      if (max_source_cycles /= NULL_I) then
        if (max_source_cycles < 1) call TLS_fatal('MAX_SOURCE_CYCLES must be > 0')
        call params%set('max-source-cycles', max_source_cycles)
      end if

      if (c_ratio /= NULL_R) then
        if (c_ratio <= 0.0_r8 .or. c_ratio > 1.0_r8) call TLS_fatal('C_RATIO must be > 0.0 and <= 1.0')
        call params%set('c-ratio', c_ratio)
      end if

      if (rel_tol /= NULL_R) call params%set('rel-tol', rel_tol)
      if (max_iter /= NULL_I) call params%set('max-iter', max_iter)
      if (print_level /= NULL_I) call params%set('print-level', print_level)

      if (td_solver_type == NULL_C) td_solver_type = 'pcg'
      select case (td_solver_type)
      case ('pcg') ! CG with Hiptmair preconditioning
        if (relax_type /= NULL_I) call params%set('relax-type', relax_type)
      case ('ams') ! Hypre AMS solver
        if (ams_cycle_type /= NULL_I) call params%set('ams-cycle-type', ams_cycle_type)
        if (ams_proj_freq /= NULL_I) call params%set('ams-proj-freq', ams_proj_freq)
      case default
        call TLS_fatal('invalid TD_SOLVER_TYPE: ' // td_solver_type)
      end select
      call params%set('td-solver-type', td_solver_type)

    end if

  end subroutine read_electromagnetics_namelist

end module electromagnetics_namelist
