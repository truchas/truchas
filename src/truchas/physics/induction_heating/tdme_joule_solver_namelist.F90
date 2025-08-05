module tdme_joule_solver_namelist

  use parameter_list_type
  implicit none
  private

  public :: read_tdme_joule_solver_namelist

contains

  subroutine read_tdme_joule_solver_namelist(lun, params)

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
    logical :: graphics_output
    namelist /tdme_joule_solver/ graphics_output

    !! Time-domain method namelist variables
    integer :: steps_per_cycle, max_source_cycles
    real(r8) :: steady_state_tol, c_ratio
    character(32) :: solver_type
    namelist /tdme_joule_solver/ steps_per_cycle, steady_state_tol, max_source_cycles, &
        c_ratio, solver_type

    !! Linear solver variables
    integer :: max_iter, print_level, ams_cycle_type, ams_proj_freq
    real(r8) :: rel_tol
    namelist /tdme_joule_solver/ rel_tol, max_iter, print_level, & ! common
        ams_cycle_type, ams_proj_freq   ! Hypre AMS

    !! Preconditioner variables (TD only)
    integer :: relax_type
    namelist /tdme_joule_solver/ relax_type

    call TLS_info('Reading tdme_joule_solver namelist ...')

    if (is_IOP) rewind(lun)

    if (is_IOP) call seek_to_namelist(lun, 'tdme_joule_solver', found, iostat=ios)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
    call broadcast(found)
    if (.not.found) call TLS_fatal('tdme_joule_solver namelist not found')

    graphics_output = .false.
    steps_per_cycle = NULL_I
    steady_state_tol = NULL_R
    max_source_cycles = NULL_I
    c_ratio = NULL_R
    solver_type = NULL_C
    rel_tol = NULL_R
    max_iter = NULL_I
    print_level = NULL_I
    ams_cycle_type = NULL_I
    ams_proj_freq = NULL_I
    relax_type = NULL_I

    if (is_IOP) read(lun,nml=tdme_joule_solver,iostat=ios,iomsg=iom)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading tdme_joule_solver namelist: ' // trim(iom))

    call broadcast(graphics_output)
    call broadcast(steps_per_cycle)
    call broadcast(steady_state_tol)
    call broadcast(max_source_cycles)
    call broadcast(c_ratio)
    call broadcast(solver_type)
    call broadcast(rel_tol)
    call broadcast(max_iter)
    call broadcast(print_level)
    call broadcast(ams_cycle_type)
    call broadcast(ams_proj_freq)
    call broadcast(relax_type)

    call params%set('graphics-output', graphics_output)

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

    if (solver_type == NULL_C) solver_type = 'pcg'
    select case (solver_type)
    case ('pcg') ! CG with Hiptmair preconditioning
      if (relax_type /= NULL_I) call params%set('relax-type', relax_type)
    case ('ams') ! Hypre AMS solver
      if (ams_cycle_type /= NULL_I) call params%set('ams-cycle-type', ams_cycle_type)
      if (ams_proj_freq /= NULL_I) call params%set('ams-proj-freq', ams_proj_freq)
    case default
      call TLS_fatal('invalid SOLVER_TYPE: ' // solver_type)
    end select
    call params%set('solver-type', solver_type)

    ! used by both solvers
    if (rel_tol /= NULL_R) call params%set('rel-tol', rel_tol)
    if (max_iter /= NULL_I) call params%set('max-iter', max_iter)
    if (print_level /= NULL_I) call params%set('print-level', print_level)

  end subroutine read_tdme_joule_solver_namelist

end module tdme_joule_solver_namelist
