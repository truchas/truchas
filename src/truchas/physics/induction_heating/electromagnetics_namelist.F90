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

    integer :: n, ios
    logical :: found
    character(128) :: iom
    type(parameter_list), pointer :: plist

    !! EM heating model namelist variables
    real(r8) :: matl_change_threshold
    character(32) :: data_mapper_kind
    logical :: use_emfd_solver, graphics_output
    namelist /electromagnetics/ matl_change_threshold, data_mapper_kind, use_emfd_solver, &
        graphics_output

    !! Time-domain method namelist variables
    integer :: steps_per_cycle, max_source_cycles
    real(r8) :: steady_state_tol, c_ratio
    character(32) :: td_solver_type
    namelist /electromagnetics/ steps_per_cycle, steady_state_tol, max_source_cycles, &
        c_ratio, td_solver_type

    !! Frequency-domain method namelist variables
    character(32) :: fd_solver_type, fd_precon_type
    namelist /electromagnetics/ fd_solver_type, fd_precon_type

    !! Linear solver variables
    integer :: max_iter, print_level, ams_cycle_type, ams_proj_freq, krylov_dim, max_vec
    real(r8) :: abs_tol, rel_tol, vec_tol
    namelist /electromagnetics/ abs_tol, rel_tol, max_iter, print_level, & ! common
        ams_cycle_type, ams_proj_freq, &   ! Hypre AMS
        krylov_dim, &       ! GMRES
        max_vec, vec_tol    ! NLK

    !! Preconditioner variables (TD only)
    integer :: relax_type
    namelist /electromagnetics/ relax_type

    !! Preconditioner variables (FD only)
    integer  :: max_ams_iter
    namelist /electromagnetics/max_ams_iter

    !! Legacy EM BC variables
    logical :: use_legacy_bc
    character :: symmetry_axis
    character(32) :: em_domain_type
    namelist /electromagnetics/ use_legacy_bc, symmetry_axis, em_domain_type

    call TLS_info('Reading ELECTROMAGNETICS namelist ...')

    if (is_IOP) rewind(lun)

    if (is_IOP) call seek_to_namelist(lun, 'electromagnetics', found, iostat=ios)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
    call broadcast(found)
    if (.not.found) call TLS_fatal('ELECTROMAGNETICS namelist not found')

    matl_change_threshold = NULL_R
    data_mapper_kind = NULL_C
    use_emfd_solver = .false.
    graphics_output = .false.

    steps_per_cycle = NULL_I
    steady_state_tol = NULL_R
    max_source_cycles = NULL_I
    c_ratio = NULL_R
    td_solver_type = NULL_C

    fd_solver_type = NULL_C
    fd_precon_type = NULL_C

    abs_tol = NULL_R
    rel_tol = NULL_R
    max_iter = NULL_I
    print_level = NULL_I
    ams_cycle_type = NULL_I
    ams_proj_freq = NULL_I
    krylov_dim = NULL_I
    max_vec = NULL_I
    vec_tol = NULL_R

    relax_type = NULL_I

    max_ams_iter = NULL_I

    use_legacy_bc = .false.
    symmetry_axis  = NULL_C
    em_domain_type = NULL_C

    if (is_IOP) read(lun,nml=electromagnetics,iostat=ios,iomsg=iom)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading ELECTROMAGNETICS namelist: ' // trim(iom))

    call broadcast(matl_change_threshold)
    call broadcast(data_mapper_kind)
    call broadcast(use_emfd_solver)
    call broadcast(graphics_output)

    call broadcast(steps_per_cycle)
    call broadcast(steady_state_tol)
    call broadcast(max_source_cycles)
    call broadcast(c_ratio)
    call broadcast(td_solver_type)

    call broadcast(fd_solver_type)
    call broadcast(fd_precon_type)

    call broadcast(abs_tol)
    call broadcast(rel_tol)
    call broadcast(max_iter)
    call broadcast(print_level)
    call broadcast(ams_cycle_type)
    call broadcast(ams_proj_freq)
    call broadcast(krylov_dim)
    call broadcast(max_vec)
    call broadcast(vec_tol)

    call broadcast(relax_type)

    call broadcast(max_ams_iter)

    call broadcast(use_legacy_bc)
    call broadcast(symmetry_axis)
    call broadcast(em_domain_type)

    !! Joule heat driver parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (matl_change_threshold /= NULL_R) then
      if (matl_change_threshold <= 0.0_r8) call TLS_fatal('MATL_CHANGE_THRESHOLD must be > 0.0')
      call params%set('matl-change-threshold', matl_change_threshold)
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

    call params%set('frequency-domain-solver', use_emfd_solver)

    if (use_emfd_solver) then ! Frequency-domain solver parameters

      plist => params%sublist('emfd-solver')
      call plist%set('graphics-output', graphics_output)

      if (abs_tol /= NULL_R) call plist%set('abs-tol', abs_tol)
      if (rel_tol /= NULL_R) call plist%set('rel-tol', rel_tol)
      if (max_iter /= NULL_I) call plist%set('max-iter', max_iter)
      if (print_level /= NULL_I) call plist%set('print-level', print_level)

      select case (fd_solver_type)
      case ('nlk')
        if (max_vec /= NULL_I) call plist%set('max-vec', max_vec)
        if (vec_tol /= NULL_R) call plist%set('vec-tol', vec_tol)
      case ('gmres')
        if (krylov_dim /= NULL_I) call plist%set('krylov-dim', krylov_dim)
      case ('minres')
      case (NULL_C)
        call TLS_fatal('FD_SOLVER_TYPE not specified')
      case default
        call TLS_fatal('invalid FD_SOLVER_TYPE: ' // fd_solver_type)
      end select
      call plist%set('solver-type', fd_solver_type)

      plist => plist%sublist('precon')

      select case (fd_precon_type)
      case ('ams')
        if (max_ams_iter /= NULL_I) call plist%set('max-iter', max_ams_iter)
        if (ams_cycle_type /= NULL_I) call plist%set('ams-cycle-type', ams_cycle_type)
      case ('hiptmair')
      case (NULL_C)
        call TLS_fatal('FD_PRECON_TYPE not specified')
      case default
        call TLS_fatal('invalid FD_PRECON_TYPE: ' // fd_precon_type)
      end select
      call plist%set('type', fd_precon_type)

    else ! Time domain solver parameters

      call params%set('graphics-output', graphics_output)

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

    !! Parameters for the legacy method of inferring boundary conditions !!!!!!!

    call params%set('use-legacy-bc', use_legacy_bc)
    if (use_legacy_bc) then
      symmetry_axis = raise_case(trim(adjustl(symmetry_axis)))
      select case (symmetry_axis)
      case ('X', 'Y', 'Z')
      case (NULL_C)
        call TLS_info('  using default value "Z" for SYMMETRY_AXIS')
        symmetry_axis = 'Z'
      case default
        call TLS_fatal('invalid SYMMETRY_AXIS: ' // trim(symmetry_axis))
      end select
      call params%set('symmetry-axis', symmetry_axis)

      em_domain_type = raise_case(trim(adjustl(em_domain_type)))
      select case (em_domain_type)
      case ('FULL_CYLINDER')
      case ('HALF_CYLINDER')
      case ('QUARTER_CYLINDER')
      case ('CYLINDER')
      case ('FRUSTUM')
      case ('VERIFICATION1')
      case (NULL_C)
        call TLS_fatal('EM_DOMAIN_TYPE must be assigned a value when USE_LEGACY_BC is true')
      case default
        call TLS_fatal('invalid EM_DOMAIN_TYPE: ' // trim(em_domain_type))
      end select
      call params%set('em-domain-type', em_domain_type)
    end if

  end subroutine read_electromagnetics_namelist

end module electromagnetics_namelist
