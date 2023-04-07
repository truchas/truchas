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

    !! Namelist variables
    logical :: use_emfd_solver, use_legacy_bc, graphics_output
    integer :: steps_per_cycle, max_source_cycles, cg_max_iter, output_level
    real(r8) :: matl_change_threshold, steady_state_tol, cg_tol, c_ratio
    character(string_len) :: data_mapper_kind, em_domain_type
    character :: symmetry_axis
    namelist /electromagnetics/ matl_change_threshold, data_mapper_kind, use_emfd_solver, &
      steps_per_cycle, steady_state_tol, max_source_cycles, cg_max_iter, cg_tol, &
      output_level, c_ratio, graphics_output, &
      use_legacy_bc, symmetry_axis, em_domain_type

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

    steps_per_cycle = NULL_I
    max_source_cycles = NULL_I
    steady_state_tol = NULL_R
    cg_max_iter = NULL_I
    cg_tol = NULL_R
    c_ratio = NULL_R

    output_level = NULL_I
    graphics_output = .false.

    use_legacy_bc = .false.
    symmetry_axis  = NULL_C
    em_domain_type = NULL_C

    if (is_IOP) read(lun,nml=electromagnetics,iostat=ios,iomsg=iom)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading ELECTROMAGNETICS namelist: ' // trim(iom))

    call broadcast(matl_change_threshold)
    call broadcast(data_mapper_kind)

    call broadcast(use_emfd_solver)

    call broadcast(steps_per_cycle)
    call broadcast(max_source_cycles)
    call broadcast(steady_state_tol)
    call broadcast(cg_max_iter)
    call broadcast(cg_tol)
    call broadcast(c_ratio)

    call broadcast(output_level)
    call broadcast(graphics_output)

    call broadcast(use_legacy_bc)
    call broadcast(symmetry_axis)
    call broadcast(em_domain_type)

    !! Joule heat driver parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (matl_change_threshold /= NULL_R) then
      if (matl_change_threshold <= 0.0_r8) call TLS_fatal('MATL_CHANGE_THRESHOLD must be > 0.0')
      call params%set('matl-change-threshold', matl_change_threshold)
    end if

    call params%set('frequency-domain-solver', use_emfd_solver)

    if (use_emfd_solver) then
      ! no parameters
      
    else ! use the time domain solver

      !! Time domain Joule heat solver parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

      if (cg_max_iter /= NULL_I) then
        if (cg_max_iter < 1) call TLS_fatal('CG_MAX_ITER must be > 0')
        call params%set('cg-max-iter', cg_max_iter)
      end if

      if (cg_tol /= NULL_R) then
        if (cg_tol <= 0.0_r8 .or. cg_tol >= 0.1_r8) call TLS_fatal('CG_TOL must be > 0.0 and < 0.1')
        call params%set('cg-tol', cg_tol)
      end if

      if (output_level /= NULL_I) call params%set('output-level', output_level)

      call params%set('graphics-output', graphics_output)
    end if

    if (steady_state_tol /= NULL_R) then
      if (steady_state_tol <= 0.0_r8) call TLS_fatal('STEADY_STATE_TOL must be > 0.0')
      call params%set('steady-state-tol', steady_state_tol)
    end if

    if (max_source_cycles /= NULL_I) then
      if (max_source_cycles < 1) call TLS_fatal('MAX_SOURCE_CYCLES must be > 0')
      call params%set('max-source-cycles', max_source_cycles)
    end if

    if (cg_max_iter /= NULL_I) then
      if (cg_max_iter < 1) call TLS_fatal('CG_MAX_ITER must be > 0')
      call params%set('cg-max-iter', cg_max_iter)
    end if

    if (cg_tol /= NULL_R) then
      if (cg_tol <= 0.0_r8 .or. cg_tol >= 0.1_r8) call TLS_fatal('CG_TOL must be > 0.0 and < 0.1')
      call params%set('cg-tol', cg_tol)
    end if

    if (c_ratio /= NULL_R) then
      if (c_ratio <= 0.0_r8 .or. c_ratio > 1.0_r8) call TLS_fatal('C_RATIO must be > 0.0 and <= 1.0')
      call params%set('c-ratio', c_ratio)
    end if

    if (output_level /= NULL_I) call params%set('output-level', output_level)

    call params%set('graphics-output', graphics_output)

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
