module electromagnetics_namelist

  use parameter_list_type
  implicit none
  private
  
  public :: read_electromagnetics_namelist
  
  type(parameter_list), allocatable, public :: params

contains

  subroutine read_electromagnetics_namelist(lun)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use string_utilities, only: i_to_c, raise_case
    use input_utilities, only: seek_to_namelist, NULL_C, NULL_I, NULL_R
    use parallel_communication, only: is_IOP, broadcast
    use truchas_logging_services

    integer, intent(in) :: lun

    integer, parameter :: string_len = 256
    integer, parameter :: MAXSV = 32

    integer :: n, ios
    logical :: found
    character(128) :: iom
    real(r8), allocatable :: src_time(:)
    real(r8), allocatable :: src_freq(:)
    real(r8), allocatable :: unif_src(:)
    real(r8) :: source_times(MAXSV-1)
    real(r8) :: source_frequency(MAXSV)
    real(r8) :: uniform_source(MAXSV)

    !! Domain type: 'FULL_CYLINDER', 'HALF_CYLINDER', or 'QUARTER_CYLINDER'
    character(string_len) :: EM_Domain_Type
    character :: Symmetry_Axis

    !! EM solver control parameters
    integer :: Steps_Per_Cycle
    integer :: Maximum_Source_Cycles
    real(r8) :: SS_Stopping_Tolerance
    integer :: Maximum_CG_Iterations
    real(r8) :: CG_Stopping_Tolerance
    real(r8) :: Num_Etasq
    real(r8) :: Material_Change_Threshold

    !! EM output control parameters
    integer :: Output_Level
    logical :: Graphics_Output
    real(r8) :: Probe_Points(3,10)
    integer :: num_probes

    namelist /electromagnetics/ EM_Domain_Type, Symmetry_Axis, &
      Source_Times, Source_Frequency, Uniform_Source, &
      Steps_Per_Cycle, Maximum_Source_Cycles, SS_Stopping_Tolerance, &
      Maximum_CG_Iterations, CG_Stopping_Tolerance, Material_Change_Threshold, &
      Num_Etasq, Output_Level, Graphics_Output, Probe_Points

    call TLS_info('Reading ELECTROMAGNETICS namelist ...')

    !! Locate the ELECTROMAGNETICS namelist
    if (is_IOP) then
      rewind(lun)
      call seek_to_namelist(lun, 'electromagnetics', found, iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
    call broadcast(found)
    if (.not.found) call TLS_fatal('ELECTROMAGNETICS namelist not found')

    !! Read the namelist
    if (is_IOP) then
      em_domain_type = NULL_C
      symmetry_axis  = NULL_C

      source_times = NULL_R
      source_frequency = NULL_R
      uniform_source = NULL_R

      steps_per_cycle = NULL_I
      maximum_source_cycles = NULL_I
      ss_stopping_tolerance = NULL_R
      maximum_cg_iterations = NULL_I
      cg_stopping_tolerance = NULL_R
      material_change_threshold = NULL_R
      num_etasq = NULL_R

      output_level = NULL_I
      graphics_output = .false.
      probe_points = NULL_R

      read(lun,nml=electromagnetics,iostat=ios,iomsg=iom)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading ELECTROMAGNETICS namelist: ' // trim(iom))

    call broadcast(em_domain_type)
    call broadcast(symmetry_axis)

    call broadcast(source_times)
    call broadcast(source_frequency)
    call broadcast(uniform_source)

    call broadcast(steps_per_cycle)
    call broadcast(maximum_source_cycles)
    call broadcast(ss_stopping_tolerance)
    call broadcast(maximum_cg_iterations)
    call broadcast(cg_stopping_tolerance)
    call broadcast(material_change_threshold)
    call broadcast(num_etasq)

    call broadcast(output_level)
    call broadcast(graphics_output)
    call broadcast(probe_points)
    
    allocate(params)

    em_domain_type = raise_case(trim(adjustl(em_domain_type)))
    select case (em_domain_type)
    case ('FULL_CYLINDER')
    case ('HALF_CYLINDER')
    case ('QUARTER_CYLINDER')
    case ('CYLINDER')
    case ('FRUSTUM')
    case ('VERIFICATION1')
    case (NULL_C)
      call TLS_fatal('EM_DOMAIN_TYPE must be assigned a value')
    case default
      call TLS_fatal('invalid EM_DOMAIN_TYPE: ' // trim(em_domain_type))
    end select
    call params%set('em-domain-type', em_domain_type)

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

    src_time = pack(source_times, mask=(source_times /= NULL_R))
    if (size(src_time) > 1) then
      n = size(src_time)
      if (any(src_time(2:n) <= src_time(:n-1))) &
          call TLS_fatal('SOURCE_TIMES values must be strictly increasing')
    end if
    call params%set('source-times', src_time) !TODO: only if not 0-sized

    src_freq = pack(source_frequency, mask=(source_frequency /= NULL_R))
    if (size(src_freq) == 0) then
      call TLS_fatal('SOURCE_FREQUENCY must be assigned a value')
    else if (any(src_freq <= 0.0_r8)) then
      call TLS_fatal('SOURCE_FREQUENCY values must be > 0.0')
    else if (size(src_freq) /= size(src_time) + 1) then
      call TLS_fatal('wrong number of values provided for SOURCE_FREQUENCY')
    end if
    call params%set('source-frequency', src_freq)

    unif_src = pack(uniform_source, mask=(uniform_source /= NULL_R))
    if (size(unif_src) > 0) then
      if (size(unif_src) /= size(src_freq)) then
        call TLS_fatal('wrong number of values provided for UNIFORM_SOURCE')
      end if
    end if
    call params%set('uiniform-source', unif_src) !TODO: only if not 0-sized

    if (steps_per_cycle == NULL_I) then
      steps_per_cycle = 20
      call TLS_info('  using default value 20 for STEPS_PER_CYCLE')
    else if (steps_per_cycle < 1) then
      call TLS_fatal('STEPS_PER_CYCLE must be > 0')
    else if (steps_per_cycle < 10) then
      call TLS_warn('for decent accuracy, STEPS_PER_CYCLE should be at least 10')
    end if
    call params%set('steps-per-cycle', steps_per_cycle)

    if (maximum_source_cycles == NULL_I) then
      maximum_source_cycles = 10
      call TLS_info('  using default value 10 for MAXIMUM_SOURCE_CYCLES')
    else if (maximum_source_cycles < 1) then
      call TLS_fatal('MAXIMUM_SOURCE_CYCLES must be > 0')
    end if
    call params%set('maximum-source-cycles', maximum_source_cycles) !TODO: optional/default

    if (ss_stopping_tolerance == NULL_R) then
      ss_stopping_tolerance = 1.0e-2_r8
      call TLS_info('  using default value 0.01 for SS_STOPPING_TOLERANCE')
    else if (ss_stopping_tolerance <= 0.0_r8) then
      call TLS_fatal('SS_STOPPING_TOLERANCE must be > 0.0')
    else if (ss_stopping_tolerance > 0.1_r8) then
      call TLS_fatal('SS_STOPPING_TOLERANCE is very loose; consider decreasing the value')
    end if
    call params%set('ss-stopping-tolerance', ss_stopping_tolerance) !TODO: optional/default

    if (maximum_cg_iterations == NULL_I) then
      call TLS_info('  using default value 500 for MAXIMUM_CG_ITERATIONS')
      maximum_cg_iterations = 500
    else if (maximum_cg_iterations < 1) then
      call TLS_fatal('MAXIMUM_CG_ITERATIONS must be > 0')
    end if
    call params%set('maximum-cg-iterations', maximum_cg_iterations) !TODO: optional/default

    if (cg_stopping_tolerance == NULL_R) then
      cg_stopping_tolerance = 1.0e-5_r8
      call TLS_info('  using default value 1.0e-5 for CG_STOPPING_TOLERANCE')
    else if (cg_stopping_tolerance <= 0.0_r8 .or. cg_stopping_tolerance >= 0.1_r8) then
      call TLS_fatal('CG_STOPPING_TOLERANCE must be > 0.0 and < 0.1')
    else if (cg_stopping_tolerance > 1.0e-4_r8) then
      call TLS_warn('CG_STOPPING_TOLERANCE is very loose and may lead to the build up of errors.')
    else if (cg_stopping_tolerance < 1.e-3_r8 * epsilon(1.0_r8)) then
      call TLS_warn('CG_STOPPING_TOLERANCE is too tight; CG iterations are unlikely to converge.')
    end if
    call params%set('cg-stopping-tolerance', cg_stopping_tolerance) !TODO: optional/default

    if (material_change_threshold == NULL_R) then
      material_change_threshold = 0.3_r8
      call TLS_info('  using default value 0.3 for MATERIAL_CHANGE_THRESHOLD')
    else if (material_change_threshold <= 0.0_r8) then
      call TLS_fatal('MATERIAL_CHANGE_THRESHOLD must be > 0.0')
    end if
    call params%set('material-change-threshold', material_change_threshold) !TODO: optional/default

    if (num_etasq == NULL_R) then
      num_etasq = 0.0_r8
    else if (num_etasq < 0.0_r8) then
      call TLS_fatal('NUM_ETASQ must be >= 0.0')
    end if
    call params%set('num-etasq', num_etasq) !TODO: optional/default

    if (output_level == NULL_I) then
      output_level = 1
    else if (output_level < 1 .or. output_level > 4) then
      call TLS_fatal('OUTPUT_LEVEL must be >= 1 and <= 4')
    end if
    call params%set('output-level', output_level) !TODO: optional/default

    call params%set('graphics-output', graphics_output)

    num_probes = 0
    do n = 1, size(probe_points,dim=2)
      if (any(probe_points(:,n) == NULL_R)) exit
      num_probes = n
    end do
    call params%set('num-probes', num_probes)
    call params%set('probe-points', probe_points(:,:num_probes))

  end subroutine read_electromagnetics_namelist

end module electromagnetics_namelist
