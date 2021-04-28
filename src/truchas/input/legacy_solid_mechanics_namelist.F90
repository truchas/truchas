module legacy_solid_mechanics_namelist

  use parameter_list_type
  use solid_mechanics_input
  implicit none
  private

  public :: read_legacy_solid_mechanics_namelist

  type(parameter_list), public :: nonlinear_params

contains

  subroutine read_legacy_solid_mechanics_namelist(lun)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use parallel_communication, only: is_IOP, broadcast
    use input_utilities, only: seek_to_namelist, NULL_R, NULL_I
    use string_utilities, only: i_to_c
    use truchas_logging_services

    integer, intent(in) :: lun

    namelist /legacy_solid_mechanics/ solid_mechanics_body_force, stress_reduced_integration, &
        contact_distance, contact_norm_trac, contact_penalty, strain_limit

    !! Parameters formerly in LINEAR_SOLVER and NONLINEAR_SOLVER
    real(r8) :: convergence_criterion, nlk_vector_tolerance
    integer :: maximum_iterations, nlk_max_vectors
    namelist /legacy_solid_mechanics/convergence_criterion, maximum_iterations, nlk_vector_tolerance, &
        nlk_max_vectors

    integer :: ios
    logical :: found
    character(80) :: iom

    call TLS_info('Reading LEGACY_SOLID_MECHANICS namelist ...')

    !! Locate the LEGACY_SOLID_MECHANICS namelist (required)
    if (is_IOP) then
      rewind(lun)
      call seek_to_namelist(lun, 'legacy_solid_mechanics', found, iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
    call broadcast(found)
    if (.not.found) call TLS_fatal('LEGACY_SOLID_MECHANICS namelist not found')

    !! Default values
    solid_mechanics_body_force = .false.
    stress_reduced_integration = .false.
    contact_distance  = 1.0e-7_r8
    contact_norm_trac = 1.0e4_r8
    contact_penalty   = 1.0e3_r8
    strain_limit = 1.0e-10_r8

    convergence_criterion = NULL_R
    maximum_iterations = NULL_I
    nlk_vector_tolerance = NULL_R
    nlk_max_vectors = NULL_I

    !! Read the namelist
    if (is_IOP) read(lun,nml=legacy_solid_mechanics,iostat=ios,iomsg=iom)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading LEGACY_SOLID_MECHANICS namelist: ' // trim(iom))

    !! Broadcast the namelist variables
    call broadcast(solid_mechanics_body_force)
    call broadcast(stress_reduced_integration)
    call broadcast(contact_distance)
    call broadcast(contact_norm_trac)
    call broadcast(contact_penalty)
    call broadcast(strain_limit)

    call broadcast(maximum_iterations)
    call broadcast(convergence_criterion)
    call broadcast(nlk_max_vectors)
    call broadcast(nlk_vector_tolerance)

    if (strain_limit < 0) call TLS_fatal('STRAIN_LIMIT must be >= 0.0')
    ! none of the other real parameters were checked :-/

    if (maximum_iterations /= NULL_I) call nonlinear_params%set('nlk-max-iter', maximum_iterations)
    if (convergence_criterion /= NULL_R) call nonlinear_params%set('nlk-tol', convergence_criterion)
    if (nlk_max_vectors /= NULL_I) call nonlinear_params%set('nlk-max-vec', nlk_max_vectors)
    if (nlk_vector_tolerance /= NULL_R) call nonlinear_params%set('nlk-vec-tol', nlk_vector_tolerance)

  end subroutine read_legacy_solid_mechanics_namelist

end module legacy_solid_mechanics_namelist
