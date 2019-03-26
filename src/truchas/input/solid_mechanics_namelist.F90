module solid_mechanics_namelist

  use solid_mechanics_input
  implicit none
  private

  public :: read_solid_mechanics_namelist

contains

  subroutine read_solid_mechanics_namelist(lun)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use parallel_communication, only: is_IOP, broadcast
    use input_utilities, only: seek_to_namelist
    use string_utilities, only: i_to_c
    use truchas_logging_services

    integer, intent(in) :: lun

    namelist /solid_mechanics/ solid_mechanics_body_force, stress_reduced_integration, &
        displacement_linear_solution, displacement_nonlinear_solution, &
        contact_distance, contact_norm_trac, contact_penalty, strain_limit

    integer :: ios
    logical :: found
    character(80) :: iom

    call TLS_info('')
    call TLS_info('Reading SOLID_MECHANICS namelist ...')

    !! Locate the SOLID_MECHANICS namelist (required)
    if (is_IOP) then
      rewind(lun)
      call seek_to_namelist(lun, 'solid_mechanics', found, iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
    call broadcast(found)
    if (.not.found) call TLS_fatal('SOLID_MECHANICS namelist not found')

    !! Default values
    solid_mechanics_body_force = .false.
    stress_reduced_integration = .false.
    displacement_linear_solution = 'default'
    displacement_nonlinear_solution = 'default'
    contact_distance  = 1.0e-7_r8
    contact_norm_trac = 1.0e4_r8
    contact_penalty   = 1.0e3_r8
    strain_limit = 1.0e-10_r8

    !! Read the FLOW namelist
    if (is_IOP) read(lun,nml=solid_mechanics,iostat=ios,iomsg=iom)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading SOLID_MECHANICS namelist: ' // trim(iom))

    !! Broadcast the namelist variables
    call broadcast(solid_mechanics_body_force)
    call broadcast(stress_reduced_integration)
    call broadcast(displacement_linear_solution)
    call broadcast(displacement_nonlinear_solution)
    call broadcast(contact_distance)
    call broadcast(contact_norm_trac)
    call broadcast(contact_penalty)
    call broadcast(strain_limit)

    if (strain_limit < 0) call TLS_fatal('STRAIN_LIMIT must be >= 0.0')
    ! none of the other real parameters were checked :-/

    call solver_init

  end subroutine read_solid_mechanics_namelist

  subroutine solver_init

    use linear_solution, only: UBIK_NK_DEFAULT, linear_solutions, Ubik_user, DEFAULT_UBIK_CONTROLS
    use nonlinear_solution, only: NK_DEFAULT, NKuser, nonlinear_solutions, DEFAULT_NK_CONTROLS
    use parameter_module, only: string_len
    use utilities_module, only: string_compare
    use truchas_logging_services

    integer :: j, k
    character(string_len) :: string
    logical :: this_string_matches

    NK_DISPLACEMENT = NK_DEFAULT

    do j = DEFAULT_NK_CONTROLS + 1, DEFAULT_NK_CONTROLS + nonlinear_solutions
      string = NKuser(j)%name
      this_string_matches = .false.
      call STRING_COMPARE (TRIM(displacement_nonlinear_solution), string, this_string_matches)
      if (this_string_matches) then
        NK_DISPLACEMENT = j
        call TLS_info('Using nonlinear solver "' // trim(displacement_nonlinear_solution) // &
                      '" for displacement nonlinear solution.')
        exit
      end if
    end do

    if (NK_DISPLACEMENT == NK_DEFAULT) then
      if (displacement_nonlinear_solution == 'default') then
        call TLS_info ('Using default nonlinear solver parameters for displacement nonlinear solution.')
      else
        call TLS_warn('Nonlinear solver "' // trim(displacement_nonlinear_solution) // &
                      '" for displacement nonlinear solution not found! &
                      & Reverting to default nonlinear solver parameters.')
      end if
    end if

    ! Find the correct linear solver for each input nonlinear solver.
    do j = DEFAULT_NK_CONTROLS + 1, DEFAULT_NK_CONTROLS + nonlinear_solutions
      string = NKuser(j)%linear_solver_name
      this_string_matches = .false.
      do k = DEFAULT_UBIK_CONTROLS + 1, DEFAULT_UBIK_CONTROLS + linear_solutions
        call string_compare(TRIM(Ubik_user(k)%name), string, this_string_matches)
        if (this_string_matches) then
          NKuser(j)%linear_solver_index = k
          call TLS_info('Using linear solver "' // trim(Ubik_user(k)%name) // &
                        '" for nonlinear solver " // trim(NKuser(j)%name) // ".')
          exit
        end if
        if (NKuser(j)%linear_solver_index == UBIK_NK_DEFAULT) then
          call TLS_info('Using default linear solver "' // trim(Ubik_user(k)%name) // &
                        '" for nonlinear solver "' // trim(NKuser(j)%name) // '".')
        end if
      end do
    end do

  end subroutine solver_init

end module solid_mechanics_namelist
