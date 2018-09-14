!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE NUMERICS_INPUT_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define various procedures for the input of global numerical
  !   algorithm flags and parameters.
  !
  !   Public Interface:
  !
  !     * call NUMERICS_INPUT ()
  !
  !       Defaults, reads, checks, and broadcasts input variables
  !       in the NUMERICS namelist.
  !
  ! Contains: NUMERICS_INPUT
  !           NUMERICS_CHECK
  !           NUMERICS_DEFAULT
  !           NUMERICS_INPUT_PARALLEL
  !
  ! Author(s): The Telluridians (telluride-info@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use truchas_logging_services
  implicit none
  private

  public :: NUMERICS_INPUT
  
CONTAINS

  SUBROUTINE NUMERICS_INPUT (lun)
    !=======================================================================
    ! Purpose(s):
    !
    !   Read NUMERICS namelist.
    !
    !=======================================================================
    use cutoffs_module,         only: alittle, cutvof
    use discrete_ops_data,      only: discrete_ops_type
    use ff_discrete_ops_data,   only: ff_discrete_ops_type
    use body_data_module,       only: body_force_implicitness, &
                                      body_force_face_method
    use fluid_data_module,      only: MinFaceFraction,  &
                                      momentum_solidify_implicitness, &
                                      mass_limiter, mass_limiter_cutoff
    use flux_volume_module,     only: flux_vol_iter_max
    use input_utilities,        only: seek_to_namelist
    use interface_module,       only: interface_topology_model
    use parallel_info_module,   only: p_info
    use porous_drag_data,       only: porous_implicitness
    use projection_data_module, only: projection_linear_solution
    use mollify,                only: interface_smoothing_length
    use time_step_module,       only: courant_number, cycle_max, cycle_number, &
                                      dt_constant, dt_grow, dt_init, dt_max,  &
                                      dt_min, t, viscous_number,&
                                      surften_number
    use viscous_data_module,    only: viscous_linear_solution,    &
                                      viscous_implicitness

    use vof_data_module,        only: volume_track_brents_method, &
                                      volume_track_iter_max,      &
                                      volume_track_subcycles,     &
                                      volume_track_iter_tol,      &
                                      volume_track_interfaces,    &
                                      interface_geometry,         &
                                      interface_area
    use solid_mechanics_input,  only: displacement_linear_solution,    &
                                      displacement_nonlinear_solution, &
                                      contact_distance,                &
                                      contact_norm_trac,               &
                                      contact_penalty, strain_limit,   &
                                      stress_reduced_integration
    use advection_data,         only: limiter_type,               &
                                      advection_order_vol,        &  
                                      advection_order_energy,     &
                                      advection_order_momentum,   &
                                      advection_order_species
    use body_data_module,       only: mechanical_energy_bound

    integer, intent(in) :: lun

    ! Argument List

    ! Local Variables
    character(128) :: message
    logical :: fatal, found
    integer :: ioerror

    ! Define NUMERICS namelist.
    namelist /NUMERICS/                                          &
         alittle, cutvof,                                        &

         courant_number,                                         &
         momentum_solidify_implicitness,                         &
         body_force_implicitness,                                &
         body_force_face_method,                                 &
         porous_implicitness,                                    &   
         viscous_number,                                         &
         viscous_implicitness,                                   & 
         mass_limiter,                                           &
         mass_limiter_cutoff,                                    &
         strain_limit,                           &
         surften_number,                                         &

         t, dt_constant,                                         &
         dt_init, dt_grow, dt_min, dt_max,                       &
         cycle_number, cycle_max,                                &

         interface_smoothing_length,                             &
         interface_topology_model,                               &

         projection_linear_solution,                             &
         displacement_linear_solution,                           &
         displacement_nonlinear_solution,                        &

         ! JSED
         viscous_linear_solution,                                &

         flux_vol_iter_max,                                      &
         volume_track_iter_tol, volume_track_iter_max,           &
         volume_track_brents_method, volume_track_interfaces,    &
         interface_geometry, interface_area,                     &
         volume_track_subcycles,                                 &

         discrete_ops_type,                                      &

         ff_discrete_ops_type,                                   &

         stress_reduced_integration,                             &
         contact_distance,                                       &
         contact_norm_trac,                                         &
         contact_penalty,                                        &

         limiter_type,                                           &
         advection_order_vol,                                    &
         advection_order_energy,                                 &
         advection_order_momentum,                               &
         advection_order_species,                                &

         mechanical_energy_bound,                                &
         MinFaceFraction

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize relevant variables.
    fatal = .false.

    ! Inform user that the NUMERICS namelist is being read.
    call TLS_info ('')
    call TLS_info (' Reading NUMERICS Namelist ...')

    ! Default variables in this namelist.
    call NUMERICS_DEFAULT

    ! Read the namelist only on the I/O Root PE.
    IO_PE_ONLY: if (p_info%IOP) then

       ! Find namelist.
       rewind lun
       call seek_to_namelist (lun, 'NUMERICS', found)
       fatal = .not.found
       message = 'NUMERICS namelist not found'

       ! Read namelist.
       if (.not. fatal) then
          read (lun, NML = NUMERICS, IOSTAT = ioerror)
          fatal = (ioerror /= 0)
          message = 'error reading NUMERICS namelist'
       end if

    end if IO_PE_ONLY

    ! Check for fatal errors before broadcasting data to other PE's.
    call TLS_fatal_if_any (fatal, 'NUMERICS_INPUT: ' // trim(message))

    ! Broadcast all variables in the NUMERICS namelist.
    call NUMERICS_INPUT_PARALLEL

    ! Check for fatal input errors.
    call NUMERICS_CHECK (fatal)
    call TLS_fatal_if_any (fatal, 'NUMERICS_INPUT: NUMERICS namelist input error!')

  END SUBROUTINE NUMERICS_INPUT

  SUBROUTINE NUMERICS_CHECK (fatal)
    !=======================================================================
    ! Purpose(s):
    !
    !   Check NUMERICS namelist.
    !
    !=======================================================================
    use input_utilities,          only: NULL_R
    use cutoffs_module,           only: cutvof
    use discrete_ops_data,        only: use_ortho_face_gradient,              &
                                        discrete_ops_type
    use ff_discrete_ops_data,     only: use_ff_ortho_face_gradient,           &
                                        ff_discrete_ops_type
    use body_data_module,         only: body_force_implicitness
    use fluid_data_module,        only: fluid_flow, MinFaceFraction,          &
                                        momentum_solidify_implicitness,       &
                                        mass_limiter, mass_limiter_cutoff
    use flux_volume_module,       only: flux_vol_iter_max
    use interface_module,         only: interface_topology_model,             &
                                        interface_model_forms,                &
                                        interface_model_default
    use linear_solution,          only: UBIK_PRESSURE_DEFAULT,                &
                                        ubik_viscous_default,                 &
                                        UBIK_NK_DEFAULT,                      &
                                        linear_solutions, Ubik_user,          &
                                        DEFAULT_UBIK_CONTROLS,                &
                                        UBIK_DISPLACEMENT_DEFAULT,            &
                                        PRECOND_NONE,                         &
                                        PRECOND_DIAGONAL
    use nonlinear_solution,       only: NK_DEFAULT, NKuser,                   &
                                        nonlinear_solutions,                  &
                                        DEFAULT_NK_CONTROLS
    use parameter_module,         only: string_len, max_topology_models
    use legacy_mesh_api,          only: ndim
    use legacy_matl_api,          only: nmat
    use porous_drag_data,         only: porous_implicitness
    use projection_data_module,   only: projection_linear_solution,           &
                                        UBIK_PRESSURE
    use time_step_module,         only: constant_dt, courant_number, dt,      &
                                        dt_constant, dt_constraint, dt_grow,  &
                                        dt_init, t, viscous_number,           &
                                        cycle_max, cycle_number, dt_min,      &
                                        dt_max, surften_number
    use utilities_module,         only: STRING_COMPARE

    ! JSED
    use viscous_data_module,    only: viscous_linear_solution,  &
                                      ubik_viscous, viscous_implicitness

    use mollify,                  only: interface_smoothing_length
    use vof_data_module,          only: volume_track_interfaces,              &
                                        interface_geometry,                   &
                                        volume_track_subcycles,               &
                                        volume_track_iter_tol,                &
                                        volume_track_iter_max
    use solid_mechanics_input,    only: solid_mechanics,                  &
                                        UBIK_DISPLACEMENT,                &
                                        displacement_linear_solution,     &
                                        displacement_nonlinear_solution,  &
                                        NK_DISPLACEMENT, strain_limit
    use advection_data,         only: limiter_type,               &
                                      advection_order_energy,     &
                                      advection_order_momentum,   &
                                      advection_order_species

    ! Argument List
    logical, intent(INOUT) :: fatal

    ! Local Variables
    logical :: strings_match, this_string_matches
    character(80) :: initialized_string
    character(string_len) :: string
    integer :: i, j, k, l
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize fatal flag.
    fatal = .false.

    ! Start time should be >= 0.0
    if (t < 0) then
       write (message, 5) t
5      format ('time = ',1pe13.5,' is invalid; must be >= 0.0!')
       call TLS_error (message)
       fatal = .true.
    end if

    ! Maximum allowed cycle number should be > 0.
    if (cycle_max <= 0) then
       write (message, 6) cycle_max
6      format ('maximum cycle number = ',i6,' is invalid; must be > 0!')
       call TLS_error (message)
       fatal = .true.
    end if

    ! Initial cycle number should be >= 0.
    if (cycle_number < 0) then
       write (message, 7) cycle_number
7      format ('cycle number = ',i6,' is invalid; must be > 0!')
       call TLS_error (message)
       fatal = .true.
    end if

    ! Maximum allowed flux volume iterations should be > 0.
    if (flux_vol_iter_max <= 0) then
       write (message, 8) flux_vol_iter_max
8      format ('flux_vol_iter_max = ',i6,' is invalid; must be > 0!')
       call TLS_error (message)
       fatal = .true.
    end if

    ! cutoff parameters
    ! Check cutvof
    if (cutvof < 0 .or. cutvof > 1) then
       write (message, 15) cutvof
15     format ('Invalid cutvof = ',1pe13.5)
       call TLS_error (message)
       fatal = .true.
    end if

    ! Mass limiter option -- try to set some sane defaults
    ! Use an exponential function over 3-decades of volume-fraction
    if (mass_limiter) then
       if (mass_limiter_cutoff <= cutvof .or. mass_limiter_cutoff > 1) then
          mass_limiter_cutoff = 1.0e3*cutvof
          if (cutvof > 1.0d-3) then
             mass_limiter_cutoff = 1
          endif
          write (message, 16) mass_limiter_cutoff
16        format ('Invalid mass_limiter_cutoff -- reset to = ',1pe13.5)
          call TLS_warn (message)
       endif
    endif

    ! timestep parameters
    ! Check dt_grow
    if (dt_grow < 1) then
       write (message, 20) dt_grow
20     format ('Invalid dt_grow = ',1pe13.5)
       call TLS_error (message)
       fatal = .true.
    end if

    ! Check dt_constant
    if (dt_constant <= 0 .and. dt_constant /= NULL_R) then
       write (message, 25) dt_constant
25     format ('dt_constant = ',1pe13.5,' is invalid; must be > 0.0!')
       call TLS_error (message)
       fatal = .true.
    end if

    ! Check dt_init
    if (dt_init <= 0) then
       write (message, 26) dt_init
26     format ('dt_init = ',1pe13.5,' is invalid; must be > 0.0!')
       call TLS_error (message)
       fatal = .true.
    end if

    ! Check dt_max
    if (dt_max <= 0) then
       write (message, 27) dt_max
27     format ('dt_max = ',1pe13.5,' is invalid; must be > 0.0!')
       call TLS_error (message)
       fatal = .true.
    end if

    ! Check dt_min
    if (dt_min <= 0) then
       write (message, 28) dt_min
28     format ('dt_min = ',1pe13.5,' is invalid; must be > 0.0!')
       call TLS_error (message)
       fatal = .true.
    end if

    ! fluid flow parameters
    ! Check Courant Number; should be 0.0 < C <= 1.0
    if (courant_number <= 0 .or. courant_number > 1) then
       write (message, 30) courant_number
30     format ('Invalid courant_number = ',1pe13.5)
       call TLS_error (message)
       fatal = .true.
    end if

    ! Check MinFaceFraction
    if (MinFaceFraction <= 0) then
       write (message, 29) MinFaceFraction
29     format ('MinFaceFraction = ',1pe13.5,' is invalid; must be > 0.0!')
       call TLS_error (message)
       fatal = .true.
    end if

    ! mollification/surface tension parameters
    ! Check interface_smoothing_length; should be > 0.0.
    if (interface_smoothing_length <= 0 ) then
       write (message, 32) interface_smoothing_length
32     format ('interface_smoothing_length = ',1pe13.5,' is invalid; must be > 0.0!')
       call TLS_error (message)
       fatal = .true.
    end if

    ! surften_number used as a coefficient in front of dt_surften
    ! 0 < surften_number <=1 
    if (surften_number <= 0 .or. surften_number > 1) then
       write (message, 33) surften_number
33     format ('Invalid surften_number = ',1pe13.5)
       call TLS_error (message)
       fatal = .true.
    end if

    ! conduction parameters

    ! momentum solidification time weighting
    ! check momentum_solidify_implicitness
    if (momentum_solidify_implicitness<0 .or. momentum_solidify_implicitness>1) then
       write (message, 41) momentum_solidify_implicitness
41     format ('Invalid momentum_solidify_implicitness = ',1pe13.5)
       call TLS_error (message)
       fatal = .true.
    end if

    ! body force time weighting
    ! check body_force_implicitness
    if (body_force_implicitness<0 .or. body_force_implicitness>1) then
       write (message, 42) body_force_implicitness
42     format ('Invalid body_force_implicitness = ',1pe13.5)
       call TLS_error (message)
       fatal = .true.
    end if

    ! porous drag parameters
    ! check porous_implicitness
    if (porous_implicitness < 0 .or. porous_implicitness > 1) then
       write (message, 43) porous_implicitness
43     format ('Invalid porous_implicitness = ',1pe13.5)
       call TLS_error (message)
       fatal = .true.
    end if

    ! viscous parameters
    ! Check Viscous Implicitness
    if (viscous_implicitness < 0 .or. viscous_implicitness > 1) then
       write (message, 45) viscous_implicitness
45     format ('Invalid viscous_implicitness = ',1pe13.5)
       call TLS_error (message)
       fatal = .true.
    end if

    ! Set Viscous Number
    if (viscous_number == NULL_R) then
      viscous_number = 0
      if (viscous_implicitness == 0) then
        if (ndim == 2) viscous_number = 0.25_r8
        if (ndim == 3) viscous_number = 1.0_r8/6.0_r8
       else if (viscous_implicitness > 0 .and. viscous_implicitness < 1) then
        if (ndim == 2) viscous_number = 0.25_r8/(1 - viscous_implicitness)
        if (ndim == 3) viscous_number = (1.0_r8/6.0_r8)/(1 - viscous_implicitness)
      end if
    else 
      if (viscous_implicitness == 0) then
        if (ndim == 2) viscous_number = MIN(0.25_r8, viscous_number)
        if (ndim == 3) viscous_number = MIN((1.0_r8/6.0_r8), viscous_number)
      else if (viscous_implicitness > 0 .and. viscous_implicitness < 1) then
        if (ndim == 2) viscous_number = MIN(0.25_r8/(1 - viscous_implicitness), viscous_number)
        if (ndim == 3) viscous_number = MIN((1.0_r8/6.0_r8)/(1 - viscous_implicitness), viscous_number)
      end if
    endif

    ! Check Viscouus Number
    if (viscous_number < 0) then
       write (message, 50) viscous_number
50     format ('Invalid viscous_number = ',1pe13.5)
       call TLS_error (message)
       fatal = .true.
    end if

    ! Check Strain Limit
    if (strain_limit < 0) then
       write (message, 55) strain_limit
55     format ('Invalid strain_limit = ',1pe13.5)
       call TLS_error (message)
       fatal = .true.
    end if

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

    ! Set the limiter type
    string = ADJUSTL(limiter_type)
    strings_match = .false.
    call STRING_COMPARE (TRIM(string), 'default', this_string_matches)
    if (this_string_matches) then
       strings_match = .true.
       limiter_type = 'Barth'
    end if
    call STRING_COMPARE (TRIM(string), 'Barth', this_string_matches)
    if (this_string_matches) then
       strings_match = .true.
    endif
    call STRING_COMPARE (TRIM(string), 'Venkat', this_string_matches)
    if (this_string_matches) then
       strings_match = .true.
    endif
    if (.not. strings_match) then
       write (message,19) TRIM(string)
19     format ('limiter_type "',a,'" is invalid!')
       call TLS_error (message)
       fatal = .true.
    end if

    ! Set discrete ops
    string = ADJUSTL(discrete_ops_type)
    strings_match = .false.
    call STRING_COMPARE (TRIM(string), 'default', this_string_matches)
    if (this_string_matches) then
       strings_match = .true.
       discrete_ops_type  = 'default'
    end if
    call STRING_COMPARE (TRIM(string), 'ortho', this_string_matches)
    if (this_string_matches) then
       use_ortho_face_gradient = .true.
       strings_match = .true.
    endif
    call STRING_COMPARE (TRIM(string), 'nonortho', this_string_matches)
    if (this_string_matches) then
       use_ortho_face_gradient = .false.
       strings_match = .true.
    endif
    if (.not. strings_match) then
       write (message,49) TRIM(string)
49     format ('Discrete_ops_type "',a,'" not valid!')
       call TLS_error (message)
       fatal = .true.
    end if

    ! Set fluid flow discrete ops
    string = ADJUSTL(ff_discrete_ops_type)
    strings_match = .false.
    call STRING_COMPARE (TRIM(string), 'default', this_string_matches)
    if (this_string_matches) then
       strings_match = .true.
       ff_discrete_ops_type  = 'default'
    end if
    call STRING_COMPARE (TRIM(string), 'ortho', this_string_matches)
    if (this_string_matches) then
       use_ff_ortho_face_gradient = .true.
       use_ortho_face_gradient = .true.
       strings_match = .true.
    endif
    call STRING_COMPARE (TRIM(string), 'nonortho', this_string_matches)
    if (this_string_matches) then
       use_ff_ortho_face_gradient = .false.
       use_ortho_face_gradient = .false.
       strings_match = .true.
    endif
    if (.not. strings_match) then
       write (message,501) TRIM(string)
501     format ('ff_discrete_ops_type "',a,'" not valid!')
       call TLS_error (message)
       fatal = .true.
    end if

    ! find out type of interface topology model (and over write with a known string)
    do l = 1,max_topology_models ! check if relation matches predefined strings
       call STRING_COMPARE (interface_topology_model, interface_model_forms(l), &
            strings_match)
       if (strings_match) exit
    end do

    ! Interface parameters...
    if (volume_track_interfaces) then
       ! Make sure we have > 1 materials before we interface track.
       ! Set volume_track_interfaces = .false. if nmat = 1.
       if (nmat <= 1) then
          write (message, 65) nmat
65        format ('Cannot track interfaces when number of materials = ', &
                  i1,'volume_track_interfaces will be set to .false.!')
          call TLS_warn (message)
          volume_track_interfaces = .false.
       end if
       ! Make sure the number of subcyles is within the allowed ranged.
       if (volume_track_subcycles <= 0 .or. volume_track_subcycles > 20) then
          write (message, 67) volume_track_subcycles
67        format ('volume_track_subcycles = ',i3,' is invalid; must be > 0 and <= 20!')
          call TLS_error (message)
          fatal = .true.
       end if
       ! Make sure the volume fraction iteration tolerance is within the allowed range.
       if (volume_track_iter_tol <= 0.0 .or. volume_track_iter_tol > 0.001) then
          write (message, 71) volume_track_iter_tol
71        format ('volume_track_iter_tol = ',1pe13.5,' is invalid; must be > 0.0 and <= 0.001!')
          call TLS_error (message)
          fatal = .true.
       end if
       ! Make sure the volume track max allowed iterations is within the allowed range.
       if (volume_track_iter_max <= 0 .or. volume_track_iter_max > 100) then
          write (message, 72) volume_track_iter_max
72        format ('volume_track_iter_max = ',i3,' is invalid; must be > 0 and <= 100!')
          call TLS_error (message)
          fatal = .true.
       end if
    end if

    call INTERFACE_MODEL_DEFAULT ()

    initialized_string = 'none'
    ! default interface_topology_model when model type is not set in input
    call STRING_COMPARE (interface_topology_model, initialized_string, strings_match)
    if (strings_match) interface_topology_model = 'least squares model' ! set default

    ! find out type of interface topology model (and over write with a known string)
    do l = 1,max_topology_models ! check if relation matches predefined strings
       call STRING_COMPARE (interface_topology_model, interface_model_forms(l), &
            strings_match)
       if (strings_match) exit
    end do

    TOPOLOGY_MODEL_SELECT: select case (l) ! overwrite for consistent string

       case(1)
          interface_topology_model = 'least squares model'
       case(2)
          interface_topology_model = 'convolution model'

       case DEFAULT ! model doesn't match predefined values: error

          write (message,23) trim(interface_topology_model)
23        format('interface_topology_model has an unknown interface topology model type: ',a)
          call TLS_error (message)
          fatal = .true.

    end select TOPOLOGY_MODEL_SELECT

    ! Check the reconstruction geometry flag.
    if (interface_geometry /= 'piecewise linear') then
       string = ADJUSTL(interface_geometry)
       strings_match = .false.
       call STRING_COMPARE (TRIM(string), 'linear', this_string_matches)
       if (this_string_matches) string = 'piecewise linear'
       strings_match = strings_match .or. this_string_matches
       call STRING_COMPARE (TRIM(string), 'plic', this_string_matches)
       if (this_string_matches) string = 'piecewise linear'
       strings_match = strings_match .or. this_string_matches
       call STRING_COMPARE (TRIM(string), 'piecewise constant', this_string_matches)
       if (this_string_matches) string = 'piecewise constant'
       strings_match = strings_match .or. this_string_matches
       call STRING_COMPARE (TRIM(string), 'slic', this_string_matches)
       if (this_string_matches) string = 'piecewise constant'
       strings_match = strings_match .or. this_string_matches
       call STRING_COMPARE (TRIM(string), 'constant', this_string_matches)
       if (this_string_matches) string = 'piecewise constant'
       strings_match = strings_match .or. this_string_matches
       if (.not. strings_match) then
          write (message,66) TRIM(string)
66        format ('Interface geometry "',a,'" not valid!')
          call TLS_error (message)
          fatal = .true.
       end if
    end if

    ! Projection linear solution. If we can't find a
    ! solution name match, keep default and warn user.
    LS_PROJ_CHECK: do i = 1, linear_solutions
       j = i + DEFAULT_UBIK_CONTROLS
       string = Ubik_user(j)%name
       this_string_matches = .false.
       call STRING_COMPARE (TRIM(projection_linear_solution), string, this_string_matches)
       if (this_string_matches) then
          UBIK_PRESSURE = j
          if (fluid_flow) then
             write (message, 77) TRIM(projection_linear_solution)
77           format ('Using linear solver "',a,'" for projection linear solution.')
             call TLS_info (message)
          end if
          exit LS_PROJ_CHECK
       end if
    end do LS_PROJ_CHECK

    if (UBIK_PRESSURE == UBIK_PRESSURE_DEFAULT) then
       if (fluid_flow) then
          if (linear_solutions == 0) then
             write (message, 78) 
78           format ('Using default linear solver parameters for projection linear solution.')
             call TLS_info (message)
          else if (linear_solutions > 0) then
             write (message, 79) TRIM(projection_linear_solution)
79           format ('Linear solver "',a,'" for projection linear solution not found!', &
                    ' Reverting to default linear solver parameters.')
             call TLS_warn (message)
          end if
       end if
    end if

    ! JSED

    ! viscous_linear_solution. If we can't find a
    ! solution name match, keep default and warn user.
    VISCOUS_CHECK: do i = 1, linear_solutions
       j = i + DEFAULT_UBIK_CONTROLS
       string = Ubik_user(j)%name
       this_string_matches = .false.
       call STRING_COMPARE (TRIM(viscous_linear_solution), string, this_string_matches)
       if (this_string_matches) then
          ubik_viscous = j
          if (fluid_flow) then
             write (message, 97) TRIM(viscous_linear_solution)
97           format ('Using linear solver "',a,'" for viscous linear solution.')
             call TLS_info (message)
          end if
          exit VISCOUS_CHECK
       end if
    end do VISCOUS_CHECK
    
    if (ubik_viscous == ubik_viscous_default) then
       if (fluid_flow) then
          if (linear_solutions == 0) then
             write (message, 98) 
98           format ('Using default linear solver parameters for viscous linear solution.')
             call TLS_info (message)
          else if (linear_solutions > 0) then
             write (message, 99) TRIM(viscous_linear_solution)
99           format ('Linear solver "',a,'" for viscous linear solution not found!', &
                     'Reverting to default linear solver parameters.')
             call TLS_warn (message)
          end if
       end if
    else
       ! This test is somewhat misplaced, and a coherent test of the 
       ! linear solver should be performed when the linear solver 
       ! namelist is being parsed. (MAC)
       if (fluid_flow) then
          if (Ubik_user(ubik_viscous)%precond /= precond_DIAGONAL .and. &
              Ubik_user(ubik_viscous)%precond /= precond_none) then 
             write (message,100) TRIM(viscous_linear_solution)
100          format('Linear solver "',a, &
                '" for viscous solution uses an invalid preconditioner')
             call TLS_error(message)
             fatal = .true.
          endif

       end if
    end if

    ! Displacement linear solution.
    DISPLACEMENT_LS_CHECK: do i = 1, linear_solutions
       j = i + DEFAULT_UBIK_CONTROLS
       string = Ubik_user(j)%name
       this_string_matches = .false.
       call STRING_COMPARE (TRIM(displacement_linear_solution), string, this_string_matches)
       if (this_string_matches) then
          UBIK_DISPLACEMENT = j
          if (solid_mechanics) then
             write (message, 91) TRIM(displacement_linear_solution)
91           format ('"',a,'" is deprecated and has no effect on the solution!', &
                      'Please remove it from the namelist')
             call TLS_warn (message)
          end if
          exit DISPLACEMENT_LS_CHECK
       end if
    end do DISPLACEMENT_LS_CHECK

    if (UBIK_DISPLACEMENT == UBIK_DISPLACEMENT_DEFAULT .and. solid_mechanics) then
       if (displacement_linear_solution == 'default') then
          write (message, 92) 
92        format ('Using default linear solver parameters for displacement linear solution.')
          call TLS_info (message)
       else
          write (message, 93) TRIM(displacement_linear_solution)
93        format ('Linear solver "',a,'" for displacement linear solution not found!', &
                  ' Reverting to default linear solver parameters.')
          call TLS_warn (message)
       end if
    end if

    DISPLACEMENT_NLS_CHECK: do i = 1, nonlinear_solutions
       j = i + DEFAULT_NK_CONTROLS
       string = NKuser(j)%name
       this_string_matches = .false.
       call STRING_COMPARE (TRIM(displacement_nonlinear_solution), string, this_string_matches)
       if (this_string_matches) then
          NK_DISPLACEMENT = j
          if (solid_mechanics) then
             write (message, 81) TRIM(displacement_nonlinear_solution)
81           format ('Using nonlinear solver "',a,'" for displacement nonlinear solution.')
             call TLS_info (message)
          end if
          exit DISPLACEMENT_NLS_CHECK
       end if
    end do DISPLACEMENT_NLS_CHECK

    if (NK_DISPLACEMENT == NK_DEFAULT .and. solid_mechanics) then
       if (displacement_nonlinear_solution == 'default') then
          write (message, 82) 
82        format ('Using default nonlinear solver parameters for displacement nonlinear solution.')
          call TLS_info (message)
       else
          write (message, 83) TRIM(displacement_nonlinear_solution)
83        format ('Nonlinear solver "',a,'" for displacement nonlinear solution not found!', &
                  ' Reverting to default nonlinear solver parameters.')
          call TLS_warn (message)
       end if
    end if

    ! Find the correct linear solver for each input nonlinear solver.
    NLSLS_MATCH: do i = 1, nonlinear_solutions

       j = i + DEFAULT_NK_CONTROLS
       string = NKuser(j)%linear_solver_name
       this_string_matches = .false.
       LS_NAME: do l = 1,linear_solutions

          k = l + DEFAULT_UBIK_CONTROLS
          call STRING_COMPARE (TRIM(Ubik_user(k)%name), string, this_string_matches)
          if (this_string_matches) then
             NKuser(j)%linear_solver_index = k
             write (message, 58) TRIM(Ubik_user(k)%name), TRIM(NKuser(j)%name)
58           format ('Using linear solver "',a,'" for nonlinear solver "',a,'".')
             call TLS_info (message)
             exit LS_NAME
          end if

          if (NKuser(j)%linear_solver_index == UBIK_NK_DEFAULT) then
             write (message, 61) TRIM(Ubik_user(k)%name), TRIM(NKuser(j)%name)
61           format ('Using default linear solver "',a,'" for nonlinear solver "',a,'".')
             call TLS_info (message)
          end if
          
       end do LS_NAME
    end do NLSLS_MATCH

    ! various checks on flags for advection schemes...

    if (volume_track_interfaces) then
       if (advection_order_energy > 1) then
          call TLS_info ('Using High Order Advection for Energy With Interfaces')
          fatal = .false.
       end if
       if (advection_order_momentum > 1) then
          call TLS_info ('Using High Order Advection for Momentum With Interfaces')
          fatal = .false.
       end if
       if (advection_order_species > 1) then
          call TLS_info ('Using High Order Advection for Species With Interfaces')
          fatal = .false.
       end if
    end if

    if (advection_order_energy < 1 .or. advection_order_energy > 2) then
       call TLS_error ('Bad order chosen for Energy advection')
       fatal = .true.
    end if
    if (advection_order_momentum < 1 .or. advection_order_momentum > 2) then
       call TLS_error ('Bad order chosen for Momentum advection')
       fatal = .true.
    end if
    if (advection_order_species < 1 .or. advection_order_species > 2) then
       call TLS_error ('Bad order chosen for Species advection')
       fatal = .true.
    end if

  END SUBROUTINE NUMERICS_CHECK

  SUBROUTINE NUMERICS_DEFAULT
    !=======================================================================
    ! Purpose(s):
    !
    !   Default NUMERICS namelist and related variables.
    !
    !=======================================================================
    use input_utilities,        only: NULL_R
    use cutoffs_module,         only: cutvof
    use discrete_ops_data,      only: discrete_ops_type
    use ff_discrete_ops_data,   only: ff_discrete_ops_type
    use body_data_module,       only: body_force_implicitness,             &
                                      body_force_face_method
    use fluid_data_module,      only: MinFaceFraction,                     &
                                      momentum_solidify_implicitness,      &
                                      mass_limiter, mass_limiter_cutoff
    use flux_volume_module,     only: flux_vol_iter_max
    use interface_module,       only: interface_topology_model
    use linear_solution,        only: UBIK_PRESSURE_DEFAULT, UBIK_DISPLACEMENT_DEFAULT, ubik_viscous_default
    use nonlinear_solution,     only: NK_DEFAULT
    use legacy_matl_api,        only: nmat
    use porous_drag_data,       only: porous_implicitness
    use projection_data_module, only: mac_projection_iterations,           &
                                      mac_projection_precond_iter,         &
                                      projection_precond_iterations,       &
                                      projection_iterations,               &
                                      projection_linear_solution,          &
                                      UBIK_PRESSURE
    use mollify,                only: interface_smoothing_length
    use time_step_module,       only: constant_dt, courant_number, cycle_max, &
                                      cycle_number, cycle_number_restart,     &
                                      dt_constant, dt_courant,                &
                                      dt_grow, dt_init, dt_max, dt_min,       &
                                      dt_viscous, t, viscous_number, &
                                      surften_number, dt_surften

    ! JSED
    use viscous_data_module,    only: viscous_linear_solution,    &
                                      ubik_viscous,               &
                                      viscous_implicitness,       &
                                      viscous_iterations

    use vof_data_module,        only: count_cases, Eps,           &
                                      volume_track_brents_method, &
                                      volume_track_subcycles,     &
                                      volume_track_iter_max,      &
                                      volume_track_iter_tol,      &
                                      volume_track_interfaces,    &
                                      interface_geometry,         &
                                      interface_area
    use solid_mechanics_input,  only: UBIK_DISPLACEMENT, displacement_linear_solution, &
                                      NK_DISPLACEMENT, displacement_nonlinear_solution, &
                                      contact_distance, contact_norm_trac, contact_penalty, &
                                      strain_limit, stress_reduced_integration
    use advection_data,         only: limiter_type,               &
                                      advection_order_vol,        &  
                                      advection_order_energy,     &
                                      advection_order_momentum,   &
                                      advection_order_species
    use body_data_module,       only: mechanical_energy_bound

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Times.
    t =  0                          ! current time

    ! Time step information.
    cycle_number         = 0        ! current cycle number
    cycle_number_restart = 0        ! restart cycle number
    cycle_max            = 1000000  ! maximum cycle number

    dt_init       = 1.0d-6    ! initial dt
    dt_min        = 1.0d-6    ! min dt
    dt_max        = 10.0            ! max dt
    dt_grow       = 1.05            ! dt growth factor

    dt_constant   = NULL_R          ! constant dt
    constant_dt   = .false.         ! constant dt flag

    dt_courant    = 1.0d10    ! courant time step
    dt_viscous    = 1.0d10    ! viscous time step
    dt_surften    = 1.0d10    ! surface tension time step

    courant_number          = 0.50      ! max courant number
    viscous_number          = NULL_R    ! max viscous number

    ! Defaults for time-weighting for flow switched to 1/2    9/19/07 (MAC)
    ! with the exception of momentum_solidy_implicitness
    ! which is forced to default to 1.0.                     12/06/07 (MAC)
    viscous_implicitness           = 0.5_r8  ! viscous implicitness (Crank-Nicolson)
    porous_implicitness            = 0.5_r8  ! porous implicitness
    momentum_solidify_implicitness = 1.0_r8       ! momentum treatment
    body_force_implicitness        = 0.5_r8  ! body force centering
    !body_force_face_method         = .false.   ! cell-centered body force
   
    ! mmfran 07/22/11
    ! The new default value for body_force_face_method is .true.
    ! back to original treatment of the body force term at cell face
    ! necessary for material ordering independent results
    ! for special case with solidification using body_force_face_method=.false. 
    ! might be better for stability reason see above comment
    body_force_face_method         = .true. 

    strain_limit            = 1.0e-10   ! minimum plastic strain increment
    surften_number          = 1       ! surface tension number

    ! Displacement linear solution parameters.
    displacement_linear_solution = 'default'       ! linear solver to use
    UBIK_DISPLACEMENT = UBIK_DISPLACEMENT_DEFAULT  ! linear solver control parameters

    ! Projection linear solution parameters.
    projection_linear_solution = 'default'   ! linear solver to use
    UBIK_PRESSURE = UBIK_PRESSURE_DEFAULT    ! linear solver control parameters

    ! JSED
    ! viscous linear solution parameters
    viscous_linear_solution = 'default'    ! linear solver to use
    ubik_viscous = ubik_viscous_default    ! linear solver control parameters

    ! NK solution parameters.
    displacement_nonlinear_solution = 'default'    ! displacement NK nonlinear solver to use
    NK_DISPLACEMENT = NK_DEFAULT

    ! Cutoff parameters.
    cutvof  = 1.0d-8 ! volume fraction cutoff

    ! Mass limiter parameters
    mass_limiter                   = .true.       ! mass limiter
    mass_limiter_cutoff            = 1.0d-5 ! cutoff for exponential mass

    ! Volume Track Parameters
    if(nmat==1) then
        volume_track_interfaces    = .false.            ! Use VOF to track interfaces?
      else
        volume_track_interfaces    = .true.
    endif
    volume_track_iter_tol      = 1.0d-8       ! Interface location tolerance.
    volume_track_subcycles     = 2                  ! Number of time subcycles.
    volume_track_iter_max      = 20                 ! Interface location max iterations.
    volume_track_brents_method = .true.             ! Use Brent's method for plane location iteration.
    flux_vol_iter_max          = 10                 ! Flux volume vertex max iterations.
    count_cases                = .false.            ! Truncation case counting flag.
    interface_geometry         = 'piecewise linear' ! Reconstruction geometry
    interface_area             =  .false.           ! Find Interface Areas
    Eps(1) = +1.0_r8                                   ! Permutation parameters.
    Eps(2) = -1.0_r8
    Eps(3) = +1.0_r8
    Eps(4) = -1.0_r8

    ! Limiter parameters.
    limiter_type = 'Barth'

    ! order, default to first.
    advection_order_vol         = 1
    advection_order_energy      = 1
    advection_order_momentum    = 1
    advection_order_species     = 1

    ! default energy bounds.
    mechanical_energy_bound     = NULL_R

    ! Surface tension parameters.
    interface_topology_model = 'least squares model'  ! topology model
    interface_smoothing_length = 0.15        ! Radius of support of convolution kernel

    ! Iteration counts.
    mac_projection_iterations     = 0
    mac_projection_precond_iter   = 0
    projection_iterations         = 0
    projection_precond_iterations = 0
    viscous_iterations            = 0

    ! Discrete_ops
    discrete_ops_type             = 'default'
    ff_discrete_ops_type          = 'default'

    ! Node based operator
    stress_reduced_integration    = .false.

    ! Contact parameters
    contact_distance              = 1.0e-7
    contact_norm_trac             = 1.0e4
    contact_penalty               = 1.0e3

    ! Projection parameters
    MinFaceFraction               = 1.0e-3

  END SUBROUTINE NUMERICS_DEFAULT

  SUBROUTINE NUMERICS_INPUT_PARALLEL
    !======================================================================
    ! Purpose(s):
    !
    !   Broadcast all elements of numerics namelist.
    !
    !======================================================================
    use cutoffs_module,         only: alittle, cutvof
    use discrete_ops_data,      only: discrete_ops_type
    use ff_discrete_ops_data,   only: ff_discrete_ops_type
    use body_data_module,       only: body_force_implicitness, &
                                      body_force_face_method
    use fluid_data_module,      only: MinFaceFraction,                      &
                                      momentum_solidify_implicitness,       &
                                      mass_limiter, mass_limiter_cutoff
    use flux_volume_module,     only: flux_vol_iter_max
    use interface_module,       only: interface_topology_model
    use parallel_info_module,   only: p_info
    use pgslib_module,          only: PGSLIB_BCAST
    use porous_drag_data,       only: porous_implicitness
    use projection_data_module, only: projection_linear_solution
    use mollify,                      only: interface_smoothing_length
    use time_step_module,       only: constant_dt, courant_number, cycle_max, &
                                      cycle_number, cycle_number_restart,     &
                                      dt_constant, dt_constraint, dt_grow,    &
                                      dt_init, dt_max, dt_min,                &
                                      t, viscous_number,      &
                                      surften_number
    use viscous_data_module,    only: viscous_implicitness, viscous_linear_solution
    use vof_data_module,        only: volume_track_brents_method, &
                                      volume_track_subcycles,     &
                                      volume_track_iter_max,      &
                                      volume_track_iter_tol,      &
                                      volume_track_interfaces,    &
                                      interface_geometry,         &
                                      interface_area
    use solid_mechanics_input,  only: displacement_linear_solution,    &
                                      displacement_nonlinear_solution, &
                                      contact_distance,                & 
                                      contact_norm_trac,               &
                                      contact_penalty,                 &
                                      strain_limit,                    &
                                      stress_reduced_integration
    use advection_data,         only: limiter_type,               &
                                      advection_order_vol,        &  
                                      advection_order_energy,     &
                                      advection_order_momentum,   &
                                      advection_order_species
    use body_data_module,       only: mechanical_energy_bound

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Broadcast the data.
    BROADCAST_VARIABLES: if (.NOT. p_info%UseGlobalServices) then
       ! Namelist variables.
       call PGSLIB_BCAST (alittle)
       call PGSLIB_BCAST (cutvof)

       call PGSLIB_BCAST (courant_number)
       call PGSLIB_BCAST (momentum_solidify_implicitness)
       call PGSLIB_BCAST (body_force_implicitness)
       call PGSLIB_BCAST (body_force_face_method)
       call PGSLIB_BCAST (porous_implicitness)
       call PGSLIB_BCAST (viscous_number)
       call PGSLIB_BCAST (viscous_implicitness)
       call PGSLIB_BCAST (mass_limiter)
       call PGSLIB_BCAST (mass_limiter_cutoff)
       call PGSLIB_BCAST (strain_limit)
       call PGSLIB_BCAST (surften_number)
       call PGSLIB_BCAST (MinFaceFraction)
       call PGSLIB_BCAST (t)
       call PGSLIB_BCAST (dt_constant)
       call PGSLIB_BCAST (dt_init)
       call PGSLIB_BCAST (dt_grow)
       call PGSLIB_BCAST (dt_min)
       call PGSLIB_BCAST (dt_max)
       call PGSLIB_BCAST (cycle_number)
       call PGSLIB_BCAST (cycle_max)

       call PGSLIB_BCAST (limiter_type)
       call PGSLIB_BCAST(advection_order_vol)
       call PGSLIB_BCAST(advection_order_energy)
       call PGSLIB_BCAST(advection_order_momentum)
       call PGSLIB_BCAST(advection_order_species)

       call PGSLIB_BCAST(mechanical_energy_bound)

       call PGSLIB_BCAST (interface_smoothing_length)
       call PGSLIB_BCAST (interface_topology_model)

       call PGSLIB_BCAST (projection_linear_solution)
       call PGSLIB_BCAST (viscous_linear_solution)
       call PGSLIB_BCAST (displacement_linear_solution)
       call PGSLIB_BCAST (displacement_nonlinear_solution)

       call PGSLIB_BCAST (discrete_ops_type)
       call PGSLIB_BCAST (ff_discrete_ops_type)

       call PGSLIB_BCAST (flux_vol_iter_max)
       call PGSLIB_BCAST (volume_track_subcycles)
       call PGSLIB_BCAST (volume_track_iter_tol)
       call PGSLIB_BCAST (volume_track_iter_max)
       call PGSLIB_BCAST (volume_track_brents_method)
       call PGSLIB_BCAST (volume_track_interfaces)
       call PGSLIB_BCAST (interface_geometry)
       call PGSLIB_BCAST (interface_area)

       call PGSLIB_BCAST (stress_reduced_integration)

       call PGSLIB_BCAST (contact_distance)
       call PGSLIB_BCAST (contact_norm_trac)
       call PGSLIB_BCAST (contact_penalty)

       ! Non-namelist variables.
       call PGSLIB_BCAST (constant_dt)
       call PGSLIB_BCAST (cycle_number_restart)
       call PGSLIB_BCAST (dt_constraint)

    end if BROADCAST_VARIABLES

  END SUBROUTINE NUMERICS_INPUT_PARALLEL

END MODULE NUMERICS_INPUT_MODULE
