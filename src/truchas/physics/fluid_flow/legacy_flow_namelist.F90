
!!
!! This collects the inputs for the original flow model from the PHYSICS and
!! NUMERICS namelists into a new LEGACY_FLOW namelist in preparation for the
!! eventual removal of the original flow model.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module legacy_flow_namelist

  use kinds, only: r8
  use parameter_module, only: maxmat
  use truchas_logging_services
  implicit none
  private

  public :: read_legacy_flow_namelist

  character(32), dimension(maxmat) :: material_priority

contains

  subroutine read_legacy_flow_namelist(lun)

    use parallel_communication, only: is_IOP, broadcast
    use string_utilities,       only: i_to_c
    use input_utilities,        only: seek_to_namelist
    use ff_discrete_ops_data,   only: ff_discrete_ops_type
    use body_data_module,       only: body_force_implicitness, body_force_face_method
    use fluid_data_module,      only: MinFaceFraction, momentum_solidify_implicitness, &
                                      mass_limiter, mass_limiter_cutoff, sound_speed,  &
                                      void_pressure, applyflow, boussinesq_approximation
    use flux_volume_module,     only: flux_vol_iter_max
    use interface_module,       only: interface_topology_model
    use surface_tension_module, only: surface_tension
    use porous_drag_data,       only: porous_flow, porous_implicitness, permeability_constant
    use projection_data_module, only: projection_linear_solution
    use mollify,                only: interface_smoothing_length
    use flow_time_step_module,  only: courant_number, viscous_number, surften_number
    use viscous_data_module,    only: inviscid, stokes, viscous_linear_solution, viscous_implicitness

    use vof_data_module,        only: volume_track_brents_method, &
                                      volume_track_iter_max,      &
                                      volume_track_subcycles,     &
                                      volume_track_iter_tol,      &
                                      volume_track_interfaces,    &
                                      interface_geometry,         &
                                      interface_area
    use advection_data,         only: limiter_type,               &
                                      advection_order_vol,        &
                                      advection_order_energy,     &
                                      advection_order_momentum,   &
                                      advection_order_species
    use body_data_module,       only: mechanical_energy_bound

    integer, intent(in) :: lun

    integer :: ios
    logical :: found
    character(80) :: iom

    ! Define LEGACY_FLOW namelist.
    namelist /legacy_flow/ inviscid, stokes, porous_flow, surface_tension, applyflow, &
        void_pressure, courant_number, momentum_solidify_implicitness, body_force_implicitness, &
        body_force_face_method, porous_implicitness, viscous_number, viscous_implicitness, &
        mass_limiter, mass_limiter_cutoff, surften_number, interface_smoothing_length, &
        interface_topology_model, projection_linear_solution, viscous_linear_solution, &
        flux_vol_iter_max, volume_track_iter_tol, volume_track_iter_max, &
        volume_track_brents_method, volume_track_interfaces, interface_geometry, interface_area, &
        volume_track_subcycles, ff_discrete_ops_type, limiter_type, advection_order_vol, &
        advection_order_energy, advection_order_momentum, advection_order_species, &
        mechanical_energy_bound, MinFaceFraction, material_priority, sound_speed, &
        permeability_constant

    call TLS_info('')
    call TLS_info('Reading LEGACY_FLOW namelist ...')

    !! Locate the LEGACY_FLOW namelist (required)
    if (is_IOP) then
      rewind(lun)
      call seek_to_namelist(lun, 'legacy_flow', found, iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
    call broadcast(found)
    if (.not.found) call TLS_fatal('LEGACY_FLOW namelist not found')

    ! Default variables in this namelist.
    inviscid = .true.
    stokes = .false.
    boussinesq_approximation = .true.
    porous_flow = .false.
    surface_tension = .false.
    void_pressure = 0
    applyflow = .false.
    call NUMERICS_DEFAULT

    !! Read the FLOW namelist
    if (is_IOP) read(lun,nml=legacy_flow,iostat=ios,iomsg=iom)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading LEGACY_FLOW namelist: ' // trim(iom))

    ! Broadcast all variables in the NUMERICS namelist.
    call broadcast(inviscid)
    call broadcast(stokes)
    call broadcast(boussinesq_approximation)
    call broadcast(porous_flow)
    call broadcast(surface_tension)
    call broadcast(void_pressure)
    call broadcast(applyflow)
    call NUMERICS_INPUT_PARALLEL

    ! Check for fatal input errors.
    call NUMERICS_CHECK

  end subroutine read_legacy_flow_namelist

  subroutine numerics_check

    use input_utilities,          only: NULL_R, NULL_I, NULL_C
    use cutoffs_module,           only: cutvof
    use discrete_ops_data,        only: use_ortho_face_gradient
    use ff_discrete_ops_data,     only: use_ff_ortho_face_gradient,           &
                                        ff_discrete_ops_type
    use body_data_module,         only: body_force_implicitness
    use fluid_data_module,        only: MinFaceFraction,          &
                                        momentum_solidify_implicitness,       &
                                        mass_limiter, mass_limiter_cutoff, sound_speed
    use flux_volume_module,       only: flux_vol_iter_max
    use volume_track_module,      only: matpri
    use interface_module,         only: interface_topology_model
    use linear_solution,          only: UBIK_PRESSURE_DEFAULT,                &
                                        ubik_viscous_default,                 &
                                        linear_solutions, Ubik_user,          &
                                        DEFAULT_UBIK_CONTROLS,                &
                                        PRECOND_NONE,                         &
                                        PRECOND_DIAGONAL
    use parameter_module,         only: nmat
    use legacy_mesh_api,          only: ndim
    use viscous_data_module,      only: inviscid, stokes
    use porous_drag_data,         only: porous_flow, porous_implicitness, permeability_constant
    use projection_data_module,   only: projection_linear_solution, UBIK_PRESSURE
    use flow_time_step_module,    only: courant_number, viscous_number, surften_number
    use string_utilities,         only: lower_case

    ! JSED
    use viscous_data_module,    only: viscous_linear_solution,  &
                                      ubik_viscous, viscous_implicitness

    use mollify,                  only: interface_smoothing_length
    use vof_data_module,          only: volume_track_interfaces,              &
                                        interface_geometry,                   &
                                        volume_track_subcycles,               &
                                        volume_track_iter_tol,                &
                                        volume_track_iter_max
    use advection_data,         only: limiter_type,               &
                                      advection_order_energy,     &
                                      advection_order_momentum,   &
                                      advection_order_species
    use flow_property_module,        only: get_material_id, material_name, have_void
    use f08_intrinsics, only: findloc

    integer :: j, m, n
    character(128) :: message

    if (stokes .and. inviscid) call TLS_fatal('STOKES and INVISCID are incompatible')

    if (flux_vol_iter_max <= 0) call TLS_fatal('FLUX_VOL_ITER_MAX must be >0')

    ! Mass limiter option -- try to set some sane defaults
    ! Use an exponential function over 3-decades of volume-fraction
    if (mass_limiter) then
       if (mass_limiter_cutoff <= cutvof .or. mass_limiter_cutoff > 1) then
          mass_limiter_cutoff = 1.0e3*cutvof
          if (cutvof > 1.0d-3) then
             mass_limiter_cutoff = 1
          endif
          write(message,'("MASS_LIMITER_CUTOFF invalid; reset to",es13.5)') mass_limiter_cutoff
          call TLS_warn(message)
       endif
    endif

    if (courant_number <= 0 .or. courant_number > 1) &
        call TLS_fatal('COURANT_NUMBER must be in (0,1]')

    if (MinFaceFraction <= 0) &
        call TLS_fatal('MINFACEFRACTION must be > 0.0')

    if (interface_smoothing_length <= 0 ) &
        call TLS_fatal('INTERFACE_SMOOTHING_LENGTH must be > 0.0')

    if (surften_number <= 0 .or. surften_number > 1) &
        call TLS_fatal('SURFTEN_NUMBER must be in (0,1]')

    if (momentum_solidify_implicitness<0 .or. momentum_solidify_implicitness>1) &
        call TLS_fatal('MOMENTUM_SOLIDIFY_IMPLICITNESS must be in [0,1]')

    if (body_force_implicitness<0 .or. body_force_implicitness>1) &
        call TLS_fatal('BODY_FORCE_IMPLICITNESS must be in [0,1]')

    if (porous_flow) then
      if (porous_implicitness < 0 .or. porous_implicitness > 1) &
          call TLS_fatal('POROUS_IMPLICITNESS must be in [0,1]')
      if (any(permeability_constant < 0)) &
          call TLS_fatal('PERMEABILITY_CONSTANT values must be >= 0.0')
    end if

    if (viscous_implicitness < 0 .or. viscous_implicitness > 1) &
        call TLS_fatal('VISCOUS_IMPLICITNESS must be in [0,1]')

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
    if (viscous_number < 0) call TLS_fatal('VISCOUS_NUMBER must be >= 0')

    ! Set the limiter type
    select case (lower_case(adjustl(limiter_type)))
    case ('default','barth')
      limiter_type = 'Barth'
    case ('venkat')
      limiter_type = 'Venkat'
    case default
      call TLS_fatal('unknown LIMITER_TYPE: ' // trim(limiter_type))
    end select

    ! Set fluid flow discrete ops
    ff_discrete_ops_type = lower_case(adjustl(ff_discrete_ops_type))
    select case (ff_discrete_ops_type)
    case ('default')
    case ('ortho')
      use_ff_ortho_face_gradient = .true.
      use_ortho_face_gradient = .true.
    case ('nonortho')
      use_ff_ortho_face_gradient = .false.
      use_ortho_face_gradient = .false.
    case default
      call TLS_fatal('unknown FF_DISCRETE_OPS_TYPE: ' // trim(ff_discrete_ops_type))
    end select

    if (volume_track_interfaces .and. nmat <= 1) then
      call TLS_warn('VOLUME_TRACK_INTERFACES disabled; only 1 material')
      volume_track_interfaces = .false.
    end if
    if (volume_track_interfaces) then
      if (volume_track_subcycles <= 0 .or. volume_track_subcycles > 20) &
          call TLS_fatal('VOLUME_TRACK_SUBCYCLES must be >= 1 and <= 20')
      if (volume_track_iter_tol <= 0.0 .or. volume_track_iter_tol > 0.001) &
          call TLS_fatal('VOLUME_TRACK_ITER_TOL must be in (0, 0.001]')
      if (volume_track_iter_max <= 0 .or. volume_track_iter_max > 100) &
          call TLS_fatal('VOLUME_TRACK_ITER_MAX must be >= 1 and <= 100')
    end if

    select case (lower_case(interface_topology_model))
    case ('none', 'least squares model')
      interface_topology_model = 'least squares model'
    case ('convolution model')
      interface_topology_model = 'convolution model'
    case default
      call TLS_fatal('unknown INTERFACE_TOPOLOGY_MODEL: ' // trim(interface_topology_model))
    end select

    ! Check the reconstruction geometry flag.
    select case (lower_case(interface_geometry))
    case ('piecewise linear', 'linear', 'plic')
      interface_geometry = 'piecewise linear'
    case ('piecewise constant', 'constant', 'slic')
      interface_geometry = 'piecewise constant'
    case default
      call TLS_fatal('unknown INTERFACE_GEOMETRY: ' // trim(interface_geometry))
    end select

    ! Projection linear solution. If we can't find a
    ! solution name match, keep default and warn user.
    do j = DEFAULT_UBIK_CONTROLS + 1, DEFAULT_UBIK_CONTROLS + linear_solutions
      if (lower_case(projection_linear_solution) == lower_case(Ubik_user(j)%name)) then
        UBIK_PRESSURE = j
        call TLS_info('Using linear solver "' // trim(projection_linear_solution) // &
                      '" for projection linear solution.')
        exit
      end if
    end do

    if (UBIK_PRESSURE == UBIK_PRESSURE_DEFAULT) then
      if (linear_solutions == 0) then
        call TLS_info('Using default linear solver parameters for projection linear solution.')
      else if (linear_solutions > 0) then
        call TLS_warn('Linear solver "' // trim(projection_linear_solution) // &
                      '" for projection linear solution not found! &
                      &Reverting to default linear solver parameters.')
      end if
    end if

    ! JSED

    ! viscous_linear_solution. If we can't find a
    ! solution name match, keep default and warn user.
    do j = DEFAULT_UBIK_CONTROLS + 1, DEFAULT_UBIK_CONTROLS + linear_solutions
      if (lower_case(viscous_linear_solution) == lower_case(Ubik_user(j)%name)) then
        ubik_viscous = j
        call TLS_info('Using linear solver "' // trim(viscous_linear_solution) // &
                      '" for viscous linear solution.')
        exit
      end if
    end do

    if (ubik_viscous == ubik_viscous_default) then
      if (linear_solutions == 0) then
        call TLS_info('Using default linear solver parameters for viscous linear solution.')
      else if (linear_solutions > 0) then
        call TLS_warn('Linear solver "' // trim(viscous_linear_solution) // &
                      '" for viscous linear solution not found! &
                      &Reverting to default linear solver parameters.')
      end if
    else
      ! This test is somewhat misplaced, and a coherent test of the
      ! linear solver should be performed when the linear solver
      ! namelist is being parsed. (MAC)
      if (Ubik_user(ubik_viscous)%precond /= precond_DIAGONAL .and. &
          Ubik_user(ubik_viscous)%precond /= precond_none) then
        call TLS_fatal('Linear solver "' // trim(viscous_linear_solution) // &
                       '" for viscous solution uses an invalid preconditioner')
      endif
    end if

    if (advection_order_energy < 1 .or. advection_order_energy > 2) &
        call TLS_fatal('ADVECTION_ORDER_ENERGY must be 1 or 2')
    if (advection_order_momentum < 1 .or. advection_order_momentum > 2) &
        call TLS_fatal('ADVECTION_ORDER_MOMENTUM must be 1 or 2')
    if (advection_order_species < 1 .or. advection_order_species > 2) &
        call TLS_fatal('ADVECTION_ORDER_SPECIES must be 1 or 2')

    ! Check material interface advection priority input
    if (all(material_priority(1:nmat) == NULL_C)) then
      call TLS_info('MATERIAL_PRIORITY not specified; defaulting to material input order')
      do m = 1, nmat
        material_priority(m) = material_name(m)
      end do
    else if (any(material_priority(1:nmat) == NULL_C)) then
      call TLS_fatal('MATERIAL_PRIORITY incompletely specified')
    else
    end if

    ! Convert to Truchas material IDs
    do m = 1, nmat
      n = get_material_id(material_priority(m))
      if (n == 0) then
        write(message,'(a,i0,2a)') &
            'unknown material name for MATERIAL_PRIORITY(', m, '): ', material_priority(m)
        call TLS_fatal(trim(message))
      end if
      matpri(m) = n
    end do

    ! Print material priorities
    call TLS_info ('')
    call TLS_info ('         Material Priorities')
    call TLS_info ('')
    call TLS_info ('         Material    Priority')
    call TLS_info ('         --------    --------')

    do m = 1,nmat
       write (message, 45) trim(material_name(m)), findloc(matpri, m)
45     format(12x,a8,7x,i2)
       call TLS_info (message)
    end do

    if (sound_speed < 0) then
      call TLS_fatal('SOUND_SPEED must be >= 0')
    else if (sound_speed > 0 .and. .not.have_void()) then
      call TLS_warn('No void material defined; SOUND_SPEED reset to 0')
      sound_speed = 0
    end if

  end subroutine numerics_check

  subroutine numerics_default

    use input_utilities,        only: NULL_R, NULL_I, NULL_C
    use ff_discrete_ops_data,   only: ff_discrete_ops_type
    use body_data_module,       only: body_force_implicitness,             &
                                      body_force_face_method
    use fluid_data_module,      only: MinFaceFraction,                     &
                                      momentum_solidify_implicitness,      &
                                      mass_limiter, mass_limiter_cutoff, sound_speed
    use flux_volume_module,     only: flux_vol_iter_max
    use interface_module,       only: interface_topology_model
    use linear_solution,        only: UBIK_PRESSURE_DEFAULT, ubik_viscous_default
    use parameter_module,       only: nmat
    use porous_drag_data,       only: porous_implicitness, permeability_constant
    use projection_data_module, only: mac_projection_iterations,           &
                                      mac_projection_precond_iter,         &
                                      projection_precond_iterations,       &
                                      projection_iterations,               &
                                      projection_linear_solution,          &
                                      UBIK_PRESSURE
    use mollify,                only: interface_smoothing_length
    use flow_time_step_module,  only: courant_number, viscous_number, surften_number

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
    use advection_data,         only: limiter_type,               &
                                      advection_order_vol,        &
                                      advection_order_energy,     &
                                      advection_order_momentum,   &
                                      advection_order_species
    use body_data_module,       only: mechanical_energy_bound


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

    surften_number          = 1       ! surface tension number

    ! Projection linear solution parameters.
    projection_linear_solution = 'default'   ! linear solver to use
    UBIK_PRESSURE = UBIK_PRESSURE_DEFAULT    ! linear solver control parameters

    ! JSED
    ! viscous linear solution parameters
    viscous_linear_solution = 'default'    ! linear solver to use
    ubik_viscous = ubik_viscous_default    ! linear solver control parameters

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
    ff_discrete_ops_type          = 'default'

    ! Projection parameters
    MinFaceFraction               = 1.0e-3

    material_priority = NULL_C
    sound_speed = 0.0_r8
    permeability_constant = 0.0_r8

  end subroutine numerics_default

  subroutine numerics_input_parallel

    use ff_discrete_ops_data,   only: ff_discrete_ops_type
    use body_data_module,       only: body_force_implicitness, &
                                      body_force_face_method
    use fluid_data_module,      only: MinFaceFraction,                      &
                                      momentum_solidify_implicitness,       &
                                      mass_limiter, mass_limiter_cutoff, sound_speed
    use flux_volume_module,     only: flux_vol_iter_max
    use interface_module,       only: interface_topology_model
    use parallel_communication, only: broadcast
    use porous_drag_data,       only: porous_implicitness, permeability_constant
    use projection_data_module, only: projection_linear_solution
    use mollify,                      only: interface_smoothing_length
    use flow_time_step_module,  only: courant_number, viscous_number, surften_number
    use viscous_data_module,    only: viscous_implicitness, viscous_linear_solution
    use vof_data_module,        only: volume_track_brents_method, &
                                      volume_track_subcycles,     &
                                      volume_track_iter_max,      &
                                      volume_track_iter_tol,      &
                                      volume_track_interfaces,    &
                                      interface_geometry,         &
                                      interface_area
    use advection_data,         only: limiter_type,               &
                                      advection_order_vol,        &
                                      advection_order_energy,     &
                                      advection_order_momentum,   &
                                      advection_order_species
    use body_data_module,       only: mechanical_energy_bound

    call broadcast(courant_number)
    call broadcast(momentum_solidify_implicitness)
    call broadcast(body_force_implicitness)
    call broadcast(body_force_face_method)
    call broadcast(porous_implicitness)
    call broadcast(viscous_number)
    call broadcast(viscous_implicitness)
    call broadcast(mass_limiter)
    call broadcast(mass_limiter_cutoff)
    call broadcast(surften_number)
    call broadcast(MinFaceFraction)

    call broadcast(limiter_type)
    call broadcast(advection_order_vol)
    call broadcast(advection_order_energy)
    call broadcast(advection_order_momentum)
    call broadcast(advection_order_species)

    call broadcast(mechanical_energy_bound)

    call broadcast(interface_smoothing_length)
    call broadcast(interface_topology_model)

    call broadcast(projection_linear_solution)
    call broadcast(viscous_linear_solution)

    call broadcast(ff_discrete_ops_type)

    call broadcast(flux_vol_iter_max)
    call broadcast(volume_track_subcycles)
    call broadcast(volume_track_iter_tol)
    call broadcast(volume_track_iter_max)
    call broadcast(volume_track_brents_method)
    call broadcast(volume_track_interfaces)
    call broadcast(interface_geometry)
    call broadcast(interface_area)

    call broadcast(material_priority)
    call broadcast(sound_speed)
    call broadcast(permeability_constant)

  end subroutine numerics_input_parallel

end module legacy_flow_namelist
