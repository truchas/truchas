!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module physics_input_module

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use truchas_logging_services
  implicit none
  private

  public :: physics_input

  logical :: electromagnetics

contains

  subroutine physics_check

    use physics_module
    use fluid_data_module,      only: applyflow
    use porous_drag_data,       only: porous_flow
    use surface_tension_module, only: surface_tension
    use viscous_data_module,    only: inviscid, stokes

    logical :: fatal

    fatal = .false.

    if (surface_tension .and. .not. legacy_flow) then
      call TLS_error ('Cannot enable SURFACE_TENSION without also enabling LEGACY_FLOW.')
      fatal = .true.
    end if

    ! Check viscous and fluid flow flags.
    FLUID_FLOW_CHECK: if (stokes) then
       if (inviscid) then
          inviscid = .false.
          call TLS_warn ('inviscid assigned .false.')
       end if
       if (.not. legacy_flow) then
          legacy_flow = .true.
          call TLS_warn ('legacy_flow assigned .true.')
       end if
    end if FLUID_FLOW_CHECK

    if (applyflow) then
       legacy_flow = .true.
       call TLS_warn ('applyflow is .true. so legacy_flow assigned .true.')
    end if

    ! Porous flow checks
    POROUS_FLOW_CHECK: if (porous_flow) then
       if (.not. legacy_flow) then
          call TLS_error ('legacy_flow must be true with the porous flow model on')
          fatal = .true.
       end if
    end if POROUS_FLOW_CHECK

    ! Check the number_of_species value.
    if (species_transport) then
      if (number_of_species <= 0) then
        call TLS_error ('NUMBER_OF_SPECIES must be > 0')
        fatal = .true.
      end if
    end if

    if (fatal) call TLS_fatal('Errors found with PHYSICS namelists variables')

  end subroutine physics_check

  subroutine physics_input(lun)

    use physics_module
    use body_data_module,       only: body_force
    use EM_data_proxy,          only: SET_EM_SIMULATION_ON_OR_OFF
    use fluid_data_module,      only: fluid_flow, void_pressure, applyflow, boussinesq_approximation
    use input_utilities,        only: seek_to_namelist, NULL_R
    use porous_drag_data,       only: porous_flow
    use surface_tension_module, only: surface_tension
    use viscous_data_module,    only: inviscid, stokes
    use solid_mechanics_input,  only: solid_mechanics, solid_mechanics_body_force
    use diffusion_solver_data,  only: ds_enabled, system_type, num_species
    use parallel_communication, only: is_IOP, broadcast
    use string_utilities, only: i_to_c

    integer, intent(in) :: lun

    logical :: found
    integer :: ios
    character(80) :: iom

    namelist /physics/ heat_transport, species_transport,  number_of_species, flow, &
                       legacy_flow, inviscid, stokes, surface_tension, porous_flow, &
                       body_force_density, prescribed_flow, applyflow, void_pressure, &
                       electromagnetics, solid_mechanics, solid_mechanics_body_force

    call TLS_info ('')
    call TLS_info ('Reading PHYSICS namelist ...')

    !! Locate the PHYSICS namelist (required)
    if (is_IOP) then
      rewind(lun)
      call seek_to_namelist(lun, 'PHYSICS', found, iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
    call broadcast(found)
    if (.not.found) call TLS_fatal('PHYSICS namelist not found')

    !! Read the PHYSICS namelist, assigning default values first
    if (is_IOP) then
      flow = .false.
      heat_transport = .false.
      species_transport = .false.
      number_of_species = 0
      electromagnetics = .false.
      solid_mechanics = .false.
      solid_mechanics_body_force = .false.
      body_force_density = 0.0_r8
      legacy_flow = .false.
      inviscid = .true.
      stokes = .false.
      boussinesq_approximation = .true.
      porous_flow = .false.
      surface_tension = .false.
      void_pressure = 0
      prescribed_flow = .false.
      applyflow = .false.
      read(lun,nml=physics,iostat=ios,iomsg=iom)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading PHYSICS namelist: ' // trim(iom))

    !! Broadcast all the namelist variables
    call broadcast(flow)
    call broadcast(void_pressure)
    call broadcast(body_force_density)
    call broadcast(legacy_flow)
    call broadcast(applyflow)
    call broadcast(heat_transport)
    call broadcast(species_transport)
    call broadcast(number_of_species)
    call broadcast(inviscid)
    call broadcast(boussinesq_approximation)
    call broadcast(porous_flow)
    call broadcast(stokes)
    call broadcast(surface_tension)
    call broadcast(solid_mechanics)
    call broadcast(solid_mechanics_body_force)
    call broadcast(electromagnetics)
    call broadcast(prescribed_flow)

    ! flow and legacy_flow are mutually exclusive
    ! flow and prescribed_flow are mutually exclusive

    ! Check for stupid input errors.
    call physics_check

    ! Configure the legacy flow solver
    fluid_flow = legacy_flow
    if (legacy_flow) body_force = body_force_density

    !NNC: temporary test code
    if (prescribed_flow .and. (fluid_flow .or. heat_transport .or. species_transport &
                                        .or. solid_mechanics .or. electromagnetics)) then
       call TLS_fatal('prescribed_flow is not compatible with any other physics')
    end if

    ! for now advection with a user-specified flow uses the flow solver
    flow = flow .or. prescribed_flow

    ! Electromagnetics is a local variable.  Now need to
    ! store the result in EM_Data_Store.
    ! Note that this handles the broadcast, also.
    call SET_EM_SIMULATION_ON_OR_OFF (electromagnetics)

    !! Configure the heat/species transport kernel.
    ds_enabled = heat_transport .or. species_transport
    if (species_transport) num_species = number_of_species
    if (heat_transport .and. species_transport) then
      system_type = 'thermal+species'
    else if (heat_transport) then
      system_type = 'thermal'
    else if (species_transport) then
      system_type = 'species'
    end if

  end subroutine physics_input

end module physics_input_module
