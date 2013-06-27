MODULE PHYSICS_INPUT_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define various procedures for the input of global physics flags.
  !
  !   Public Interface:
  !
  !     * call PHYSICS_INPUT ()
  !
  !       Defaults, reads, checks, and broadcasts input variables
  !       in the PHYSICS namelist.
  !
  ! Contains: PHYSICS_INPUT
  !           PHYSICS_CHECK
  !           PHYSICS_DEFAULT
  !           PHYSICS_INPUT_PARALLEL
  !
  !=======================================================================
  use truchas_logging_services
  implicit none
  private
 
  ! Public Subroutines
  public :: PHYSICS_INPUT
 
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  
  logical :: electromagnetics
 
CONTAINS
 
  SUBROUTINE PHYSICS_CHECK (fatal)
    !=======================================================================
    ! Purpose:
    !
    !   Check PHYSICS namelist
    !=======================================================================
    use physics_module
    use fluid_data_module,      only: fluid_flow, applyflow
    use porous_drag_data,       only: porous_flow
    use surface_tension_module, only: surface_tension
    use viscous_data_module,    only: inviscid, stokes
 
    ! argument list
    logical, intent(INOUT) :: fatal
 
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    
    if (surface_tension .and. .not. fluid_flow) then
      call TLS_error ('Cannot enable SURFACE_TENSION without also enabling FLUID_FLOW.')
      fatal = .true.
    end if

    ! Check viscous and fluid flow flags.
    FLUID_FLOW_CHECK: if (stokes) then
       if (inviscid) then
          inviscid = .false.
          call TLS_warn ('inviscid assigned .false.')
       end if
       if (.not. fluid_flow) then
          fluid_flow = .true.
          call TLS_warn ('fluid_flow assigned .true.')
       end if
    end if FLUID_FLOW_CHECK

    if (applyflow) then
       fluid_flow = .true.
       call TLS_warn ('applyflow is .true. so fluid_flow assigned .true.')
    end if
 
    ! Porous flow checks
    POROUS_FLOW_CHECK: if (porous_flow) then
       if (.not. fluid_flow) then
          call TLS_error ('fluid_flow must be true with the porous flow model on')
          fatal = .true.
       end if
    end if POROUS_FLOW_CHECK
 
    if (count([heat_transport, species_transport, heat_species_transport]) > 1) then
      call TLS_error ('At most one of HEAT_TRANSPORT, SPECIES_TRANSPORT, and &
                      &HEAT_SPECIES_TRANSPORT can be enabled.')
    end if

    ! Check the number_of_species value.
    if (species_transport .or. heat_species_transport) then
      if (number_of_species <= 0) then
        call TLS_error ('NUMBER_OF_SPECIES must be > 0')
        fatal = .true.
      end if
    end if

  END SUBROUTINE PHYSICS_CHECK
 
  SUBROUTINE PHYSICS_DEFAULT ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Default PHYSICS namelist and related variables.
    !
    !=======================================================================
    use physics_module
    use body_data_module,       only: Body_Force
    use fluid_data_module,      only: fluid_flow, qin, qout, void_pressure, &
                                      boussinesq_approximation
    use porous_drag_data,       only: porous_flow
    use surface_tension_module, only: surface_tension
    use viscous_data_module,    only: inviscid, stokes
    use solid_mechanics_data,   only: solid_mechanics, solid_mechanics_body_force

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
 
    ! Namelist variables . . .
    ! Heat/species transport
    heat_transport = .false.
    species_transport = .false.
    heat_species_transport = .false.
    number_of_species = 0

    ! solid_mechanics
    solid_mechanics = .false.
    solid_mechanics_body_force = .false.

    ! Fluid flow flags and parameters.
    Body_Force               =  0   ! Body force
    fluid_flow               = .true.  ! Fluid flow
    inviscid                 = .true.  ! Inviscid flow
    stokes                   = .false. ! Stokes flow
    boussinesq_approximation = .true.  ! Boussinesq approximation
    porous_flow              = .false. ! Porous flow model
    surface_tension          = .false.
    
    electromagnetics = .false.
 
     ! Inflow/Outflow volumes.
    qin  = 0
    qout = 0
 
    ! Void pressure.
    void_pressure = 0

  END SUBROUTINE PHYSICS_DEFAULT
 
  SUBROUTINE PHYSICS_INPUT (lun)
    !=======================================================================
    ! Purpose(s):
    !
    !   Read the PHYSICS namelist.
    !
    !=======================================================================
    use physics_module
    use body_data_module,       only: Body_Force
    use EM_data_proxy,          only: SET_EM_SIMULATION_ON_OR_OFF
    use fluid_data_module,      only: fluid_flow, void_pressure, applyflow
    use input_utilities,        only: seek_to_namelist
    use parallel_info_module,   only: p_info
    use porous_drag_data,       only: porous_flow
    use surface_tension_module, only: surface_tension
    use viscous_data_module,    only: inviscid, stokes
    use solid_mechanics_data,   only: solid_mechanics, solid_mechanics_body_force
    use diffusion_solver_data,  only: ds_enabled, system_type, num_species
    
    integer, intent(in) :: lun

    ! Local Variables
    logical :: fatal, found, use_phys_defaults
    integer :: ioerror

    ! Define PHYSICS namelist
    Namelist /PHYSICS/ heat_transport,                 &
                       species_transport,              &
                       heat_species_transport,         &
                       number_of_species,              &
                       fluid_flow,                     &
                       inviscid,                       &
                       Body_Force,                     &
                       stokes,                         &
                       surface_tension,                &
                       porous_flow,                    &
                       solid_mechanics,                &
                       solid_mechanics_body_force,     &
                       void_pressure,                  &
                       electromagnetics,               &
                       applyflow

    ! The following option is left here for possible future use (MAC)
    ! As discussed 5/3/05 this is being removed to avoid confusion regarding
    ! options for heat transfer, etc.  By default, this option is set to
    ! ".true"
    !                  boussinesq_approximation,       &

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
 
    ! Initialize flags.
    use_phys_defaults = .false.
    fatal = .false.
 
    ! Inform user that the PHYSICS namelist is being read.
    call TLS_info ('')
    call TLS_info ('Reading PHYSICS namelist ...')
 
    ! Default variables in the PHYSICS namelist.
    call PHYSICS_DEFAULT ()
 
    ! Find and read the namelist. The actual reading is only done
    ! on the IO PE; the results are then broadcast to other PEs.
    IO_PE_ONLY: if (p_info%IOP) then
 
       ! Search for the PHYSICS namelist.
       rewind lun
       call seek_to_namelist (lun, 'PHYSICS', found)
       use_phys_defaults = .not.found
 
       ! Didn't find it; warn user.
       if (use_phys_defaults) then
          call TLS_warn ('PHYSICS namelist not found!  Using PHYSICS defaults.')
       end if
 
       ! Read namelist if not using defaults.
       if (.not. use_phys_defaults) read (lun, NML = physics, IOSTAT = ioerror)
       fatal = (ioerror /= 0)
 
       ! Error reading namelist; fatal.
       if (fatal) then
          call TLS_error ('Error reading PHYSICS namelist!')
       end if
 
    endif IO_PE_ONLY
 
    ! Broadcast all items in physics namelist to other PE's,
    ! even if we didn't read anything in.  This is a bit of
    ! extra work, but the code is cleaner. New PHYSICS namelist
    ! variables must be added to the PGSLib_BCAST call list
    ! in PHYSICS_INPUT_PARALLEL.
    call PHYSICS_INPUT_PARALLEL ()
 
    ! Check for stupid input errors.
    call PHYSICS_CHECK (fatal)
 
    ! Halt all PE's if fatal error
    call TLS_fatal_if_any (fatal, 'PHYSICS_INPUT: physics namelist input error')
 
    ! Electromagnetics is a local variable.  Now need to
    ! store the result in EM_Data_Store.
    ! Note that this handles the broadcast, also.
    call SET_EM_SIMULATION_ON_OR_OFF (electromagnetics)

    !! Configure the heat/species transport kernel.
    if (heat_species_transport) then
      ds_enabled = .true.
      system_type = 'thermal+species'
      num_species = number_of_species
    else if (heat_transport) then
      ds_enabled = .true.
      system_type = 'thermal'
    else if (species_transport) then
      ds_enabled = .true.
      system_type = 'species'
      num_species = number_of_species
    else
      ds_enabled = .false.
    end if
 
  END SUBROUTINE PHYSICS_INPUT
 
  SUBROUTINE PHYSICS_INPUT_PARALLEL ()
    !======================================================================
    ! Purpose(s):
    !
    !  Broadcast PHYSICS namelist variables to all other PEs.
    !
    !======================================================================
    use physics_module
    use body_data_module,       only: Body_Force
    use fluid_data_module,      only: fluid_flow, boussinesq_approximation, &
                                      void_pressure, applyflow
    use parallel_info_module,   only: p_info
    use pgslib_module,          only: PGSLib_BCAST
    use porous_drag_data,       only: porous_flow
    use surface_tension_module, only: surface_tension
    use viscous_data_module,    only: inviscid, stokes
    use solid_mechanics_data,   only: solid_mechanics, solid_mechanics_body_force
 
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
 
    ! Broadcast Data
    BROADCAST_VARIABLES: if (.NOT. p_info%UseGlobalServices) then

       call PGSLib_BCAST (void_pressure)
       call PGSLib_BCAST (Body_Force)
       call PGSLib_BCAST (fluid_flow)
       call PGSLib_BCAST (applyflow)
       call PGSLib_BCAST (heat_transport)
       call PGSLib_BCAST (species_transport)
       call PGSLib_BCAST (heat_species_transport)
       call PGSLib_BCAST (number_of_species)
       call PGSLib_BCAST (inviscid)
       call PGSLib_BCAST (boussinesq_approximation)
       call PGSLib_BCAST (porous_flow)
       call PGSLib_BCAST (stokes)
       call PGSLib_BCAST (surface_tension)
       call PGSLib_BCAST (solid_mechanics)
       call PGSLib_BCAST (solid_mechanics_body_force)
       call PGSLib_BCAST (electromagnetics)

    end if BROADCAST_VARIABLES
 
  END SUBROUTINE PHYSICS_INPUT_PARALLEL
 
END MODULE PHYSICS_INPUT_MODULE
