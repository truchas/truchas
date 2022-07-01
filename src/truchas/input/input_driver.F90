!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module input_driver
  !=======================================================================
  ! Purpose:
  !
  !   driver routine that reads the input file
  !
  ! Contains: READ_INPUT
  !
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov)
  !            Bryan Lally (lally@lanl.gov)
  !
  !=======================================================================
  implicit none
  private

  public :: read_input

contains

  subroutine read_input (infile, title)
    !=======================================================================
    ! Purpose:
    !
    !   driver for reading input file
    !=======================================================================
    use probe_namelist,            only: read_probe_namelists
    use EM_input,                  only: read_em_input
    use body_input_module,         only: interfaces_input
    use numerics_input_module,     only: numerics_input
    use outputs_input_module,      only: outputs_input
    use parallel_communication,    only: is_IOP, broadcast
    use physics_input_module,      only: physics_input
    use restart_variables,         only: restart, read_restart_namelist
    use restart_driver,            only: open_restart_file
    use EM_data_proxy,             only: em_is_on
    use mesh_manager,              only: read_truchas_mesh_namelists
    use diffusion_solver_data,     only: ds_enabled, heat_eqn
    use diffusion_solver,          only: read_ds_namelists
    use evaporation_namelist,      only: read_evaporation_namelist
    use ustruc_driver,             only: read_microstructure_namelist
    use flow_driver,               only: read_flow_namelists
    use solid_mechanics_driver,    only: read_solid_mechanics_namelists
    use physical_constants,        only: read_physical_constants
    use function_namelist,         only: read_function_namelists
    use vfunction_namelist,        only: read_vfunction_namelists
    use material_namelist,         only: read_material_namelists
    use simulation_event_queue,    only: read_simulation_control_namelist
    use toolpath_driver,           only: read_toolpath_namelists
    use toolhead_namelist,         only: read_toolhead_namelists
    use physics_module,            only: heat_transport, flow, solid_mechanics
    use advection_velocity_namelist, only: read_advection_velocity_namelist
    use body_namelist,             only: read_body_namelists
    use truchas_logging_services
    use truchas_timers
    use string_utilities, only: i_to_c

    use material_model_driver,  only: init_material_model

    character(*), intent(in)  :: infile
    character(*), intent(out) :: title

    integer :: ios, lun
    character(128) :: iom

    call start_timer ('Input')
    call TLS_info ('')
    call TLS_info ('Opening input file "' // trim(infile) // '"')

    ! open input file
    if (is_IOP) then
      open(newunit=lun,file=trim(infile),status='old',position='rewind',action='read',iostat=ios,iomsg=iom)
    else
      lun = -1
    end if
    call broadcast (ios)
    if (ios /= 0) call TLS_fatal ('error opening input file: ' // trim(iom))

    ! read first line as title information
    if (is_IOP) read(lun,'(a)',iostat=ios) title
    call broadcast (ios)
    if (ios /= 0) call TLS_fatal ('error reading initial title line: iostat=' // i_to_c(ios))
    call broadcast (title)
     
    call read_function_namelists (lun)
    call read_physical_constants (lun)
    call read_toolpath_namelists (lun)
    call read_toolhead_namelists (lun)
    call read_vfunction_namelists (lun)

    ! read current physics data
    call physics_input (lun)

    call read_material_namelists(lun)
    call init_material_model  ! later namelists will query the model

    ! read output specifications
    call outputs_input (lun)

    call read_simulation_control_namelist (lun)

    if (restart) then
      call read_restart_namelist (lun)
      call open_restart_file () ! NB: reads ncells and nnodes used later by mesh_sizes.
    end if

    ! Read the MESH and ALTMESH namelists: used to initialize MESH_MANAGER (new mesh)
    call read_truchas_mesh_namelists (lun)

    ! read volume fraction data
    call interfaces_input (lun)
    call read_body_namelists (lun)

    ! read numerical options data
    call numerics_input (lun)

    if (flow) call read_flow_namelists(lun)

    ! read namelists for solid mechanics options
    if (solid_mechanics) call read_solid_mechanics_namelists (lun)

    ! Read Electromagnetics
    if (em_is_on()) call read_em_input (lun)

    ! Read diffusion solver namelists
    if (ds_enabled) then
      if (heat_transport) then
        call read_evaporation_namelist (lun)
      end if
      call read_ds_namelists (lun)
      if (heat_eqn) call read_microstructure_namelist (lun)
    end if

    ! Read probe information.
    call read_probe_namelists(lun)

    if (is_IOP) close(lun)

    call TLS_info ('Input file "' // trim(infile) // '" closed')
    call stop_timer('Input')

  end subroutine read_input

end module input_driver
