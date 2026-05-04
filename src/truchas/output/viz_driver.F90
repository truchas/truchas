!!
!! VIZ_DRIVER
!!
!! Public driver facade for VTKHDF visualization output. This module owns the
!! single process-wide viz_manager instance and exposes the TVO_* procedures
!! called from the existing Truchas execution flow. It is intentionally thin:
!! input has already been packed into output_control%params, and all stream,
!! provider, registry, file, and timestep work is delegated to viz_manager
!! and the objects it owns.

module viz_driver

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use viz_manager_type
  use viz_provider_class
  implicit none
  private

  public :: TVO_init
  public :: TVO_schedule_events
  public :: TVO_open_and_write_mesh
  public :: TVO_register_temporal_data
  public :: TVO_write_timestep
  public :: TVO_write_stream_timestep
  public :: TVO_close

  type(viz_manager) :: manager

contains

  subroutine TVO_init
    use output_control, only: params
    use truchas_logging_services
    integer :: stat
    character(:), allocatable :: errmsg
    call manager%init(params, stat, errmsg)
    if (stat /= 0) call TLS_fatal(errmsg)
  end subroutine

  subroutine TVO_schedule_events(event_queue)
    use sim_event_queue_type
    class(sim_event_queue), intent(inout) :: event_queue
    call manager%schedule_events(event_queue)
  end subroutine

  subroutine TVO_open_and_write_mesh(mesh)
    use unstr_mesh_type
    use truchas_logging_services
    type(unstr_mesh), intent(in) :: mesh
    integer :: stat
    character(:), allocatable :: errmsg
    call manager%open_and_write_mesh(mesh, stat, errmsg)
    if (stat /= 0) call TLS_fatal(errmsg)
  end subroutine

  subroutine TVO_register_temporal_data
    use truchas_logging_services
    integer :: stat
    character(:), allocatable :: errmsg
    type(viz_provider_box), allocatable :: catalog(:)
    call get_viz_provider_catalog(catalog)
    call manager%register_temporal_data(catalog, stat, errmsg)
    if (stat /= 0) call TLS_fatal(errmsg)
  end subroutine

  subroutine TVO_write_timestep(t)
    real(r8), intent(in) :: t
    call manager%write_timestep(t)
  end subroutine

  subroutine TVO_write_stream_timestep(stream_id, t)
    integer, intent(in) :: stream_id
    real(r8), intent(in) :: t
    call manager%write_stream_timestep(stream_id, t)
  end subroutine

  subroutine TVO_close
    call manager%close
  end subroutine

  subroutine get_viz_provider_catalog(providers)

    use diffusion_solver_viz_provider_type, only: diffusion_solver_viz_provider
    use em_heat_viz_provider_type, only: em_heat_viz_provider
    use flow_viz_provider_type, only: flow_viz_provider
    use solid_mechanics_viz_provider_type, only: solid_mechanics_viz_provider
    use thermal_state_viz_provider_type, only: thermal_state_viz_provider
    use ustruc_viz_provider_type, only: ustruc_viz_provider
    use vfrac_viz_provider_type, only: vfrac_viz_provider

    type(viz_provider_box), allocatable, intent(out) :: providers(:)

    allocate(providers(7))
    allocate(thermal_state_viz_provider :: providers(1)%p)
    allocate(diffusion_solver_viz_provider :: providers(2)%p)
    allocate(vfrac_viz_provider :: providers(3)%p)
    allocate(em_heat_viz_provider :: providers(4)%p)
    allocate(flow_viz_provider :: providers(5)%p)
    allocate(solid_mechanics_viz_provider :: providers(6)%p)
    allocate(ustruc_viz_provider :: providers(7)%p)

  end subroutine get_viz_provider_catalog

end module viz_driver
