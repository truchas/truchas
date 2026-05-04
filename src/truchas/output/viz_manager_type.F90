!!
!! VIZ_MANAGER_TYPE
!!
!! Top-level coordinator for VTKHDF visualization output. A viz_manager owns
!! the set of output streams and the active provider objects. It is the layer
!! that connects user-selected field names to package-specific provider states: active packages register their available
!! fields into a scratch registry, stream selections are resolved against that
!! registry, and per-stream provider states are created with either "all fields"
!! or the selected subset for that provider.
!!
!! The manager does not write mesh geometry or field values directly.  It
!! delegates stream/file mechanics to viz_stream_state and package-specific
!! data registration/writing to viz_provider_state implementations.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NOTES
!!
!! viz_manager lifecycle:
!!
!!   MANAGER_CLOSED
!!     -> init
!!   MANAGER_INITIALIZED
!!     -> schedule_events
!!     -> open_and_write_mesh
!!   MANAGER_MESH_WRITTEN
!!     -> schedule_events
!!     -> register_temporal_data
!!   MANAGER_TEMPORAL_REGISTERED
!!     -> schedule_events
!!     -> write_timestep (zero or more times)
!!
!! close is valid from any state and returns the manager to MANAGER_CLOSED.
!!

#include "f90_assert.fpp"

module viz_manager_type

  use iso_fortran_env, only: r8 => real64
  use viz_stream_state_type
  use viz_field_registry_type, only: viz_field_registry
  use viz_provider_class
  implicit none
  private

  integer, parameter :: MANAGER_CLOSED = 0
  integer, parameter :: MANAGER_INITIALIZED = 1
  integer, parameter :: MANAGER_MESH_WRITTEN = 2
  integer, parameter :: MANAGER_TEMPORAL_REGISTERED = 3

  type, public :: viz_manager
    private
    type(viz_stream_state), allocatable :: streams(:)
    type(viz_provider_box), allocatable :: providers(:)
    integer :: nprovider = 0
    integer :: state = MANAGER_CLOSED
  contains
    procedure :: init
    procedure :: schedule_events
    procedure :: open_and_write_mesh
    procedure :: register_temporal_data
    procedure :: write_timestep
    procedure :: write_stream_timestep
    procedure :: close
    procedure, private :: register_active_providers
    procedure, private :: append_provider
  end type

contains

  subroutine init(this, params, stat, errmsg)

    use parameter_list_type

    class(viz_manager), intent(out) :: this
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    type(parameter_list), pointer :: streams, stream
    type(parameter_list_iterator) :: piter

    streams => params%sublist('streams')
    piter = parameter_list_iterator(streams, sublists_only=.true.)

    n = piter%count()
    allocate(this%streams(n))

    n = 0
    do while (.not.piter%at_end())
      n = n + 1
      stream => piter%sublist()
      call this%streams(n)%init(stream, stat, errmsg)
      if (stat /= 0) return
      call piter%next
    end do

    this%state = MANAGER_INITIALIZED
    stat = 0

  end subroutine init

  subroutine schedule_events(this, event_queue, t_final)
    use sim_event_queue_type
    class(viz_manager), intent(inout) :: this
    class(sim_event_queue), intent(inout) :: event_queue
    real(r8), intent(out), optional :: t_final
    integer :: i
    real(r8) :: max_t_final
    INSIST(this%state >= MANAGER_INITIALIZED)
    max_t_final = -huge(1.0_r8)
    do i = 1, size(this%streams)
      call this%streams(i)%schedule_stream_events(event_queue, stream_id=i, t_final=t_final)
      if (present(t_final)) max_t_final = max(t_final, max_t_final)
    end do
    if (present(t_final)) t_final = max_t_final
  end subroutine

  subroutine open_and_write_mesh(this, mesh, stat, errmsg)
    use unstr_mesh_type
    class(viz_manager), intent(inout) :: this
    type(unstr_mesh), intent(in) :: mesh
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer :: i
    INSIST(this%state == MANAGER_INITIALIZED)
    do i = 1, size(this%streams)
      call this%streams(i)%open_and_write_mesh(mesh, stat, errmsg)
      if (stat /= 0) then
        call this%close()
        return
      end if
    end do
    this%state = MANAGER_MESH_WRITTEN
    stat = 0
  end subroutine

  subroutine register_temporal_data(this, catalog, stat, errmsg)
    class(viz_manager), intent(inout) :: this
    type(viz_provider_box), intent(in) :: catalog(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    type(viz_field_registry) :: registry
    integer :: i
    INSIST(this%state == MANAGER_MESH_WRITTEN)
    call this%register_active_providers(catalog, registry)
    do i = 1, size(this%streams)
      call this%streams(i)%configure_provider_states(registry, this%providers)
      call this%streams(i)%register_temporal_data()
    end do
    this%state = MANAGER_TEMPORAL_REGISTERED
    stat = 0
  end subroutine

  subroutine register_active_providers(this, catalog, registry)
    class(viz_manager), intent(inout) :: this
    type(viz_provider_box), intent(in) :: catalog(:)
    type(viz_field_registry), intent(inout) :: registry
    integer :: i
    do i = 1, size(catalog)
      associate (provider => catalog(i)%p)
        if (.not.provider%is_active()) cycle
        this%nprovider = this%nprovider + 1
        call provider%register_fields(registry, this%nprovider)
        call this%append_provider(provider)
      end associate
    end do
  end subroutine

  subroutine append_provider(this, provider)
    class(viz_manager), intent(inout) :: this
    class(viz_provider), intent(in) :: provider
    type(viz_provider_box), allocatable :: tmp(:)
    integer :: i, n
    if (allocated(this%providers)) then
      n = size(this%providers) + 1
      allocate(tmp(n))
      do i = 1, size(this%providers)
        call move_alloc(this%providers(i)%p, tmp(i)%p)
      end do
      call move_alloc(tmp, this%providers)
    else
      n = 1
      allocate(this%providers(n))
    end if
    allocate(this%providers(n)%p, source=provider)
  end subroutine

  subroutine write_timestep(this, t)
    class(viz_manager), intent(inout) :: this
    real(r8), intent(in) :: t
    integer :: i
    INSIST(this%state == MANAGER_TEMPORAL_REGISTERED)
    do i = 1, size(this%streams)
      call this%streams(i)%write_timestep(t)
    end do
  end subroutine

  subroutine write_stream_timestep(this, stream_id, t)
    class(viz_manager), intent(inout) :: this
    integer, intent(in) :: stream_id
    real(r8), intent(in) :: t
    INSIST(this%state == MANAGER_TEMPORAL_REGISTERED)
    INSIST(stream_id >= 1 .and. stream_id <= size(this%streams))
    call this%streams(stream_id)%write_timestep(t)
  end subroutine

  subroutine close(this)
    class(viz_manager), intent(inout) :: this
    call clear_runtime_state()
    call clear_configuration_state(this)
  contains
    subroutine clear_runtime_state()
      integer :: i
      if (allocated(this%streams)) then
        do i = 1, size(this%streams)
          call this%streams(i)%close
        end do
      end if
      if (allocated(this%providers)) deallocate(this%providers)
      this%nprovider = 0
    end subroutine
    subroutine clear_configuration_state(this)
      type(viz_manager), intent(out) :: this
    end subroutine
  end subroutine

end module viz_manager_type
