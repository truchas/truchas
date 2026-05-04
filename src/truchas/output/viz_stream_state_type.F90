!!
!! VIZ_STREAM_STATE_TYPE
!!
!! Runtime state for one VTKHDF output stream.
!!
!! A stream owns one VTKHDF file, its output schedule, mesh-block metadata, any
!! moving-block state, and the provider states that write fields into that file.
!! It is responsible for file creation, mesh-block creation, temporal dataset
!! registration, timestep begin/end calls, moving-geometry updates, and dispatch
!! to each provider state.  It does not decide which providers are active; that
!! policy belongs to viz_manager_type. Given the active provider set and field
!! registry, it assembles the stream-local provider states needed for this
!! stream's configured field selection.
!!
!! Multiple streams are supported by this type, but the current namelist bridge
!! instantiates only the "main" stream.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NOTES
!!
!!   STREAM_CLOSED
!!     -> init
!!   STREAM_INITIALIZED
!!     -> schedule_stream_events
!!     -> open_and_write_mesh
!!   STREAM_MESH_WRITTEN
!!     -> schedule_stream_events
!!     -> configure_provider_states
!!     -> register_temporal_data
!!   STREAM_TEMPORAL_REGISTERED
!!     -> schedule_stream_events
!!     -> write_timestep (zero or more times)
!!
!! close is valid from any state and returns the stream to STREAM_CLOSED.
!!

#include "f90_assert.fpp"

module viz_stream_state_type

  use,intrinsic :: iso_fortran_env, only: int8, r8 => real64
  use parameter_list_type
  use toolpath_type
  use parallel_communication, only: comm, is_IOP, nPE, this_PE
  use truchas_logging_services
  use unstr_mesh_type
  use vtkhdf_mb_file_type
  use viz_field_registry_type
  use viz_provider_class
  use viz_provider_state_class
  implicit none
  private

  integer, parameter :: STREAM_CLOSED = 0
  integer, parameter :: STREAM_INITIALIZED = 1
  integer, parameter :: STREAM_MESH_WRITTEN = 2
  integer, parameter :: STREAM_TEMPORAL_REGISTERED = 3
  integer, parameter :: FIELD_SELECTION_DEFAULT = 1
  integer, parameter :: FIELD_SELECTION_ALL = 2
  integer, parameter :: FIELD_SELECTION_SELECTED = 3
  character(*), parameter :: PROCESS_IDS_FIELD_NAME = 'ProcessIds'
  character(*), parameter :: PROCESS_ID_INPUT_FIELD_NAME = 'process_id'

  type :: vtkhdf_block_data
    type(vtkhdf_block_handle) :: handle
    type(vtkhdf_cell_data_handle) :: process_ids
    type(vtkhdf_point_data_handle) :: global_point_ids
    integer, allocatable :: cells(:), nodes(:)
    real(r8), allocatable :: x0(:,:)
    logical :: moving = .false.
    logical :: process_ids_written = .false.
    logical :: global_point_ids_written = .false.
  end type

  type, public :: viz_stream_state
    private
    character(:), allocatable :: filename
    real(r8), allocatable :: times(:)
    integer, allocatable :: subintervals(:)
    character(:), allocatable :: field_names(:)
    integer :: field_selection_mode = FIELD_SELECTION_DEFAULT
    logical :: write_process_ids = .true.
    integer, allocatable :: part(:)
    type(toolpath), allocatable :: part_path
    type(vtkhdf_mb_file) :: file
    type(vtkhdf_block_data), allocatable :: block_data(:)
    type(viz_provider_state_box), allocatable :: providers(:)
    real(r8), allocatable :: last_moving_translation(:)
    real(r8) :: last_write_time = huge(1.0_r8)
    integer :: state = STREAM_CLOSED
  contains
    procedure :: init
    procedure :: configure_provider_states
    procedure :: schedule_stream_events
    procedure :: open_and_write_mesh
    procedure :: register_temporal_data
    procedure :: write_timestep
    procedure :: close
  end type viz_stream_state

contains

  subroutine init(this, params, stat, errmsg)

    class(viz_stream_state), intent(out) :: this
    type(parameter_list), intent(inout), target :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: plist

    call params%get('filename', this%filename, stat, errmsg)
    if (stat /= 0) return

    plist => params%sublist('schedule')
    call validate_schedule(plist, this%times, this%subintervals, stat, errmsg)
    if (stat /= 0) return

    call validate_move_blocks(params, this%part, this%part_path, stat, errmsg)
    if (stat /= 0) return

    call validate_field_selection(params, this%field_selection_mode, this%field_names, stat, errmsg)
    if (stat /= 0) return
    if (allocated(this%field_names)) then
      this%write_process_ids = stream_writes_process_ids(this%field_selection_mode, &
          this%field_names)
    else
      this%write_process_ids = stream_writes_process_ids(this%field_selection_mode)
    end if
    this%write_process_ids = this%write_process_ids .and. nPE > 1
    this%state = STREAM_INITIALIZED

  end subroutine init


  subroutine configure_provider_states(this, registry, providers)

    class(viz_stream_state), intent(inout) :: this
    type(viz_field_registry), intent(in) :: registry
    type(viz_provider_box), intent(in) :: providers(:)

    integer :: i, j, nblock, provider_id
    class(viz_provider_state), allocatable :: state
    integer, allocatable :: field_ids(:)
    character(:), allocatable :: missing_names(:)
    character(:), allocatable :: provider_field_names(:)
    type(viz_provider_field_selection), allocatable :: provider_fields(:)

    INSIST(this%state == STREAM_MESH_WRITTEN)
    INSIST(.not.allocated(this%providers))
    nblock = size(this%block_data)

    select case (this%field_selection_mode)
    case (FIELD_SELECTION_DEFAULT)
      call TLS_info('OUTPUTS FIELDS not specified; writing all available visualization fields. ' // &
          'Specify FIELDS = ''all'' to preserve this behavior in future versions.')
      call create_all_provider_states()
      return
    case (FIELD_SELECTION_ALL)
      call create_all_provider_states()
      return
    case (FIELD_SELECTION_SELECTED)
      continue
    case default
      INSIST(.false.)
    end select

    ! Resolve selected field names, warn on missing names, and group by the
    ! active provider array index recorded when each field was registered.
    call provider_selected_field_names(this%field_names, provider_field_names)
    call registry%resolve_selection(provider_field_names, field_ids, missing_names)
    do i = 1, size(missing_names)
      call TLS_warn('requested visualization field "' // trim(missing_names(i)) // &
          '" is unavailable and will be omitted')
    end do
    call registry%group_fields(field_ids, provider_fields)

    do j = 1, size(provider_fields)
      provider_id = provider_fields(j)%provider_id
      ! provider_id indexes PROVIDERS directly; it is not a stable provider id.
      call providers(provider_id)%p%create_state(nblock, state, &
          provider_fields(j)%field_names)
      call append_provider_state(state)
    end do

  contains

    subroutine create_all_provider_states()
      do j = 1, size(providers)
        call providers(j)%p%create_state(nblock, state)
        call append_provider_state(state)
      end do
    end subroutine

    subroutine append_provider_state(state)
      class(viz_provider_state), allocatable, intent(inout) :: state
      type(viz_provider_state_box), allocatable :: tmp(:)
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
      call state%configure_blocks(this%block_data%moving)
      call move_alloc(state, this%providers(n)%p)
    end subroutine

  end subroutine configure_provider_states


  subroutine validate_field_selection(plist, selection_mode, field_names, stat, errmsg)

    type(parameter_list), intent(inout) :: plist
    integer, intent(out) :: selection_mode
    character(:), allocatable, intent(out) :: field_names(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: i, nall

    stat = 0
    selection_mode = FIELD_SELECTION_DEFAULT

    if (.not.plist%is_parameter('fields')) return

    call plist%get('fields', field_names, stat, errmsg)
    if (stat /= 0) return

    nall = 0
    do i = 1, size(field_names)
      if (trim(field_names(i)) == 'all') nall = nall + 1
    end do

    if (nall > 0) then
      if (nall /= 1 .or. size(field_names) /= 1) then
        stat = 1
        errmsg = 'OUTPUTS FIELDS value ''all'' may not be combined with other field names'
        return
      end if
      selection_mode = FIELD_SELECTION_ALL
      deallocate(field_names)
      return
    end if

    selection_mode = FIELD_SELECTION_SELECTED

  end subroutine validate_field_selection


  logical function stream_writes_process_ids(selection_mode, field_names)

    integer, intent(in) :: selection_mode
    character(*), intent(in), optional :: field_names(:)

    select case (selection_mode)
    case (FIELD_SELECTION_DEFAULT, FIELD_SELECTION_ALL)
      stream_writes_process_ids = .true.
    case (FIELD_SELECTION_SELECTED)
      INSIST(present(field_names))
      stream_writes_process_ids = any_process_ids_field(field_names)
    case default
      INSIST(.false.)
    end select

  end function stream_writes_process_ids


  subroutine provider_selected_field_names(field_names, provider_field_names)

    character(*), intent(in) :: field_names(:)
    character(:), allocatable, intent(out) :: provider_field_names(:)

    integer :: i, n

    n = count(.not.is_process_ids_field(field_names))
    allocate(character(max_provider_field_name_len(field_names)) :: provider_field_names(n))
    n = 0
    do i = 1, size(field_names)
      if (is_process_ids_field(field_names(i))) cycle
      n = n + 1
      provider_field_names(n) = trim(field_names(i))
    end do

  end subroutine provider_selected_field_names


  pure logical function any_process_ids_field(field_names)

    character(*), intent(in) :: field_names(:)

    integer :: i

    any_process_ids_field = .false.
    do i = 1, size(field_names)
      if (is_process_ids_field(field_names(i))) then
        any_process_ids_field = .true.
        return
      end if
    end do

  end function any_process_ids_field


  pure elemental logical function is_process_ids_field(field_name)

    character(*), intent(in) :: field_name

    select case (trim(field_name))
    case (PROCESS_ID_INPUT_FIELD_NAME)
      is_process_ids_field = .true.
    case default
      is_process_ids_field = .false.
    end select

  end function is_process_ids_field


  pure integer function max_provider_field_name_len(field_names) result(len)

    character(*), intent(in) :: field_names(:)

    integer :: i

    len = 0
    do i = 1, size(field_names)
      if (is_process_ids_field(field_names(i))) cycle
      len = max(len, len_trim(field_names(i)))
    end do

  end function max_provider_field_name_len


  subroutine validate_schedule(plist, times, subintervals, stat, errmsg)

    type(parameter_list), intent(inout) :: plist
    real(r8), allocatable, intent(out) :: times(:)
    integer, allocatable, intent(out) :: subintervals(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n, j

    call plist%get('times', times, stat, errmsg)
    if (stat /= 0) return

    n = size(times)
    if (n < 2) then
      stat = 1
      errmsg = 'times must have length at least 2'
      return
    end if

    do j = 1, n-1
      if (times(j+1) <= times(j)) then
        stat = 1
        errmsg = 'times must be strictly increasing'
        return
      end if
    end do

    if (plist%is_parameter('subintervals')) then
      call plist%get('subintervals', subintervals, stat, errmsg)
      if (stat /= 0) return
      if (size(subintervals) /= n-1) then
        stat = 1
        errmsg = 'size(subintervals) must equal size(times)-1'
        return
      end if
      do j = 1, n-1
        if (subintervals(j) <= 0) then
          stat = 1
          errmsg = 'subintervals must be positive'
          return
        end if
      end do
    else
      allocate(subintervals(n-1), source=1)
    end if

  end subroutine validate_schedule


  subroutine validate_move_blocks(plist, part, part_path, stat, errmsg)

    use toolpath_driver, only: alloc_toolpath

    type(parameter_list), intent(inout) :: plist
    integer, allocatable, intent(out) :: part(:)
    type(toolpath), allocatable, intent(out) :: part_path
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    character(:), allocatable :: string

    if (.not.plist%is_parameter('move-block-ids')) then
      if (plist%is_parameter('move-toolpath-name')) then
        stat = 1
        errmsg = 'move-toolpath-name requires move-block-ids'
      else
        stat = 0
      end if
      return
    end if

    call plist%get('move-block-ids', part, stat, errmsg)
    if (stat /= 0) return
    if (.not.plist%is_parameter('move-toolpath-name')) then
      stat = 1
      errmsg = 'move-toolpath-name not specified'
      return
    end if

    call plist%get('move-toolpath-name', string, stat, errmsg)
    if (stat /= 0) return
    call alloc_toolpath(part_path, string, stat, errmsg)
    if (stat /= 0) then
      errmsg = 'unable to create toolpath: ' // errmsg
      return
    end if

  end subroutine validate_move_blocks


  subroutine schedule_stream_events(this, event_queue, stream_id, t_final)

    use sim_event_queue_type
    use viz_output_action_type

    class(viz_stream_state), intent(in) :: this
    class(sim_event_queue),  intent(inout) :: event_queue
    integer, intent(in) :: stream_id
    real(r8), intent(out), optional :: t_final

    integer :: i, j
    real(r8) :: t0, t1, dt

    INSIST(this%state >= STREAM_INITIALIZED)
    do j = 1, size(this%subintervals)
      t0 = this%times(j)
      t1 = this%times(j+1)
      dt = (t1 - t0) / this%subintervals(j)
      do i = 0, this%subintervals(j)-1
        call event_queue%add_event(t0 + i*dt, viz_output_action(stream_id=stream_id))
      end do
    end do
    call event_queue%add_event(t1, viz_output_action(stream_id=stream_id))

    if (present(t_final)) t_final = t1

  end subroutine schedule_stream_events


  subroutine open_and_write_mesh(this, mesh, stat, errmsg)

    use string_utilities, only: i_to_c

    class(viz_stream_state), intent(inout) :: this
    type(unstr_mesh), intent(in) :: mesh
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(vtkhdf_block_handle) :: block
    real(r8), allocatable :: x(:,:)
    integer, allocatable :: xcnode(:), cnode(:), block_cells(:), block_nodes(:)
    integer(int8), allocatable :: types(:)
    integer(kind(mesh%cell_set_mask)) :: bitmask
    integer :: n
    logical :: moving

    INSIST(this%state == STREAM_INITIALIZED)
    if (allocated(this%part)) then
      do n = 1, size(this%part)
        if (any(this%part(n) == mesh%cell_set_id)) cycle
        stat = 1
        errmsg = 'move-block-ids contains unknown mesh block id ' // i_to_c(this%part(n))
        return
      end do
    end if

    call this%file%create(this%filename, comm, stat, errmsg)
    if (stat /= 0) return

    allocate(this%block_data(size(mesh%cell_set_name)))
    do n = 1, size(mesh%cell_set_name)
      associate (name => mesh%cell_set_name(n)%s)
        moving = .false.
        if (allocated(this%part)) moving = any(this%part == mesh%cell_set_id(n))
        bitmask = 0
        bitmask = ibset(bitmask, pos=n)
        call mesh%get_mesh_block(bitmask, x, xcnode, cnode, &
            block_cells, block_nodes)
        types = get_vtk_cell_types(xcnode)
        if (moving) then
          block = this%file%add_block(name, mode=UG_MOVING_MESH)
          call this%file%write_mesh_topology(block, cnode, xcnode, types)
          call move_alloc(x, this%block_data(n)%x0)
        else
          block = this%file%add_block(name, mode=UG_FIXED_MESH)
          call this%file%write_mesh(block, x, cnode, xcnode, types)
          if (this%write_process_ids) &
              call this%file%write_cell_data(block, PROCESS_IDS_FIELD_NAME, &
                  spread(this_PE, dim=1, ncopies=size(block_cells)))
          call this%file%write_point_data(block, 'vtkGlobalPointIds', &
              mesh%node_imap%global_index(block_nodes))
        end if
        this%block_data(n)%handle = block
        this%block_data(n)%moving = moving
        this%block_data(n)%process_ids_written = .false.
        this%block_data(n)%global_point_ids_written = .false.
        call move_alloc(block_cells, this%block_data(n)%cells)
        call move_alloc(block_nodes, this%block_data(n)%nodes)
      end associate
    end do

    if (allocated(this%providers)) then
      do n = 1, size(this%providers)
        call this%providers(n)%p%configure_blocks(this%block_data%moving)
      end do
    end if

    this%state = STREAM_MESH_WRITTEN
    if (allocated(this%last_moving_translation)) deallocate(this%last_moving_translation)
    stat = 0

  end subroutine open_and_write_mesh


  subroutine register_temporal_data(this)

    class(viz_stream_state), intent(inout) :: this
    integer :: n, p

    INSIST(this%state == STREAM_MESH_WRITTEN)

    do n = 1, size(this%block_data)
      associate (b => this%block_data(n), block => this%block_data(n)%handle)
        if (b%moving) then
          if (this%write_process_ids) &
              b%process_ids = this%file%register_temporal_cell_data(block, &
                  PROCESS_IDS_FIELD_NAME, 0)
          b%global_point_ids = this%file%register_temporal_point_data(block, &
              'vtkGlobalPointIds', 0)
        end if
        if (allocated(this%providers)) then
          do p = 1, size(this%providers)
            call this%providers(p)%p%register_mesh_block_temporal_data(this%file, n, &
                block, b%cells, b%nodes)
          end do
        end if
      end associate
    end do
    if (allocated(this%providers)) then
      do p = 1, size(this%providers)
        call this%providers(p)%p%define_provider_blocks(this%file)
      end do
    end if
    this%state = STREAM_TEMPORAL_REGISTERED

  end subroutine register_temporal_data


  subroutine write_timestep(this, t)

    use mesh_manager, only: unstr_mesh_ptr

    class(viz_stream_state), intent(inout) :: this
    real(r8), intent(in) :: t
    integer :: n, p
    real(r8) :: r(3)
    real(r8), allocatable :: x(:,:)
    type(unstr_mesh), pointer :: mesh
    logical :: write_moving_geometry

    INSIST(this%state == STREAM_TEMPORAL_REGISTERED)
    if (t == this%last_write_time) return

    write_moving_geometry = .false.
    if (any(this%block_data%moving)) then
      INSIST(allocated(this%part_path))
      call this%part_path%set_segment(t)
      call this%part_path%get_position(t, r)
      if (allocated(this%last_moving_translation)) then
        write_moving_geometry = any(r /= this%last_moving_translation)
      else
        write_moving_geometry = .true.
      end if
      mesh => unstr_mesh_ptr('MAIN')
      INSIST(associated(mesh))
    end if

    call this%file%start_time_step(t)
    do n = 1, size(this%block_data)
      associate (b => this%block_data(n))
        if (b%moving) then
          if (write_moving_geometry) then
            x = b%x0
            x(1,:) = x(1,:) + r(1)
            x(2,:) = x(2,:) + r(2)
            x(3,:) = x(3,:) + r(3)
            call this%file%write_mesh_geometry(b%handle, x)
          end if
          if (this%write_process_ids .and. .not.b%process_ids_written) then
            call this%file%write_cell_data(b%handle, b%process_ids, &
                spread(this_PE, dim=1, ncopies=size(b%cells)))
            b%process_ids_written = .true.
          end if
          if (.not.b%global_point_ids_written) then
            call this%file%write_point_data(b%handle, b%global_point_ids, &
                mesh%node_imap%global_index(b%nodes))
            b%global_point_ids_written = .true.
          end if
        end if
      end associate
    end do
    if (any(this%block_data%moving) .and. write_moving_geometry) then
      if (.not.allocated(this%last_moving_translation)) &
          allocate(this%last_moving_translation(3))
      this%last_moving_translation = r
    end if
    do n = 1, size(this%block_data)
      associate (b => this%block_data(n))
        if (allocated(this%providers)) then
          do p = 1, size(this%providers)
            call this%providers(p)%p%write_mesh_block_timestep(this%file, n, b%handle, &
                b%cells, b%nodes)
          end do
        end if
      end associate
    end do
    if (allocated(this%providers)) then
      do p = 1, size(this%providers)
        call this%providers(p)%p%write_provider_blocks_timestep(this%file)
      end do
    end if
    call this%file%finalize_time_step()
    call this%file%flush()
    this%last_write_time = t

  end subroutine write_timestep

  ! Today close is a full teardown. The runtime/configuration split is kept
  ! explicit so a future continuation-file path can retain configuration
  ! while only clearing the runtime state associated with the current file.
  subroutine close(this)
    class(viz_stream_state), intent(inout) :: this
    call clear_runtime_state()
    call clear_configuration_state(this)
  contains
    subroutine clear_runtime_state()
      call this%file%close()
      if (allocated(this%block_data)) deallocate(this%block_data)
      if (allocated(this%providers)) deallocate(this%providers)
      if (allocated(this%last_moving_translation)) &
          deallocate(this%last_moving_translation)
    end subroutine
    subroutine clear_configuration_state(this)
      type(viz_stream_state), intent(out) :: this
    end subroutine
  end subroutine close

  function get_vtk_cell_types(xcnode) result(types)
    use vtkhdf_vtk_cell_types
    integer, intent(in) :: xcnode(:)
    integer(int8), allocatable :: types(:)
    integer :: j, nnode
    allocate(types(size(xcnode)-1))
    do j = 1, size(types)
      nnode = xcnode(j+1) - xcnode(j)
      select case (nnode)
      case (4)
        types(j) = VTK_TETRA
      case (5)
        types(j) = VTK_PYRAMID
      case (6)
        types(j) = VTK_WEDGE
      case (8)
        types(j) = VTK_HEXAHEDRON
      case default
        INSIST(.false.)
      end select
    end do
  end function get_vtk_cell_types

end module viz_stream_state_type
