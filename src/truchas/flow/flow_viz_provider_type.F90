!!
!! FLOW_VIZ_PROVIDER_TYPE
!!
!! Flow package adapter for the VTKHDF provider layer.
!!
!! This provider advertises flow fields when flow is enabled and creates a
!! stream-local state object holding per-block VTKHDF handles plus the resolved
!! output policy.  The provider state owns the VTKHDF registration/write
!! logic and accesses the flow package only through its specific data views.
!!

module flow_viz_provider_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use viz_field_selection_type
  use viz_field_registry_type
  use viz_provider_class
  use viz_provider_state_class
  use vtkhdf_mb_file_type
  implicit none
  private

  type :: block_data
    type(vtkhdf_cell_data_handle) :: velocity
    type(vtkhdf_cell_data_handle) :: pressure
  end type

  type :: output_plan
    logical :: velocity = .false.
    logical :: pressure = .false.
  end type

  type, extends(viz_provider), public :: flow_viz_provider
  contains
    procedure :: is_active
    procedure :: register_fields
    procedure :: create_state
  end type

  type, extends(viz_provider_state), public :: flow_viz_provider_state
    private
    type(block_data), allocatable :: blocks(:)
    type(output_plan) :: output
  contains
    procedure :: register_mesh_block_temporal_data
    procedure :: write_mesh_block_timestep
  end type

contains

  logical function is_active(this)
    use flow_driver, only: flow_enabled
    class(flow_viz_provider), intent(in) :: this
    is_active = flow_enabled()
  end function

  subroutine register_fields(this, registry, provider_id)
    class(flow_viz_provider), intent(in) :: this
    type(viz_field_registry), intent(inout) :: registry
    integer, intent(in) :: provider_id
    if (.not.this%is_active()) return
    call registry%register_field('velocity', provider_id)
    call registry%register_field('pressure', provider_id)
  end subroutine

  subroutine create_state(this, nblock, state, field_names)
    class(flow_viz_provider), intent(in) :: this
    integer, intent(in) :: nblock
    class(viz_provider_state), allocatable, intent(out) :: state
    character(*), intent(in), optional :: field_names(:)
    type(viz_field_selection) :: fields
    type(flow_viz_provider_state), allocatable :: new_state
    allocate(new_state)
    call fields%init(field_names)
    call build_output_plan(fields, new_state%output)
    allocate(new_state%blocks(nblock))
    call move_alloc(new_state, state)
  end subroutine


  subroutine build_output_plan(fields, output)

    type(viz_field_selection), intent(in) :: fields
    type(output_plan), intent(out) :: output

    integer :: i

    if (fields%write_all()) then
      output%velocity = .true.
      output%pressure = .true.
      return
    end if

    if (.not.fields%has_selected_fields()) return

    do i = 1, size(fields%names)
      select case (fields%names(i))
      case ('velocity')
        output%velocity = .true.
      case ('pressure')
        output%pressure = .true.
      end select
    end do

  end subroutine build_output_plan


  subroutine register_mesh_block_temporal_data(this, file, iblock, block, block_cells, block_nodes)

    class(flow_viz_provider_state), intent(inout) :: this
    type(vtkhdf_mb_file), intent(inout) :: file
    integer, intent(in) :: iblock
    type(vtkhdf_block_handle), intent(in) :: block
    integer, intent(in) :: block_cells(:), block_nodes(:)

    real(r8) :: vector(3)

    associate (handle => this%blocks(iblock))
      if (this%output%velocity) &
          handle%velocity = file%register_temporal_cell_data(block, 'velocity', vector)
      if (this%output%pressure) &
          handle%pressure = file%register_temporal_cell_data(block, 'pressure', 0.0_r8)
    end associate

  end subroutine register_mesh_block_temporal_data


  subroutine write_mesh_block_timestep(this, file, iblock, block, block_cells, block_nodes)

    use flow_driver, only: flow_vel_cc_view, flow_P_cc_view

    class(flow_viz_provider_state), intent(inout) :: this
    type(vtkhdf_mb_file), intent(inout) :: file
    integer, intent(in) :: iblock
    type(vtkhdf_block_handle), intent(in) :: block
    integer, intent(in) :: block_cells(:), block_nodes(:)

    real(r8), pointer :: p(:), v(:,:)

    associate (handle => this%blocks(iblock))
      if (this%output%velocity) then
        v => flow_vel_cc_view()
        call file%write_cell_data(block, handle%velocity, v(:,block_cells))
      end if
      if (this%output%pressure) then
        p => flow_P_cc_view()
        call file%write_cell_data(block, handle%pressure, p(block_cells))
      end if
    end associate

  end subroutine write_mesh_block_timestep

end module flow_viz_provider_type
