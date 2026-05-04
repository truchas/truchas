!!
!! THERMAL_STATE_VIZ_PROVIDER_TYPE
!!
!! Thermal-state adapter for the VTKHDF provider layer.
!!
!! Temperature always exists as simulation state and enthalpy is always
!! available as derived cell state. This provider owns those fields for VTKHDF
!! output independently of whether the heat-transfer solver is active. For
!! fixed blocks with inactive heat transfer, the fields are written once as
!! static cell data. For moving blocks, VTKHDF requires temporal datasets, so
!! inactive heat transfer writes those datasets once on the first timestep.
!!

module thermal_state_viz_provider_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use diffusion_solver_data, only: heat_eqn
  use viz_field_selection_type
  use viz_field_variability_type, only: VIZ_FIELD_CONSTANT, VIZ_FIELD_DYNAMIC, &
      VIZ_FIELD_OUTPUT_STATIC, VIZ_FIELD_OUTPUT_TEMPORAL_CONSTANT, &
      VIZ_FIELD_OUTPUT_TEMPORAL_DYNAMIC, get_viz_field_output_mode
  use viz_field_registry_type
  use viz_provider_class
  use viz_provider_state_class
  use vtkhdf_mb_file_type
  implicit none
  private

  type :: block_data
    private
    type(vtkhdf_cell_data_handle) :: temperature
    type(vtkhdf_cell_data_handle) :: enthalpy
    integer :: temperature_output_mode = VIZ_FIELD_OUTPUT_STATIC
    integer :: enthalpy_output_mode = VIZ_FIELD_OUTPUT_STATIC
    logical :: temperature_written = .false.
    logical :: enthalpy_written = .false.
    logical :: requires_temporal_output = .false.
  end type

  type :: output_plan
    logical :: temperature = .false.
    logical :: enthalpy = .false.
  end type

  type, extends(viz_provider), public :: thermal_state_viz_provider
  contains
    procedure :: is_active
    procedure :: register_fields
    procedure :: create_state
  end type

  type, extends(viz_provider_state) :: thermal_state_viz_provider_state
    private
    type(block_data), allocatable :: blocks(:)
    type(output_plan) :: output
    logical :: heat_eqn_active = .false.
  contains
    procedure :: configure_blocks
    procedure :: register_mesh_block_temporal_data
    procedure :: write_mesh_block_timestep
  end type

contains

  logical function is_active(this)
    class(thermal_state_viz_provider), intent(in) :: this
    is_active = .true.
  end function

  subroutine register_fields(this, registry, provider_id)
    class(thermal_state_viz_provider), intent(in) :: this
    type(viz_field_registry), intent(inout) :: registry
    integer, intent(in) :: provider_id
    call registry%register_field('temperature', provider_id)
    call registry%register_field('enthalpy', provider_id)
  end subroutine

  subroutine create_state(this, nblock, state, field_names)
    class(thermal_state_viz_provider), intent(in) :: this
    integer, intent(in) :: nblock
    class(viz_provider_state), allocatable, intent(out) :: state
    character(*), intent(in), optional :: field_names(:)
    type(viz_field_selection) :: fields
    type(thermal_state_viz_provider_state), allocatable :: new_state
    allocate(new_state)
    call fields%init(field_names)
    call build_output_plan(fields, new_state%output)
    new_state%heat_eqn_active = heat_eqn
    allocate(new_state%blocks(nblock))
    call move_alloc(new_state, state)
  end subroutine

  subroutine build_output_plan(fields, output)

    type(viz_field_selection), intent(in) :: fields
    type(output_plan), intent(out) :: output

    integer :: i

    if (fields%write_all()) then
      output%temperature = .true.
      output%enthalpy = .true.
      return
    end if

    if (.not.fields%has_selected_fields()) return

    do i = 1, size(fields%names)
      select case (trim(fields%names(i)))
      case ('temperature')
        output%temperature = .true.
      case ('enthalpy')
        output%enthalpy = .true.
      end select
    end do

  end subroutine build_output_plan

  subroutine configure_blocks(this, block_requires_temporal_output)
    class(thermal_state_viz_provider_state), intent(inout) :: this
    logical, intent(in) :: block_requires_temporal_output(:)
    integer :: iblock
    do iblock = 1, size(this%blocks)
      this%blocks(iblock)%requires_temporal_output = block_requires_temporal_output(iblock)
    end do
  end subroutine

  subroutine register_mesh_block_temporal_data(this, file, iblock, block, block_cells, block_nodes)

    use zone_module, only: zone

    class(thermal_state_viz_provider_state), intent(inout) :: this
    type(vtkhdf_mb_file), intent(inout) :: file
    integer, intent(in) :: iblock
    type(vtkhdf_block_handle), intent(in) :: block
    integer, intent(in) :: block_cells(:), block_nodes(:)

    integer :: variability

    if (this%output%temperature) then
      variability = merge(VIZ_FIELD_DYNAMIC, VIZ_FIELD_CONSTANT, this%heat_eqn_active)
      this%blocks(iblock)%temperature_output_mode = get_viz_field_output_mode(variability, &
          this%blocks(iblock)%requires_temporal_output)
      select case (this%blocks(iblock)%temperature_output_mode)
      case (VIZ_FIELD_OUTPUT_STATIC)
        call file%write_cell_data(block, 'temperature', zone(block_cells)%temp)
        this%blocks(iblock)%temperature_written = .true.
      case default
        this%blocks(iblock)%temperature = file%register_temporal_cell_data(block, 'temperature', 0.0_r8)
      end select
    end if

    if (this%output%enthalpy) then
      variability = merge(VIZ_FIELD_DYNAMIC, VIZ_FIELD_CONSTANT, this%heat_eqn_active)
      this%blocks(iblock)%enthalpy_output_mode = get_viz_field_output_mode(variability, &
          this%blocks(iblock)%requires_temporal_output)
      select case (this%blocks(iblock)%enthalpy_output_mode)
      case (VIZ_FIELD_OUTPUT_STATIC)
        call file%write_cell_data(block, 'enthalpy', zone(block_cells)%enthalpy)
        this%blocks(iblock)%enthalpy_written = .true.
      case default
        this%blocks(iblock)%enthalpy = file%register_temporal_cell_data(block, 'enthalpy', 0.0_r8)
      end select
    end if

  end subroutine register_mesh_block_temporal_data

  subroutine write_mesh_block_timestep(this, file, iblock, block, block_cells, block_nodes)

    use zone_module, only: zone

    class(thermal_state_viz_provider_state), intent(inout) :: this
    type(vtkhdf_mb_file), intent(inout) :: file
    integer, intent(in) :: iblock
    type(vtkhdf_block_handle), intent(in) :: block
    integer, intent(in) :: block_cells(:), block_nodes(:)

    if (this%output%temperature) then
      select case (this%blocks(iblock)%temperature_output_mode)
      case (VIZ_FIELD_OUTPUT_TEMPORAL_CONSTANT)
        if (.not.this%blocks(iblock)%temperature_written) then
          call file%write_cell_data(block, this%blocks(iblock)%temperature, zone(block_cells)%temp)
          this%blocks(iblock)%temperature_written = .true.
        end if
      case (VIZ_FIELD_OUTPUT_TEMPORAL_DYNAMIC)
        call file%write_cell_data(block, this%blocks(iblock)%temperature, zone(block_cells)%temp)
        this%blocks(iblock)%temperature_written = .true.
      end select
    end if

    if (this%output%enthalpy) then
      select case (this%blocks(iblock)%enthalpy_output_mode)
      case (VIZ_FIELD_OUTPUT_TEMPORAL_CONSTANT)
        if (.not.this%blocks(iblock)%enthalpy_written) then
          call file%write_cell_data(block, this%blocks(iblock)%enthalpy, zone(block_cells)%enthalpy)
          this%blocks(iblock)%enthalpy_written = .true.
        end if
      case (VIZ_FIELD_OUTPUT_TEMPORAL_DYNAMIC)
        call file%write_cell_data(block, this%blocks(iblock)%enthalpy, zone(block_cells)%enthalpy)
        this%blocks(iblock)%enthalpy_written = .true.
      end select
    end if

  end subroutine write_mesh_block_timestep

end module thermal_state_viz_provider_type
