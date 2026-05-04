!!
!! EM_HEAT_VIZ_PROVIDER_TYPE
!!
!! Electromagnetic-heating adapter for the VTKHDF provider layer.
!!
!! This provider exposes electromagnetic heating output when EM heat is active.
!! The stream-local state records VTKHDF handles and the last EM heat generation
!! written for each block.  The EM heat driver remains responsible for computing
!! the heat source and exposing a read-only view of it.

module em_heat_viz_provider_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use em_heat_driver, only: em_heat_enabled, em_heat_generation, em_heat_ptr
  use viz_field_selection_type
  use viz_field_registry_type
  use viz_provider_class
  use viz_provider_state_class
  use vtkhdf_mb_file_type
  implicit none
  private

  type :: em_heat_viz_block_data
    type(vtkhdf_cell_data_handle) :: joule_power
    integer :: last_written_generation = -1
  end type

  type, extends(viz_provider), public :: em_heat_viz_provider
  contains
    procedure :: is_active
    procedure :: register_fields
    procedure :: create_state
  end type

  type, extends(viz_provider_state) :: em_heat_viz_provider_state
    private
    type(em_heat_viz_block_data), allocatable :: blocks(:)
    type(viz_field_selection) :: fields
    logical :: write_em_heat = .true.
  contains
    procedure :: register_mesh_block_temporal_data
    procedure :: write_mesh_block_timestep
    procedure, private :: write_em_heat_block
  end type

contains

  logical function is_active(this)
    class(em_heat_viz_provider), intent(in) :: this
    is_active = em_heat_enabled()
  end function

  subroutine register_fields(this, registry, provider_id)
    class(em_heat_viz_provider), intent(in) :: this
    type(viz_field_registry), intent(inout) :: registry
    integer, intent(in) :: provider_id
    if (.not.this%is_active()) return
    call registry%register_field('em_heat', provider_id)
  end subroutine

  subroutine create_state(this, nblock, state, field_names)
    class(em_heat_viz_provider), intent(in) :: this
    integer, intent(in) :: nblock
    class(viz_provider_state), allocatable, intent(out) :: state
    character(*), intent(in), optional :: field_names(:)
    type(em_heat_viz_provider_state), allocatable :: new_state
    allocate(new_state)
    call new_state%fields%init(field_names)
    new_state%write_em_heat = field_is_selected('em_heat')
    allocate(new_state%blocks(nblock))
    call move_alloc(new_state, state)
  contains
    logical function field_is_selected(name)
      character(*), intent(in) :: name
      field_is_selected = new_state%fields%write_all()
      if (.not.field_is_selected .and. new_state%fields%has_selected_fields()) &
          field_is_selected = any(new_state%fields%names == name)
    end function
  end subroutine

  subroutine register_mesh_block_temporal_data(this, file, iblock, block, block_cells, block_nodes)
    class(em_heat_viz_provider_state), intent(inout) :: this
    type(vtkhdf_mb_file), intent(inout) :: file
    integer, intent(in) :: iblock
    type(vtkhdf_block_handle), intent(in) :: block
    integer, intent(in) :: block_cells(:), block_nodes(:)
    if (.not.this%write_em_heat) return
    this%blocks(iblock)%joule_power = file%register_temporal_cell_data(block, 'em_heat', 0.0_r8)
    this%blocks(iblock)%last_written_generation = -1
  end subroutine

  subroutine write_mesh_block_timestep(this, file, iblock, block, block_cells, block_nodes)
    class(em_heat_viz_provider_state), intent(inout) :: this
    type(vtkhdf_mb_file), intent(inout) :: file
    integer, intent(in) :: iblock
    type(vtkhdf_block_handle), intent(in) :: block
    integer, intent(in) :: block_cells(:), block_nodes(:)
    if (.not.this%write_em_heat) return
    call this%write_em_heat_block(file, iblock, block, block_cells)
  end subroutine

  subroutine write_em_heat_block(this, file, iblock, block, block_cells)
    class(em_heat_viz_provider_state), intent(inout) :: this
    type(vtkhdf_mb_file), intent(inout) :: file
    integer, intent(in) :: iblock
    type(vtkhdf_block_handle), intent(in) :: block
    integer, intent(in) :: block_cells(:)
    integer :: generation
    real(r8), pointer :: q(:)
    generation = em_heat_generation()
    if (this%blocks(iblock)%last_written_generation == generation) return
    q => em_heat_ptr()
    call file%write_cell_data(block, this%blocks(iblock)%joule_power, q(block_cells))
    this%blocks(iblock)%last_written_generation = generation
  end subroutine

end module em_heat_viz_provider_type
