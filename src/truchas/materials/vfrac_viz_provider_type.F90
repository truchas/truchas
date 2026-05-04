!!
!! VFRAC_VIZ_PROVIDER_TYPE
!!
!! Volume-fraction material adapter for the VTKHDF provider layer.
!!
!! Volume-fraction fields are named from the configured material phases, so this provider
!! registers one selectable field per phase when the run has multiple
!! materials.  Its stream-local state resolves selected field names into
!! `(phase_id, field_name)` pairs once, then uses legacy material accessors to write the selected phase
!! fractions on each mesh block.

module vfrac_viz_provider_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use legacy_matl_api, only: gather_vof_cells, nmat
  use material_model_driver, only: matl_model
  use viz_field_registry_type
  use viz_provider_class
  use viz_provider_state_class
  use vtkhdf_mb_file_type
  implicit none
  private

  type :: block_data
    type(vtkhdf_cell_data_handle), allocatable :: phase(:)
  end type

  type :: output_field
    integer :: phase_id = 0
    character(:), allocatable :: field_name
  end type

  type, extends(viz_provider), public :: vfrac_viz_provider
  contains
    procedure :: is_active
    procedure :: register_fields
    procedure :: create_state
  end type

  type, extends(viz_provider_state), public :: vfrac_viz_provider_state
    private
    type(block_data), allocatable :: blocks(:)
    type(output_field), allocatable :: fields(:)
  contains
    procedure :: register_mesh_block_temporal_data
    procedure :: write_mesh_block_timestep
  end type

contains

  logical function is_active(this)
    class(vfrac_viz_provider), intent(in) :: this
    is_active = (nmat > 1)
  end function

  subroutine register_fields(this, registry, provider_id)
    class(vfrac_viz_provider), intent(in) :: this
    type(viz_field_registry), intent(inout) :: registry
    integer, intent(in) :: provider_id
    integer :: m
    if (.not.this%is_active()) return
    call registry%register_field('vfrac', provider_id)
    do m = 1, nmat
      call registry%register_field(vfrac_field_name(m), provider_id)
    end do
  end subroutine

  subroutine create_state(this, nblock, state, field_names)

    class(vfrac_viz_provider), intent(in) :: this
    integer, intent(in) :: nblock
    class(viz_provider_state), allocatable, intent(out) :: state
    character(*), intent(in), optional :: field_names(:)
    type(vfrac_viz_provider_state), allocatable :: new_state
    integer :: iblock

    allocate(new_state)
    call init_output_fields(new_state%fields, field_names)
    allocate(new_state%blocks(nblock))
    do iblock = 1, nblock
      allocate(new_state%blocks(iblock)%phase(size(new_state%fields)))
    end do
    call move_alloc(new_state, state)

  end subroutine create_state

  subroutine init_output_fields(fields, selected_field_names)

    integer :: i, m, n
    type(output_field), allocatable :: tmp(:)
    type(output_field), allocatable, intent(out) :: fields(:)
    character(*), intent(in), optional :: selected_field_names(:)

    if (.not.present(selected_field_names)) then
      allocate(fields(nmat))
      do m = 1, nmat
        fields(m)%phase_id = m
        fields(m)%field_name = vfrac_field_name(m)
      end do
      return
    end if

    allocate(tmp(size(selected_field_names)))
    n = 0
    do i = 1, size(selected_field_names)
      if (trim(selected_field_names(i)) == 'vfrac') then
        deallocate(tmp)
        allocate(fields(nmat))
        do m = 1, nmat
          fields(m)%phase_id = m
          fields(m)%field_name = vfrac_field_name(m)
        end do
        return
      end if
      do m = 1, nmat
        if (trim(selected_field_names(i)) == vfrac_field_name(m)) then
          n = n + 1
          tmp(n)%phase_id = m
          tmp(n)%field_name = trim(selected_field_names(i))
          exit
        end if
      end do
    end do

    allocate(fields(n))
    fields = tmp(:n)

  end subroutine init_output_fields

  function vfrac_field_name(m) result(name)
    integer, intent(in) :: m
    character(:), allocatable :: name
    name = 'vfrac_' // trim(matl_model%phase_name(m))
  end function

  subroutine register_mesh_block_temporal_data(this, file, iblock, block, block_cells, block_nodes)

    class(vfrac_viz_provider_state), intent(inout) :: this
    type(vtkhdf_mb_file), intent(inout) :: file
    integer, intent(in) :: iblock
    type(vtkhdf_block_handle), intent(in) :: block
    integer, intent(in) :: block_cells(:), block_nodes(:)

    integer :: i

    associate (handle => this%blocks(iblock))
      do i = 1, size(this%fields)
        handle%phase(i) = file%register_temporal_cell_data(block, this%fields(i)%field_name, 0.0_r8)
      end do
    end associate

  end subroutine register_mesh_block_temporal_data

  subroutine write_mesh_block_timestep(this, file, iblock, block, block_cells, block_nodes)

    class(vfrac_viz_provider_state), intent(inout) :: this
    type(vtkhdf_mb_file), intent(inout) :: file
    integer, intent(in) :: iblock
    type(vtkhdf_block_handle), intent(in) :: block
    integer, intent(in) :: block_cells(:), block_nodes(:)

    integer :: i
    real(r8), allocatable :: vof(:)

    allocate(vof(size(block_cells)))
    associate (handle => this%blocks(iblock))
      do i = 1, size(this%fields)
        ! NB: VOF data is only available for on-process cells. If mesh-block
        ! output ever includes ghost cells, this gather path must be revisited.
        call gather_vof_cells(this%fields(i)%phase_id, block_cells, vof)
        call file%write_cell_data(block, handle%phase(i), vof)
      end do
    end associate

  end subroutine write_mesh_block_timestep

end module vfrac_viz_provider_type
