!!
!! DIFFUSION_SOLVER_VIZ_PROVIDER_TYPE
!!
!! Diffusion/species package adapter for the VTKHDF provider layer.
!!
!! This provider registers species field names available for the current
!! diffusion-solver configuration. Its stream-local state owns per-block
!! VTKHDF handles, the resolved output policy, and the package-specific
!! species VTKHDF registration/write logic.
!!

module diffusion_solver_viz_provider_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use viz_field_selection_type
  use viz_field_registry_type
  use viz_provider_class
  use viz_provider_state_class
  use vtkhdf_mb_file_type
  use diffusion_solver, only: ds_get_phi_view
  use diffusion_solver_data, only: num_species
  use string_utilities, only: i_to_c
  implicit none
  private

  type :: block_data
    private
    type(vtkhdf_cell_data_handle), allocatable :: species(:)
  end type

  type :: output_plan
    logical, allocatable :: species(:)
  end type

  type, extends(viz_provider), public :: diffusion_solver_viz_provider
  contains
    procedure :: is_active
    procedure :: register_fields
    procedure :: create_state
  end type

  type, extends(viz_provider_state) :: diffusion_solver_viz_provider_state
    private
    type(block_data), allocatable :: blocks(:)
    type(output_plan) :: output
  contains
    procedure :: register_mesh_block_temporal_data
    procedure :: write_mesh_block_timestep
  end type

contains

  logical function is_active(this)
    class(diffusion_solver_viz_provider), intent(in) :: this
    is_active = (num_species > 0)
  end function

  subroutine register_fields(this, registry, provider_id)
    class(diffusion_solver_viz_provider), intent(in) :: this
    type(viz_field_registry), intent(inout) :: registry
    integer, intent(in) :: provider_id
    integer :: m
    if (.not.this%is_active()) return
    do m = 1, num_species
      call registry%register_field('phi_'//i_to_c(m), provider_id)
    end do
  end subroutine

  subroutine create_state(this, nblock, state, field_names)
    class(diffusion_solver_viz_provider), intent(in) :: this
    integer, intent(in) :: nblock
    class(viz_provider_state), allocatable, intent(out) :: state
    character(*), intent(in), optional :: field_names(:)
    type(viz_field_selection) :: fields
    type(diffusion_solver_viz_provider_state), allocatable :: new_state
    integer :: iblock
    allocate(new_state)
    call fields%init(field_names)
    call build_output_plan(fields, new_state%output)
    allocate(new_state%blocks(nblock))
    do iblock = 1, nblock
      allocate(new_state%blocks(iblock)%species(num_species))
    end do
    call move_alloc(new_state, state)
  end subroutine

  subroutine build_output_plan(fields, output)

    type(viz_field_selection), intent(in) :: fields
    type(output_plan), intent(out) :: output

    integer :: i, m
    character(:), allocatable :: name

    allocate(output%species(num_species), source=.false.)

    if (fields%write_all()) then
      output%species = .true.
      return
    end if

    if (.not.fields%has_selected_fields()) return

    do i = 1, size(fields%names)
      name = trim(fields%names(i))

      do m = 1, num_species
        if (name == 'phi_'//i_to_c(m)) then
          output%species(m) = .true.
          exit
        end if
      end do
    end do

  end subroutine build_output_plan

  subroutine register_mesh_block_temporal_data(this, file, iblock, block, block_cells, block_nodes)

    class(diffusion_solver_viz_provider_state), intent(inout) :: this
    type(vtkhdf_mb_file), intent(inout) :: file
    integer, intent(in) :: iblock
    type(vtkhdf_block_handle), intent(in) :: block
    integer, intent(in) :: block_cells(:), block_nodes(:)

    integer :: m

    associate (handle => this%blocks(iblock))
      do m = 1, size(handle%species)
        if (.not.this%output%species(m)) cycle
        handle%species(m) = file%register_temporal_cell_data(block, 'phi_'//i_to_c(m), 0.0_r8)
      end do
    end associate

  end subroutine

  subroutine write_mesh_block_timestep(this, file, iblock, block, block_cells, block_nodes)

    class(diffusion_solver_viz_provider_state), intent(inout) :: this
    type(vtkhdf_mb_file), intent(inout) :: file
    integer, intent(in) :: iblock
    type(vtkhdf_block_handle), intent(in) :: block
    integer, intent(in) :: block_cells(:), block_nodes(:)

    integer :: m
    real(r8), pointer :: phi(:)

    associate (handle => this%blocks(iblock))
      do m = 1, size(handle%species)
        if (.not.this%output%species(m)) cycle
        call ds_get_phi_view(m, phi)
        call file%write_cell_data(block, handle%species(m), phi(block_cells))
      end do
    end associate

  end subroutine

end module diffusion_solver_viz_provider_type
