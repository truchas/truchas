!!
!! SOLID_MECHANICS_VIZ_PROVIDER_TYPE
!!
!! Solid-mechanics adapter for the VTKHDF provider layer.
!!
!! This provider registers solid-mechanics volume fields and contact/interface
!! fields when the package is active.  Its stream-local state stores per-block
!! dataset handles, interface output handles, and a cached output plan derived
!! from the selected field names.  All VTKHDF-specific registration and write
!! logic lives in this provider module, with data acquired from the
!! solid-mechanics driver.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module solid_mechanics_viz_provider_type

  use,intrinsic :: iso_fortran_env, only: int8, r8 => real64
  use bitfield_type, only: bitfield
  use solid_mechanics_driver, only: solid_mechanics_enabled, &
      solid_mechanics_get_point_fields, solid_mechanics_get_cell_fields, &
      solid_mechanics_get_contact_set_ids, solid_mechanics_point_field_request, &
      solid_mechanics_point_field_values, solid_mechanics_cell_field_request, &
      solid_mechanics_cell_field_values, solid_mechanics_viscoplasticity_enabled
  use vtkhdf_mb_file_type, only: vtkhdf_mb_file, vtkhdf_block_handle, &
      vtkhdf_cell_data_handle, vtkhdf_point_data_handle, UG_FIXED_MESH
  use viz_field_selection_type
  use viz_field_registry_type
  use viz_provider_class
  use viz_provider_state_class
  implicit none
  private

  type :: block_data
    type(vtkhdf_point_data_handle) :: displacement
    type(vtkhdf_cell_data_handle) :: strain
    type(vtkhdf_cell_data_handle) :: thermal_strain
    type(vtkhdf_cell_data_handle) :: stress
    type(vtkhdf_cell_data_handle) :: rotation
    type(vtkhdf_cell_data_handle) :: plastic_strain
    type(vtkhdf_cell_data_handle) :: plastic_strain_rate
  end type

  type :: interface_block_data
    integer :: setid = 0
    type(vtkhdf_block_handle) :: handle
    type(vtkhdf_point_data_handle) :: gap_displacement
    type(vtkhdf_point_data_handle) :: gap_normal_traction
    integer, allocatable :: nodes(:)
  end type

  type :: output_data
    type(interface_block_data), allocatable :: interface(:)
  end type

  type :: output_plan
    logical :: displacement         = .false.
    logical :: strain               = .false.
    logical :: thermal_strain       = .false.
    logical :: stress               = .false.
    logical :: rotation             = .false.
    logical :: plastic_strain       = .false.
    logical :: plastic_strain_rate  = .false.
    logical :: gap_displacement     = .false.
    logical :: gap_normal_traction  = .false.
  contains
    procedure :: needs_any_cell_output
    procedure :: needs_interface_output
  end type

  type, extends(viz_provider), public :: solid_mechanics_viz_provider
  contains
    procedure :: is_active
    procedure :: register_fields
    procedure :: create_state
  end type

  type, extends(viz_provider_state), public :: solid_mechanics_viz_provider_state
    private
    type(block_data), allocatable :: blocks(:)
    type(output_data) :: output
    type(output_plan) :: output_plan
  contains
    procedure :: register_mesh_block_temporal_data
    procedure :: define_provider_blocks
    procedure :: write_mesh_block_timestep
    procedure :: write_provider_blocks_timestep
  end type

contains

  logical function is_active(this)
    class(solid_mechanics_viz_provider), intent(in) :: this
    is_active = solid_mechanics_enabled()
  end function

  subroutine register_fields(this, registry, provider_id)

    class(solid_mechanics_viz_provider), intent(in) :: this
    type(viz_field_registry), intent(inout) :: registry
    integer, intent(in) :: provider_id

    if (.not.this%is_active()) return

    call registry%register_field('displacement', provider_id)
    call registry%register_field('strain', provider_id)
    call registry%register_field('thermal_strain', provider_id)
    call registry%register_field('stress', provider_id)
    call registry%register_field('rotation', provider_id)
    call registry%register_field('contact_normal_displ', provider_id)
    call registry%register_field('contact_pressure', provider_id)
    if (solid_mechanics_viscoplasticity_enabled()) then
      call registry%register_field('plastic_strain', provider_id)
      call registry%register_field('plastic_strain_rate', provider_id)
    end if

  end subroutine register_fields

  subroutine create_state(this, nblock, state, field_names)

    class(solid_mechanics_viz_provider), intent(in) :: this
    integer, intent(in) :: nblock
    class(viz_provider_state), allocatable, intent(out) :: state
    character(*), intent(in), optional :: field_names(:)
    type(viz_field_selection) :: fields
    type(solid_mechanics_viz_provider_state), allocatable :: new_state

    allocate(new_state)
    call fields%init(field_names)
    call build_output_plan(fields, new_state%output_plan)
    allocate(new_state%blocks(nblock))
    call move_alloc(new_state, state)

  end subroutine create_state

  subroutine register_mesh_block_temporal_data(this, file, iblock, block, block_cells, &
      block_nodes)

    class(solid_mechanics_viz_provider_state), intent(inout) :: this
    type(vtkhdf_mb_file), intent(inout) :: file
    integer, intent(in) :: iblock
    type(vtkhdf_block_handle), intent(in) :: block
    integer, intent(in) :: block_cells(:), block_nodes(:)
    real(r8) :: vector(3), tensor(6)

    associate (handle => this%blocks(iblock))
      if (this%output_plan%displacement) &
          handle%displacement = file%register_temporal_point_data(block, 'displacement', vector)
      if (this%output_plan%strain) &
          handle%strain = file%register_temporal_cell_data(block, 'strain', tensor)
      if (this%output_plan%thermal_strain) &
          handle%thermal_strain = file%register_temporal_cell_data(block, 'thermal_strain', tensor)
      if (this%output_plan%stress) &
          handle%stress = file%register_temporal_cell_data(block, 'stress', tensor)
      if (this%output_plan%rotation) &
          handle%rotation = file%register_temporal_cell_data(block, 'rotation', 0.0_r8)
      if (this%output_plan%plastic_strain) &
          handle%plastic_strain = file%register_temporal_cell_data(block, 'plastic_strain', tensor)
      if (this%output_plan%plastic_strain_rate) &
          handle%plastic_strain_rate = file%register_temporal_cell_data(block, &
              'plastic_strain_rate', 0.0_r8)
    end associate

  end subroutine register_mesh_block_temporal_data

  subroutine write_mesh_block_timestep(this, file, iblock, block, block_cells, block_nodes)

    class(solid_mechanics_viz_provider_state), intent(inout) :: this
    type(vtkhdf_mb_file), intent(inout) :: file
    integer, intent(in) :: iblock
    type(vtkhdf_block_handle), intent(in) :: block
    integer, intent(in) :: block_cells(:), block_nodes(:)
    type(solid_mechanics_point_field_request) :: point_request
    type(solid_mechanics_point_field_values) :: point_values
    type(solid_mechanics_cell_field_request) :: cell_request
    type(solid_mechanics_cell_field_values) :: cell_values

    if (this%output_plan%displacement) then
      point_request%displacement = .true.
      call solid_mechanics_get_point_fields(block_nodes, point_request, point_values)
    end if

    if (this%output_plan%needs_any_cell_output()) then
      cell_request%total_strain = this%output_plan%strain
      cell_request%thermal_strain = this%output_plan%thermal_strain
      cell_request%stress = this%output_plan%stress
      cell_request%rotation = this%output_plan%rotation
      cell_request%plastic_strain = this%output_plan%plastic_strain
      cell_request%plastic_strain_rate = this%output_plan%plastic_strain_rate
      call solid_mechanics_get_cell_fields(block_cells, cell_request, cell_values)
    end if

    if (allocated(cell_values%total_strain)) &
      call to_vtk_symmetric_tensor_order(cell_values%total_strain)
    if (allocated(cell_values%thermal_strain)) &
      call to_vtk_symmetric_tensor_order(cell_values%thermal_strain)
    if (allocated(cell_values%stress)) &
      call to_vtk_symmetric_tensor_order(cell_values%stress)
    if (allocated(cell_values%plastic_strain)) &
      call to_vtk_symmetric_tensor_order(cell_values%plastic_strain)

    associate (handle => this%blocks(iblock))
      if (this%output_plan%displacement) &
        call file%write_point_data(block, handle%displacement, point_values%displacement)
      if (this%output_plan%strain) &
        call file%write_cell_data(block, handle%strain, cell_values%total_strain)
      if (this%output_plan%thermal_strain) &
        call file%write_cell_data(block, handle%thermal_strain, cell_values%thermal_strain)
      if (this%output_plan%stress) &
        call file%write_cell_data(block, handle%stress, cell_values%stress)
      if (this%output_plan%rotation) &
        call file%write_cell_data(block, handle%rotation, cell_values%rotation)
      if (this%output_plan%plastic_strain) &
        call file%write_cell_data(block, handle%plastic_strain, cell_values%plastic_strain)
      if (this%output_plan%plastic_strain_rate) &
        call file%write_cell_data(block, handle%plastic_strain_rate, &
            cell_values%plastic_strain_rate)
    end associate

  end subroutine write_mesh_block_timestep

  ! Provider-defined gap/interface block output hooks.

  subroutine define_provider_blocks(this, file)

    use mesh_manager, only: unstr_mesh_ptr
    use string_utilities, only: i_to_c
    use unstr_mesh_type, only: unstr_mesh

    class(solid_mechanics_viz_provider_state), intent(inout) :: this
    type(vtkhdf_mb_file), intent(inout) :: file
    type(unstr_mesh), pointer :: mesh
    type(vtkhdf_block_handle) :: block
    type(bitfield) :: bitmask
    real(r8), allocatable :: x(:,:)
    integer, allocatable :: setids(:), xfnode(:), fnode(:), faces(:), nodes(:)
    integer(int8), allocatable :: types(:)
    integer :: n, stat
    character(:), allocatable :: errmsg

    if (.not.this%output_plan%needs_interface_output()) return

    call solid_mechanics_get_contact_set_ids(setids)
    if (size(setids) == 0) return

    mesh => unstr_mesh_ptr('MAIN')
    INSIST(associated(mesh))

    allocate(this%output%interface(size(setids)))
    do n = 1, size(setids)
      call mesh%get_link_set_bitmask([setids(n)], bitmask, stat, errmsg)
      INSIST(stat == 0)
      call mesh%get_link_set_mesh(bitmask, x, xfnode, fnode, faces, nodes)
      types = get_vtk_face_types(xfnode)
      block = file%add_block('gap-interface-' // i_to_c(setids(n)), mode=UG_FIXED_MESH)
      call file%write_mesh(block, x, fnode, xfnode, types)
      this%output%interface(n)%setid = setids(n)
      this%output%interface(n)%handle = block
      if (this%output_plan%gap_displacement) &
          this%output%interface(n)%gap_displacement = file%register_temporal_point_data(block, &
              'contact_normal_displ', 0.0_r8)
      if (this%output_plan%gap_normal_traction) &
          this%output%interface(n)%gap_normal_traction = file%register_temporal_point_data(block, &
              'contact_pressure', 0.0_r8)
      call move_alloc(nodes, this%output%interface(n)%nodes)
    end do

  end subroutine define_provider_blocks

  subroutine write_provider_blocks_timestep(this, file)

    class(solid_mechanics_viz_provider_state), intent(inout) :: this
    type(vtkhdf_mb_file), intent(inout) :: file
    type(solid_mechanics_point_field_request) :: request
    type(solid_mechanics_point_field_values) :: values
    integer :: n

    if (.not.this%output_plan%needs_interface_output()) return
    if (.not.allocated(this%output%interface)) return

    request%gap_displacement = this%output_plan%gap_displacement
    request%gap_normal_traction = this%output_plan%gap_normal_traction
    do n = 1, size(this%output%interface)
      call solid_mechanics_get_point_fields(this%output%interface(n)%nodes, request, values)
      if (this%output_plan%gap_displacement) &
          call file%write_point_data(this%output%interface(n)%handle, &
              this%output%interface(n)%gap_displacement, values%gap_displacement)
      if (this%output_plan%gap_normal_traction) &
          call file%write_point_data(this%output%interface(n)%handle, &
              this%output%interface(n)%gap_normal_traction, values%gap_normal_traction)
    end do

  end subroutine write_provider_blocks_timestep

  subroutine to_vtk_symmetric_tensor_order(tensor)

    real(r8), intent(inout) :: tensor(:,:)
    real(r8) :: tmp(size(tensor, dim=2))

    INSIST(size(tensor, dim=1) == 6)
    tmp = tensor(5,:)
    tensor(5,:) = tensor(6,:)
    tensor(6,:) = tmp

  end subroutine to_vtk_symmetric_tensor_order

  function get_vtk_face_types(xfnode) result(types)

    use vtkhdf_vtk_cell_types, only: VTK_TRIANGLE, VTK_QUAD

    integer, intent(in) :: xfnode(:)
    integer(int8), allocatable :: types(:)
    integer :: j, nnode

    allocate(types(size(xfnode)-1))
    do j = 1, size(types)
      nnode = xfnode(j+1) - xfnode(j)
      select case (nnode)
      case (3)
        types(j) = VTK_TRIANGLE
      case (4)
        types(j) = VTK_QUAD
      case default
        INSIST(.false.)
      end select
    end do

  end function get_vtk_face_types

  logical function needs_any_cell_output(this)

    class(output_plan), intent(in) :: this

    needs_any_cell_output = this%strain              .or. this%thermal_strain .or. &
                            this%stress              .or. this%rotation       .or. &
                            this%plastic_strain      .or. &
                            this%plastic_strain_rate

  end function needs_any_cell_output

  logical function needs_interface_output(this)

    class(output_plan), intent(in) :: this

    needs_interface_output = this%gap_displacement .or. this%gap_normal_traction

  end function needs_interface_output

  subroutine build_output_plan(fields, plan)

    type(viz_field_selection), intent(in) :: fields
    type(output_plan), intent(out) :: plan

    integer :: i
    character(:), allocatable :: name

    if (fields%write_all()) then
      plan%displacement = .true.
      plan%strain = .true.
      plan%thermal_strain = .true.
      plan%stress = .true.
      plan%rotation = .true.
      plan%gap_displacement = .true.
      plan%gap_normal_traction = .true.
      if (solid_mechanics_viscoplasticity_enabled()) then
        plan%plastic_strain = .true.
        plan%plastic_strain_rate = .true.
      end if
      return
    end if

    if (.not.fields%has_selected_fields()) return

    do i = 1, size(fields%names)
      name = trim(fields%names(i))

      select case (name)
      case ('displacement')
        plan%displacement = .true.
      case ('strain')
        plan%strain = .true.
      case ('thermal_strain')
        plan%thermal_strain = .true.
      case ('stress')
        plan%stress = .true.
      case ('rotation')
        plan%rotation = .true.
      case ('contact_normal_displ')
        plan%gap_displacement = .true.
      case ('contact_pressure')
        plan%gap_normal_traction = .true.
      case ('plastic_strain')
        INSIST(solid_mechanics_viscoplasticity_enabled())
        plan%plastic_strain = .true.
      case ('plastic_strain_rate')
        INSIST(solid_mechanics_viscoplasticity_enabled())
        plan%plastic_strain_rate = .true.
      case default
        INSIST(.false.)
      end select
    end do

  end subroutine build_output_plan

end module solid_mechanics_viz_provider_type
