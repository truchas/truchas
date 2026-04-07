!!
!! SOLID_MECHANICS_VTKHDF_OUTPUT
!!
!! VTKHDF graphics output for solid mechanics fields.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module solid_mechanics_vtkhdf_output

  use,intrinsic :: iso_fortran_env, only: int8, r8 => real64
  use bitfield_type, only: bitfield
  use solid_mechanics_driver, only: solid_mechanics_enabled, solid_mechanics_get_point_fields, &
      solid_mechanics_get_cell_fields, solid_mechanics_get_contact_set_ids, &
      solid_mechanics_point_field_request, solid_mechanics_point_field_values, &
      solid_mechanics_cell_field_request, solid_mechanics_cell_field_values, &
      solid_mechanics_viscoplasticity_enabled
  use vtkhdf_mb_file_type, only: vtkhdf_mb_file, vtkhdf_block_handle, &
      vtkhdf_cell_data_handle, vtkhdf_point_data_handle, UG_FIXED_MESH
  implicit none
  private

  public :: solid_mechanics_vtkhdf_block_data
  public :: solid_mechanics_vtkhdf_init_output
  public :: solid_mechanics_vtkhdf_register_temporal_data
  public :: solid_mechanics_vtkhdf_write_output
  public :: solid_mechanics_vtkhdf_output_step

  !! The solid-mechanics VTKHDF data registered for one mesh block.
  type :: solid_mechanics_vtkhdf_block_data
    private
    type(vtkhdf_point_data_handle) :: displacement
    type(vtkhdf_cell_data_handle) :: strain
    type(vtkhdf_cell_data_handle) :: thermal_strain
    type(vtkhdf_cell_data_handle) :: stress
    type(vtkhdf_cell_data_handle) :: rotation
    type(vtkhdf_cell_data_handle) :: plastic_strain
    type(vtkhdf_cell_data_handle) :: plastic_strain_rate
  end type

  type :: solid_mechanics_vtkhdf_interface_block_data
    private
    integer :: setid = 0
    type(vtkhdf_block_handle) :: handle
    type(vtkhdf_point_data_handle) :: gap_displacement
    type(vtkhdf_point_data_handle) :: gap_normal_traction
    integer, allocatable :: nodes(:)
  end type

  type, public :: solid_mechanics_vtkhdf_output_data
    type(solid_mechanics_vtkhdf_interface_block_data), allocatable :: interface(:)
  end type

contains

  subroutine solid_mechanics_vtkhdf_init_output(file, data)

    use mesh_manager, only: unstr_mesh_ptr
    use string_utilities, only: i_to_c
    use unstr_mesh_type, only: unstr_mesh

    type(vtkhdf_mb_file), intent(inout) :: file
    type(solid_mechanics_vtkhdf_output_data), intent(out) :: data

    type(unstr_mesh), pointer :: mesh
    type(vtkhdf_block_handle) :: block
    type(bitfield) :: bitmask
    real(r8), allocatable :: x(:,:)
    integer, allocatable :: setids(:), xfnode(:), fnode(:), faces(:), nodes(:)
    integer(int8), allocatable :: types(:)
    integer :: n, stat
    character(:), allocatable :: errmsg

    if (.not.solid_mechanics_enabled()) return

    call solid_mechanics_get_contact_set_ids(setids)
    if (size(setids) == 0) return

    mesh => unstr_mesh_ptr('MAIN')
    INSIST(associated(mesh))

    allocate(data%interface(size(setids)))
    do n = 1, size(setids)
      call mesh%get_link_set_bitmask([setids(n)], bitmask, stat, errmsg)
      INSIST(stat == 0)
      call mesh%get_link_set_mesh(bitmask, x, xfnode, fnode, faces, nodes)
      types = get_vtk_face_types(xfnode)
      block = file%add_block('gap-interface-' // i_to_c(setids(n)), mode=UG_FIXED_MESH)
      call file%write_mesh(block, x, fnode, xfnode, types)
      data%interface(n)%setid = setids(n)
      data%interface(n)%handle = block
      data%interface(n)%gap_displacement = file%register_temporal_point_data(block, &
          'Gap Displacement', 0.0_r8)
      data%interface(n)%gap_normal_traction = file%register_temporal_point_data(block, &
          'Gap Normal Traction', 0.0_r8)
      call move_alloc(nodes, data%interface(n)%nodes)
    end do

  end subroutine solid_mechanics_vtkhdf_init_output

  subroutine solid_mechanics_vtkhdf_register_temporal_data(file, block, data)

    type(vtkhdf_mb_file), intent(inout) :: file
    type(vtkhdf_block_handle), intent(in) :: block
    type(solid_mechanics_vtkhdf_block_data), intent(inout) :: data

    real(r8) :: vector(3), tensor(6)

    if (.not.solid_mechanics_enabled()) return

    data%displacement = file%register_temporal_point_data(block, 'Displacement', vector)
    data%strain = file%register_temporal_cell_data(block, 'epsilon', tensor)
    data%thermal_strain = file%register_temporal_cell_data(block, 'epstherm', tensor)
    data%stress = file%register_temporal_cell_data(block, 'sigma', tensor)
    data%rotation = file%register_temporal_cell_data(block, 'Rotation', 0.0_r8)

    if (solid_mechanics_viscoplasticity_enabled()) then
      data%plastic_strain = file%register_temporal_cell_data(block, 'e_plastic', tensor)
      data%plastic_strain_rate = file%register_temporal_cell_data(block, 'epsdot', 0.0_r8)
    end if

  end subroutine solid_mechanics_vtkhdf_register_temporal_data

  subroutine solid_mechanics_vtkhdf_output_step(file, block, block_cells, block_nodes, &
      block_data)

    type(vtkhdf_mb_file), intent(inout) :: file
    type(vtkhdf_block_handle), intent(in) :: block
    integer, intent(in) :: block_cells(:), block_nodes(:)
    type(solid_mechanics_vtkhdf_block_data), intent(inout) :: block_data

    type(solid_mechanics_point_field_request) :: point_request
    type(solid_mechanics_point_field_values) :: point_values
    type(solid_mechanics_cell_field_request) :: cell_request
    type(solid_mechanics_cell_field_values) :: cell_values

    if (.not.solid_mechanics_enabled()) return

    point_request%displacement = .true.
    call solid_mechanics_get_point_fields(block_nodes, point_request, point_values)

    cell_request%total_strain = .true.
    cell_request%thermal_strain = .true.
    cell_request%stress = .true.
    cell_request%rotation = .true.
    if (solid_mechanics_viscoplasticity_enabled()) then
      cell_request%plastic_strain = .true.
      cell_request%plastic_strain_rate = .true.
    end if
    call solid_mechanics_get_cell_fields(block_cells, cell_request, cell_values)

    call to_vtk_symmetric_tensor_order(cell_values%total_strain)
    call to_vtk_symmetric_tensor_order(cell_values%thermal_strain)
    call to_vtk_symmetric_tensor_order(cell_values%stress)
    if (allocated(cell_values%plastic_strain)) &
        call to_vtk_symmetric_tensor_order(cell_values%plastic_strain)

    call file%write_point_data(block, block_data%displacement, &
        point_values%displacement)

    call file%write_cell_data(block, block_data%strain, cell_values%total_strain)
    call file%write_cell_data(block, block_data%thermal_strain, &
        cell_values%thermal_strain)
    call file%write_cell_data(block, block_data%stress, cell_values%stress)
    call file%write_cell_data(block, block_data%rotation, cell_values%rotation)

    if (solid_mechanics_viscoplasticity_enabled()) then
      call file%write_cell_data(block, block_data%plastic_strain, &
          cell_values%plastic_strain)
      call file%write_cell_data(block, block_data%plastic_strain_rate, &
          cell_values%plastic_strain_rate)
    end if

  end subroutine solid_mechanics_vtkhdf_output_step

  subroutine solid_mechanics_vtkhdf_write_output(file, data)

    type(vtkhdf_mb_file), intent(inout) :: file
    type(solid_mechanics_vtkhdf_output_data), intent(inout) :: data

    type(solid_mechanics_point_field_request) :: request
    type(solid_mechanics_point_field_values) :: values
    integer :: n

    if (.not.solid_mechanics_enabled()) return
    if (.not.allocated(data%interface)) return

    request%gap_displacement = .true.
    request%gap_normal_traction = .true.
    do n = 1, size(data%interface)
      call solid_mechanics_get_point_fields(data%interface(n)%nodes, request, values)
      call file%write_point_data(data%interface(n)%handle, data%interface(n)%gap_displacement, &
          values%gap_displacement)
      call file%write_point_data(data%interface(n)%handle, &
          data%interface(n)%gap_normal_traction, values%gap_normal_traction)
    end do

  end subroutine solid_mechanics_vtkhdf_write_output

  !! Reorder a Truchas tensor (xx, yy, zz, xy, xz, yz) to VTK tensor order.
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

end module solid_mechanics_vtkhdf_output
