!!
!! TRUCHAS_VTKHDF_OUTPUT
!!
!! High-level VTKHDF graphics output.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module truchas_vtkhdf_output

  use,intrinsic :: iso_fortran_env, only: int8, real64
  use diffusion_solver_vtkhdf_output, only: diffusion_solver_vtkhdf_block_data
  use em_heat_driver, only: em_heat_vtkhdf_block_data
  use flow_vtkhdf_output, only: flow_vtkhdf_block_data
  use parallel_communication, only: comm, is_IOP, this_PE
  use solid_mechanics_vtkhdf_output, only: solid_mechanics_vtkhdf_block_data, &
      solid_mechanics_vtkhdf_output_data
  use truchas_logging_services, only: TLS_panic
  use unstr_mesh_type, only: unstr_mesh
  use ustruc_driver, only: ustruc_vtkhdf_block_data
  use vtkhdf_mb_file_type, only: vtkhdf_mb_file, vtkhdf_block_handle, &
      vtkhdf_cell_data_handle, vtkhdf_point_data_handle, UG_FIXED_MESH, &
      UG_MOVING_MESH
  implicit none
  private

  public :: TVO_write_mesh, TVO_register_temporal_data, TVO_write_timestep, TVO_close

  type :: vtkhdf_block_data
    type(vtkhdf_block_handle) :: handle
    type(vtkhdf_cell_data_handle) :: process_ids
    type(vtkhdf_point_data_handle) :: global_point_ids
    type(vtkhdf_cell_data_handle), allocatable :: vof(:)
    type(diffusion_solver_vtkhdf_block_data) :: diffusion
    type(em_heat_vtkhdf_block_data) :: em_heat
    type(flow_vtkhdf_block_data) :: flow
    type(solid_mechanics_vtkhdf_block_data) :: solid_mechanics
    type(ustruc_vtkhdf_block_data) :: ustruc
    integer, allocatable :: cells(:), nodes(:)
    real(real64), allocatable :: x0(:,:)
    logical :: moving = .false.
  end type vtkhdf_block_data

  type(vtkhdf_mb_file), save :: file
  type(vtkhdf_block_data), allocatable, save :: block_data(:)
  type(solid_mechanics_vtkhdf_output_data), save :: solid_mechanics_output
  logical, save :: file_is_open = .false.
  logical, save :: temporal_data_registered = .false.
  real(real64), allocatable, save :: last_moving_translation(:)

contains

  subroutine TVO_write_mesh(mesh)

    use output_control, only: part, part_path
    use truchas_env, only: output_file_name, overwrite_output

    type(unstr_mesh), intent(in) :: mesh

    type(vtkhdf_block_handle) :: block
    character(:), allocatable :: errmsg
    real(real64), allocatable :: x(:,:)
    integer, allocatable :: xcnode(:), cnode(:), block_cells(:), block_nodes(:)
    integer(int8), allocatable :: types(:)
    integer(kind(mesh%cell_set_mask)) :: bitmask
    integer :: n, stat
    logical :: exists, moving

    if (is_IOP .and. .not. overwrite_output) then
      inquire(file=output_file_name('vtkhdf'), exist=exists)
      if (exists) then
        call TLS_panic("must specify `-f` flag to overwrite `" // &
            output_file_name('vtkhdf')//"`")
      end if
    endif

    call file%create(output_file_name('vtkhdf'), comm, stat, errmsg)
    if (stat /= 0) then
      if (allocated(errmsg)) then
        call TLS_panic('failed to create VTKHDF graphics file: ' // errmsg)
      else
        call TLS_panic('failed to create VTKHDF graphics file')
      end if
    end if

    allocate(block_data(size(mesh%cell_set_name)))
    do n = 1, size(mesh%cell_set_name)
      associate (name => mesh%cell_set_name(n)%s)
        moving = .false.
        if (allocated(part_path) .and. allocated(part)) &
            moving = any(part == mesh%cell_set_id(n))
        bitmask = 0
        bitmask = ibset(bitmask, pos=n)
        call mesh%get_mesh_block(bitmask, x, xcnode, cnode, &
            block_cells, block_nodes)
        types = get_vtk_cell_types(xcnode)
        if (moving) then
          block = file%add_block(name, mode=UG_MOVING_MESH)
          call file%write_mesh_topology(block, cnode, xcnode, types)
          call move_alloc(x, block_data(n)%x0)
        else
          block = file%add_block(name, mode=UG_FIXED_MESH)
          call file%write_mesh(block, x, cnode, xcnode, types)
          call file%write_cell_data(block, 'ProcessIds', &
              spread(this_PE, dim=1, ncopies=size(block_cells)))
          call file%write_point_data(block, 'vtkGlobalPointIds', &
              mesh%node_imap%global_index(block_nodes))
        end if
        block_data(n)%handle = block
        block_data(n)%moving = moving
        call move_alloc(block_cells, block_data(n)%cells)
        call move_alloc(block_nodes, block_data(n)%nodes)
      end associate
    end do

    file_is_open = .true.
    temporal_data_registered = .false.
    if (allocated(last_moving_translation)) deallocate(last_moving_translation)

  end subroutine TVO_write_mesh

  subroutine TVO_register_temporal_data

    use,intrinsic :: iso_fortran_env, only: int8, r8 => real64
    use diffusion_solver_vtkhdf_output, only: &
        diffusion_solver_vtkhdf_register_temporal_data
    use em_heat_driver, only: em_heat_vtkhdf_register_temporal_data
    use flow_vtkhdf_output, only: flow_vtkhdf_register_temporal_data
    use legacy_matl_api, only: nmat
    use material_model_driver, only: matl_model
    use solid_mechanics_vtkhdf_output, only: &
        solid_mechanics_vtkhdf_init_output, solid_mechanics_vtkhdf_register_temporal_data
    use ustruc_driver, only: ustruc_vtkhdf_register_temporal_data

    integer :: m, n

    if (.not.file_is_open .or. temporal_data_registered) return

    do n = 1, size(block_data)
      associate (b => block_data(n), block => block_data(n)%handle)
        if (b%moving) then
          b%process_ids = file%register_temporal_cell_data(block, 'ProcessIds', 0)
          b%global_point_ids = file%register_temporal_point_data(block, &
              'vtkGlobalPointIds', 0)
        end if
        call diffusion_solver_vtkhdf_register_temporal_data(file, block, b%cells, &
            b%moving, b%diffusion)
        if (nmat > 1) then
          allocate(b%vof(nmat))
          do m = 1, nmat
            b%vof(m) = file%register_temporal_cell_data(block, &
                'VOF_' // matl_model%phase_name(m), 0.0_r8)
          end do
        end if
        call em_heat_vtkhdf_register_temporal_data(file, block, b%em_heat)
        call flow_vtkhdf_register_temporal_data(file, block, b%flow)
        call solid_mechanics_vtkhdf_register_temporal_data(file, block, b%solid_mechanics)
        call ustruc_vtkhdf_register_temporal_data(file, block, b%ustruc)
      end associate
    end do
    call solid_mechanics_vtkhdf_init_output(file, solid_mechanics_output)
    temporal_data_registered = .true.

  end subroutine TVO_register_temporal_data

  subroutine TVO_write_timestep

    use time_step_module, only: t
    use zone_module, only: zone
    use diffusion_solver_vtkhdf_output, only: diffusion_solver_vtkhdf_output_step
    use em_heat_driver, only: em_heat_vtkhdf_output
    use flow_vtkhdf_output, only: flow_vtkhdf_output_step
    use legacy_matl_api, only: nmat, gather_vof
    use mesh_manager, only: unstr_mesh_ptr
    use output_control, only: part_path
    use solid_mechanics_driver, only: solid_mechanics_enabled
    use solid_mechanics_vtkhdf_output, only: solid_mechanics_vtkhdf_output_step, &
        solid_mechanics_vtkhdf_write_output
    use ustruc_driver, only: ustruc_vtkhdf_output

    integer :: m, n
    real(real64) :: r(3)
    real(real64), allocatable :: vof(:), x(:,:)
    type(unstr_mesh), pointer :: mesh
    logical :: write_moving_geometry

    if (.not.file_is_open) return
    INSIST(temporal_data_registered)

    write_moving_geometry = .false.
    if (any(block_data%moving)) then
      INSIST(allocated(part_path))
      call part_path%get_position(t, r)
      if (allocated(last_moving_translation)) then
        write_moving_geometry = any(r /= last_moving_translation)
      else
        write_moving_geometry = .true.
      end if
      mesh => unstr_mesh_ptr('MAIN')
      INSIST(associated(mesh))
    end if

    call file%start_time_step(t)
    do n = 1, size(block_data)
      associate (b => block_data(n))
        if (b%moving) then
          if (write_moving_geometry) then
            x = b%x0
            x(1,:) = x(1,:) + r(1)
            x(2,:) = x(2,:) + r(2)
            x(3,:) = x(3,:) + r(3)
            call file%write_mesh_geometry(b%handle, x)
          end if
          call file%write_cell_data(b%handle, b%process_ids, &
              spread(this_PE, dim=1, ncopies=size(b%cells)))
          call file%write_point_data(b%handle, b%global_point_ids, &
              mesh%node_imap%global_index(b%nodes))
        end if
        call diffusion_solver_vtkhdf_output_step(file, b%handle, b%cells, &
            b%moving, b%diffusion)
      end associate
    end do
    if (any(block_data%moving) .and. write_moving_geometry) then
      if (.not.allocated(last_moving_translation)) allocate(last_moving_translation(3))
      last_moving_translation = r
    end if
    if (nmat > 1) then
      allocate(vof(size(zone)))
      do m = 1, nmat
        call gather_vof(m, vof)
        do n = 1, size(block_data)
          associate (b => block_data(n))
            call file%write_cell_data(b%handle, b%vof(m), vof(b%cells))
          end associate
        end do
      end do
    end if
    do n = 1, size(block_data)
      associate (b => block_data(n))
        call em_heat_vtkhdf_output(file, b%handle, b%cells, b%em_heat)
      end associate
    end do
    do n = 1, size(block_data)
      associate (b => block_data(n))
        call flow_vtkhdf_output_step(file, b%handle, b%cells, b%flow)
      end associate
    end do
    if (solid_mechanics_enabled()) then
      do n = 1, size(block_data)
        associate (b => block_data(n))
          call solid_mechanics_vtkhdf_output_step(file, b%handle, b%cells, b%nodes, &
              b%solid_mechanics)
        end associate
      end do
      call solid_mechanics_vtkhdf_write_output(file, solid_mechanics_output)
    end if
    do n = 1, size(block_data)
      associate (b => block_data(n))
        call ustruc_vtkhdf_output(file, b%handle, b%cells, b%ustruc)
      end associate
    end do
    call file%finalize_time_step()
    call file%flush()

  end subroutine TVO_write_timestep

  subroutine TVO_close
    if (.not.file_is_open) return
    call file%close()
    if (allocated(block_data)) deallocate(block_data)
    if (allocated(solid_mechanics_output%interface)) deallocate(solid_mechanics_output%interface)
    file_is_open = .false.
    temporal_data_registered = .false.
    if (allocated(last_moving_translation)) deallocate(last_moving_translation)
  end subroutine

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

end module truchas_vtkhdf_output
