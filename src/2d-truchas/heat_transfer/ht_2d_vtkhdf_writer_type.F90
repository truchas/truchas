!!
!! HT_2D_VTKHDF_WRITER_TYPE
!!
!! This module defines the VTKHDF writer for the two-dimensional heat-transfer
!! simulation. It writes the fixed unstructured mesh, mesh-associated global
!! and pedigree identifiers, ghost-cell metadata, and temporal cell enthalpy
!! and temperature data to a VTKHDF multiblock file.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

module ht_2d_vtkhdf_writer_type

  use,intrinsic :: iso_fortran_env, only: int8, r8 => real64
  use unstr_2d_mesh_type
  use vtkhdf_mb_file_type, only: vtkhdf_mb_file, vtkhdf_block_handle, &
      vtkhdf_cell_data_handle, UG_FIXED_MESH
  implicit none
  private

  type, public :: ht_2d_vtkhdf_writer
    private
    type(unstr_2d_mesh), pointer :: mesh => null()
    type(vtkhdf_mb_file) :: file
    type(vtkhdf_block_handle) :: block
    type(vtkhdf_cell_data_handle) :: enthalpy
    type(vtkhdf_cell_data_handle) :: temperature
    logical :: is_open = .false.
  contains
    procedure :: open
    procedure :: write_solution
    procedure :: close
  end type ht_2d_vtkhdf_writer

contains

  subroutine open(this, mesh, stat, errmsg)

    use parallel_communication, only: comm
    use vtkhdf_vtk_cell_types, only: VTK_TRIANGLE, VTK_QUAD

    class(ht_2d_vtkhdf_writer), intent(out) :: this
    type(unstr_2d_mesh), target, intent(in) :: mesh
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8), allocatable :: x(:,:)
    integer, allocatable :: xcnode(:), cnode(:), global_cell_ids(:), global_node_ids(:)
    integer :: j, nnode
    integer(int8), allocatable :: types(:), cell_ghost_type(:), node_ghost_type(:)

    call this%file%create('out.vtkhdf', comm, stat, errmsg)
    if (stat /= 0) return

    this%mesh => mesh
    this%block = this%file%add_block('main', mode=UG_FIXED_MESH)
    allocate(x(3, mesh%nnode))
    x = 0.0_r8
    x(:2,:) = mesh%x(:, :mesh%nnode)
    xcnode = mesh%cstart
    cnode = mesh%cnode
    allocate(types(mesh%ncell))
    do j = 1, mesh%ncell
      nnode = xcnode(j+1) - xcnode(j)
      select case (nnode)
      case (3)
        types(j) = VTK_TRIANGLE
      case (4)
        types(j) = VTK_QUAD
      case default
        stat = 1
        errmsg = 'unsupported 2D cell topology'
        call this%file%close()
        return
      end select
    end do
    call this%file%write_mesh(this%block, x, cnode, xcnode, types)
    global_cell_ids = [(mesh%cell_imap%global_index(j), j=1,mesh%ncell)]
    call this%file%write_cell_data(this%block, 'GlobalCellIds', &
        global_cell_ids, attribute='GlobalIds')
    call this%file%write_cell_data(this%block, 'ExternalCellIds', mesh%xcell, &
        attribute='PedigreeIds')
    allocate(cell_ghost_type(mesh%ncell), source=0_int8)
    cell_ghost_type(mesh%ncell_onP+1:) = 1_int8
    call this%file%write_cell_data(this%block, 'vtkGhostType', cell_ghost_type)
    global_node_ids = [(mesh%node_imap%global_index(j), j=1,mesh%nnode)]
    call this%file%write_point_data(this%block, 'GlobalNodeIds', &
        global_node_ids, attribute='GlobalIds')
    call this%file%write_point_data(this%block, 'ExternalNodeIds', mesh%xnode, &
        attribute='PedigreeIds')
    allocate(node_ghost_type(mesh%nnode), source=0_int8)
    node_ghost_type(mesh%nnode_onP+1:) = 1_int8
    call this%file%write_point_data(this%block, 'vtkGhostType', node_ghost_type)

    this%enthalpy = this%file%register_temporal_cell_data(this%block, 'H', 0.0_r8)
    this%temperature = this%file%register_temporal_cell_data(this%block, 'T', 0.0_r8)
    this%is_open = .true.

  end subroutine open

  subroutine write_solution(this, time, enthalpy, temperature)

    class(ht_2d_vtkhdf_writer), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(in) :: enthalpy(:), temperature(:)
    real(r8), allocatable :: H(:), T(:)

    call this%file%start_time_step(time)
    allocate(H(this%mesh%ncell), T(this%mesh%ncell))
    H(:this%mesh%ncell_onP) = enthalpy
    T(:this%mesh%ncell_onP) = temperature
    call this%mesh%cell_imap%gather_offp(H)
    call this%mesh%cell_imap%gather_offp(T)
    call this%file%write_cell_data(this%block, this%enthalpy, H)
    call this%file%write_cell_data(this%block, this%temperature, T)
    call this%file%finalize_time_step()
    call this%file%flush()

  end subroutine write_solution

  subroutine close(this)
    class(ht_2d_vtkhdf_writer), intent(inout) :: this
    if (this%is_open) call this%file%close()
    this%is_open = .false.
    nullify(this%mesh)
  end subroutine close

end module ht_2d_vtkhdf_writer_type
