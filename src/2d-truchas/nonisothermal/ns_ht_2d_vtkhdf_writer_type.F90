!!
!! NS_HT_2D_VTKHDF_WRITER_TYPE
!!
!! This module defines NS_HT_2D_VTKHDF_WRITER, which writes the fixed mesh
!! and coupled two-dimensional Navier--Stokes/thermal solution to VTKHDF.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

module ns_ht_2d_vtkhdf_writer_type

  use,intrinsic :: iso_fortran_env, only: int8, r8 => real64
  use simulation_environment_type
  use unstr_2d_mesh_type
  use vtkhdf_mb_file_type, only: vtkhdf_mb_file, vtkhdf_block_handle, vtkhdf_cell_data_handle, UG_FIXED_MESH
  implicit none
  private

  type, public :: ns_ht_2d_vtkhdf_writer
    private
    type(unstr_2d_mesh), pointer :: mesh => null()
    type(vtkhdf_mb_file) :: file
    type(vtkhdf_block_handle) :: block
    type(vtkhdf_cell_data_handle) :: pressure, velocity, enthalpy, temperature
    logical :: is_open = .false.
  contains
    procedure :: open
    procedure :: write_solution
    procedure :: close
  end type

contains

  subroutine open(this, env, mesh, stat, errmsg)
    use vtkhdf_vtk_cell_types, only: VTK_TRIANGLE, VTK_QUAD

    class(ns_ht_2d_vtkhdf_writer), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(unstr_2d_mesh), target, intent(in) :: mesh
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8), allocatable :: x(:,:)
    integer, allocatable :: cnode(:), xcnode(:), global_cell_ids(:), global_node_ids(:)
    integer :: c, nnode
    integer(int8), allocatable :: cell_type(:), cell_ghost_type(:), node_ghost_type(:)
    real(r8) :: scalar_mold, vector_mold(3)

    call this%file%create(trim(env%output_dir)//'/out.vtkhdf', env%comm%mpi_val, stat, errmsg)
    if (stat /= 0) return
    this%mesh => mesh
    this%block = this%file%add_block('main', mode=UG_FIXED_MESH)
    allocate(x(3,mesh%nnode), cell_type(mesh%ncell))
    x = 0.0_r8
    x(:2,:) = mesh%x(:,:mesh%nnode)
    xcnode = mesh%cstart
    cnode = mesh%cnode
    do c = 1, mesh%ncell
      nnode = xcnode(c+1) - xcnode(c)
      select case (nnode)
      case (3); cell_type(c) = VTK_TRIANGLE
      case (4); cell_type(c) = VTK_QUAD
      case default
        stat = 1
        errmsg = 'unsupported 2D cell topology'
        call this%file%close()
        return
      end select
    end do
    call this%file%write_mesh(this%block, x, cnode, xcnode, cell_type)
    global_cell_ids = [(mesh%cell_imap%global_index(c), c=1,mesh%ncell)]
    call this%file%write_cell_data(this%block, 'GlobalCellIds', global_cell_ids, attribute='GlobalIds')
    call this%file%write_cell_data(this%block, 'ExternalCellIds', mesh%xcell, attribute='PedigreeIds')
    allocate(cell_ghost_type(mesh%ncell), source=0_int8)
    cell_ghost_type(mesh%ncell_onP+1:) = 1_int8
    call this%file%write_cell_data(this%block, 'vtkGhostType', cell_ghost_type)
    global_node_ids = [(mesh%node_imap%global_index(c), c=1,mesh%nnode)]
    call this%file%write_point_data(this%block, 'GlobalNodeIds', global_node_ids, attribute='GlobalIds')
    call this%file%write_point_data(this%block, 'ExternalNodeIds', mesh%xnode, attribute='PedigreeIds')
    allocate(node_ghost_type(mesh%nnode), source=0_int8)
    node_ghost_type(mesh%nnode_onP+1:) = 1_int8
    call this%file%write_point_data(this%block, 'vtkGhostType', node_ghost_type)
    this%pressure = this%file%register_temporal_cell_data(this%block, 'pressure', scalar_mold)
    this%velocity = this%file%register_temporal_cell_data(this%block, 'velocity', vector_mold)
    this%enthalpy = this%file%register_temporal_cell_data(this%block, 'H', scalar_mold)
    this%temperature = this%file%register_temporal_cell_data(this%block, 'T', scalar_mold)
    this%is_open = .true.
  end subroutine

  subroutine write_solution(this, time, pressure, velocity, enthalpy, temperature)
    class(ns_ht_2d_vtkhdf_writer), intent(inout) :: this
    real(r8), intent(in) :: time, pressure(:), velocity(:,:), enthalpy(:), temperature(:)
    real(r8), allocatable :: p(:), v(:,:), H(:), T(:)

    allocate(p(this%mesh%ncell), v(3,this%mesh%ncell), H(this%mesh%ncell), T(this%mesh%ncell))
    p = pressure(:this%mesh%ncell)
    v = 0.0_r8
    v(:2,:) = velocity(:2,:this%mesh%ncell)
    H(:this%mesh%ncell_onP) = enthalpy
    T(:this%mesh%ncell_onP) = temperature
    call this%mesh%cell_imap%gather_offp(H)
    call this%mesh%cell_imap%gather_offp(T)
    call this%file%start_time_step(time)
    call this%file%write_cell_data(this%block, this%pressure, p)
    call this%file%write_cell_data(this%block, this%velocity, v)
    call this%file%write_cell_data(this%block, this%enthalpy, H)
    call this%file%write_cell_data(this%block, this%temperature, T)
    call this%file%finalize_time_step()
    call this%file%flush()
  end subroutine

  subroutine close(this)
    class(ns_ht_2d_vtkhdf_writer), intent(inout) :: this
    if (this%is_open) call this%file%close()
    this%is_open = .false.
    nullify(this%mesh)
  end subroutine

end module ns_ht_2d_vtkhdf_writer_type
