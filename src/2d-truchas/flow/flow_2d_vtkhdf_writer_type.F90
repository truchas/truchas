!!
!! FLOW_2D_VTKHDF_WRITER_TYPE
!!
!! This module writes the mesh-associated state of a two-dimensional flow
!! simulation to a VTKHDF unstructured-grid file.  The mesh and its identifier data
!! are static; cell pressure and velocity are temporal datasets.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

module flow_2d_vtkhdf_writer_type

  use,intrinsic :: iso_fortran_env, only: int8, r8 => real64
  use simulation_environment_type
  use unstr_2d_mesh_type
  use vtkhdf_ug_file_type, only: vtkhdf_ug_file, vtkhdf_cell_data_handle, UG_FIXED_MESH
  implicit none
  private

  type, public :: flow_2d_vtkhdf_writer
    private
    type(unstr_2d_mesh), pointer :: mesh => null()
    type(vtkhdf_ug_file) :: file
    type(vtkhdf_cell_data_handle) :: pressure, velocity
    logical :: is_open = .false.
  contains
    procedure :: open
    procedure :: write_solution
    procedure :: close
  end type

contains

  subroutine open(this, env, mesh, stat, errmsg)
    use vtkhdf_vtk_cell_types, only: VTK_TRIANGLE, VTK_QUAD

    class(flow_2d_vtkhdf_writer), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(unstr_2d_mesh), target, intent(in) :: mesh
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8), allocatable :: x(:,:)
    integer, allocatable :: xcnode(:), cnode(:), global_cell_ids(:), global_node_ids(:)
    integer :: c, nnode
    integer(int8), allocatable :: types(:), cell_ghost_type(:), node_ghost_type(:)
    real(r8) :: vector_mold(3), scalar_mold

    call this%file%create(trim(env%output_dir)//'/out.vtkhdf', &
        env%comm%mpi_val, stat, errmsg, mode=UG_FIXED_MESH)
    if (stat /= 0) return
    this%mesh => mesh
    allocate(x(3, mesh%nnode))
    x = 0.0_r8
    x(:2,:) = mesh%x(:, :mesh%nnode)
    xcnode = mesh%cstart
    cnode = mesh%cnode
    allocate(types(mesh%ncell))
    do c = 1, mesh%ncell
      nnode = xcnode(c+1) - xcnode(c)
      select case (nnode)
      case (3)
        types(c) = VTK_TRIANGLE
      case (4)
        types(c) = VTK_QUAD
      case default
        stat = 1
        errmsg = 'unsupported 2D cell topology'
        call this%file%close()
        return
      end select
    end do
    call this%file%write_mesh(x, cnode, xcnode, types)

    global_cell_ids = [(mesh%cell_imap%global_index(c), c=1,mesh%ncell)]
    call this%file%write_cell_data('GlobalCellIds', global_cell_ids, attribute='GlobalIds')
    call this%file%write_cell_data('ExternalCellIds', mesh%xcell, attribute='PedigreeIds')
    allocate(cell_ghost_type(mesh%ncell), source=0_int8)
    cell_ghost_type(mesh%ncell_onP+1:) = 1_int8
    call this%file%write_cell_data('vtkGhostType', cell_ghost_type)
    global_node_ids = [(mesh%node_imap%global_index(c), c=1,mesh%nnode)]
    call this%file%write_point_data('GlobalNodeIds', global_node_ids, attribute='GlobalIds')
    call this%file%write_point_data('ExternalNodeIds', mesh%xnode, attribute='PedigreeIds')
    allocate(node_ghost_type(mesh%nnode), source=0_int8)
    node_ghost_type(mesh%nnode_onP+1:) = 1_int8
    call this%file%write_point_data('vtkGhostType', node_ghost_type)

    this%pressure = this%file%register_temporal_cell_data('pressure', scalar_mold)
    this%velocity = this%file%register_temporal_cell_data('velocity', vector_mold)
    this%is_open = .true.
  end subroutine


  subroutine write_solution(this, time, pressure, velocity)
    class(flow_2d_vtkhdf_writer), intent(inout) :: this
    real(r8), intent(in) :: time, pressure(:), velocity(:,:)

    real(r8), allocatable :: p(:), v(:,:)

    allocate(p(this%mesh%ncell), v(3,this%mesh%ncell))
    p = pressure(:this%mesh%ncell)
    v = 0.0_r8
    v(:2,:) = velocity(:2,:this%mesh%ncell)
    call this%file%start_time_step(time)
    call this%file%write_cell_data(this%pressure, p)
    call this%file%write_cell_data(this%velocity, v)
    call this%file%finalize_time_step()
    call this%file%flush()
  end subroutine


  subroutine close(this)
    class(flow_2d_vtkhdf_writer), intent(inout) :: this

    if (this%is_open) call this%file%close()
    this%is_open = .false.
    nullify(this%mesh)
  end subroutine

end module flow_2d_vtkhdf_writer_type
