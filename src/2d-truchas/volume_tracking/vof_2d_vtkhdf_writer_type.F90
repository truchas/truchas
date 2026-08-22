!!
!! VOF_2D_VTKHDF_WRITER_TYPE
!!
!! This module defines a VTKHDF writer for the two-dimensional volume
!! tracking tests. It writes the complete local unstructured mesh, including
!! ghost cells and adjacent nodes, together with mesh identifiers and one
!! temporal scalar volume-fraction dataset for each material slot.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

module vof_2d_vtkhdf_writer_type

  use,intrinsic :: iso_fortran_env, only: int8, r8 => real64
  use simulation_environment_type
  use unstr_2d_mesh_type
  use vtkhdf_mb_file_type, only: vtkhdf_mb_file, vtkhdf_block_handle, &
      vtkhdf_cell_data_handle, UG_FIXED_MESH
  implicit none
  private

  type, public :: vof_2d_vtkhdf_writer
    private
    type(unstr_2d_mesh), pointer :: mesh => null()
    type(vtkhdf_mb_file) :: file
    type(vtkhdf_block_handle) :: block
    type(vtkhdf_cell_data_handle), allocatable :: volume_fraction(:)
    logical :: is_open = .false.
  contains
    procedure :: open
    procedure :: write_solution
    procedure :: close
  end type vof_2d_vtkhdf_writer

contains

  subroutine open(this, env, mesh, nmat, stat, errmsg)

    use vtkhdf_vtk_cell_types, only: VTK_TRIANGLE, VTK_QUAD

    class(vof_2d_vtkhdf_writer), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(unstr_2d_mesh), target, intent(in) :: mesh
    integer, intent(in) :: nmat
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    character(32) :: name
    real(r8), allocatable :: x(:,:)
    integer, allocatable :: xcnode(:), cnode(:), global_cell_ids(:), global_node_ids(:)
    integer :: j, nnode
    integer(int8), allocatable :: types(:), cell_ghost_type(:), node_ghost_type(:)

    stat = 0
    errmsg = ''
    if (nmat < 1) then
      stat = 1
      errmsg = 'volume-tracking VTKHDF output requires at least one material'
      return
    end if

    call this%file%create(trim(env%output_dir)//'/out.vtkhdf', &
        env%comm%mpi_val, stat, errmsg)
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

    allocate(this%volume_fraction(nmat))
    do j = 1, nmat
      write(name, '(a,i0)') 'volume-fraction-', j
      this%volume_fraction(j) = &
          this%file%register_temporal_cell_data(this%block, trim(name), 0.0_r8)
    end do
    this%is_open = .true.

  end subroutine open


  subroutine write_solution(this, time, volume_fraction)

    class(vof_2d_vtkhdf_writer), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(in) :: volume_fraction(:,:)

    real(r8), allocatable :: v(:)
    integer :: j

    if (size(volume_fraction, 1) < size(this%volume_fraction) .or. &
        size(volume_fraction, 2) < this%mesh%ncell_onP) return

    allocate(v(this%mesh%ncell))
    call this%file%start_time_step(time)
    do j = 1, size(this%volume_fraction)
      v(:this%mesh%ncell_onP) = volume_fraction(j,:this%mesh%ncell_onP)
      call this%mesh%cell_imap%gather_offp(v)
      call this%file%write_cell_data(this%block, this%volume_fraction(j), v)
    end do
    call this%file%finalize_time_step()
    call this%file%flush()

  end subroutine write_solution


  subroutine close(this)

    class(vof_2d_vtkhdf_writer), intent(inout) :: this

    if (this%is_open) call this%file%close()
    this%is_open = .false.
    if (allocated(this%volume_fraction)) deallocate(this%volume_fraction)
    nullify(this%mesh)

  end subroutine close

end module vof_2d_vtkhdf_writer_type
