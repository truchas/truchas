module ht_2d_vtkhdf_output

  use,intrinsic :: iso_fortran_env, only: int8, r8 => real64
  use unstr_2d_mesh_type
  use vtkhdf_mb_file_type, only: vtkhdf_mb_file, vtkhdf_block_handle, &
      vtkhdf_cell_data_handle, UG_FIXED_MESH
  implicit none
  private

  type, public :: ht_2d_vtkhdf_writer
    private
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
    type(unstr_2d_mesh), intent(in) :: mesh
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8), allocatable :: x(:,:)
    integer, allocatable :: xcnode(:), cnode(:)
    integer :: j, nnode
    integer(int8), allocatable :: types(:)

    call this%file%create('out.vtkhdf', comm, stat, errmsg)
    if (stat /= 0) return

    this%block = this%file%add_block('main', mode=UG_FIXED_MESH)
    allocate(x(3, mesh%nnode))
    x = 0.0_r8
    x(:2,:) = mesh%x(:, :mesh%nnode)
    xcnode = mesh%cstart(:mesh%ncell_onP+1)
    cnode = mesh%cnode(:mesh%cstart(mesh%ncell_onP+1)-1)
    allocate(types(mesh%ncell_onP))
    do j = 1, mesh%ncell_onP
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

    this%enthalpy = this%file%register_temporal_cell_data(this%block, 'H', 0.0_r8)
    this%temperature = this%file%register_temporal_cell_data(this%block, 'T', 0.0_r8)
    this%is_open = .true.

  end subroutine open

  subroutine write_solution(this, time, enthalpy, temperature)

    class(ht_2d_vtkhdf_writer), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(in) :: enthalpy(:), temperature(:)

    call this%file%start_time_step(time)
    call this%file%write_cell_data(this%block, this%enthalpy, enthalpy)
    call this%file%write_cell_data(this%block, this%temperature, temperature)
    call this%file%finalize_time_step()
    call this%file%flush()

  end subroutine write_solution

  subroutine close(this)
    class(ht_2d_vtkhdf_writer), intent(inout) :: this
    if (this%is_open) call this%file%close()
    this%is_open = .false.
  end subroutine close

end module ht_2d_vtkhdf_output
