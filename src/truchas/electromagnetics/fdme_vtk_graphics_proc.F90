!!
!! FDME_VTK_GRAPHICS_PROC
!!
!! This module provides a procedure for writing a vtkhdf-format graphics file
!! of the solution to the frequency-domain Maxwell equations computed by an
!! FDME_SOLVER object and associated derived quantities.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module fdme_vtk_graphics_proc

  use,intrinsic :: iso_fortran_env, only: r8 => real64, int8
  use fdme_solver_type
  use vtkhdf_file_type
  use parallel_communication
  use string_utilities, only: i_to_c
  implicit none
  private

  public :: fdme_vtk_graphics

contains

  !TODO: Give user some control over what fields are output
  subroutine fdme_vtk_graphics(solver, filename, stat, errmsg)

    type(fdme_solver), intent(in) :: solver
    character(*), intent(in) :: filename
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(vtkhdf_file) :: viz_file
    real(r8), allocatable :: x(:,:)
    integer, allocatable :: xcnode(:), cnode(:), block_cells(:), block_nodes(:)
    integer(int8), allocatable :: types(:)
    real(r8), allocatable :: scalar(:)
    complex(r8), allocatable :: zvector(:,:), zscalar(:)
    integer, allocatable :: g_iscalar(:)
    real(r8) :: q(solver%mesh%ncell)
    complex(r8) :: e(3,solver%mesh%ncell), h(3,solver%mesh%ncell), dd(solver%mesh%nnode)
    integer :: j, n, bitmask, ncell, nnode

    call viz_file%create(filename, comm, stat, errmsg)
    if (stat /= 0) return

    call solver%get_heat_source(q)
    call solver%get_cell_efield(e)
    call solver%get_cell_hfield(h)
    call solver%get_div_dfield(dd)
    call solver%mesh%node_imap%gather_offp(dd)

    do j = 1, size(solver%mesh%cell_set_name)
      associate (name => solver%mesh%cell_set_name(j)%s)
        call viz_file%create_block(name, stat, errmsg)
        if (stat /= 0) return
        bitmask = ibset(0,pos=j)
        call solver%mesh%get_mesh_block(bitmask, x, xcnode, cnode, block_cells, block_nodes)
        types = spread(VTK_TETRA, dim=1, ncopies=size(xcnode)-1)
        call viz_file%write_block_mesh(name, x, cnode, xcnode, types, stat, errmsg)
        INSIST(stat == 0)

        ncell = size(block_cells) ! number of local cells in block

        scalar = q(block_cells)
        call viz_file%write_cell_dataset(name, 'Q_EM', scalar, stat, errmsg)
        INSIST(stat == 0)

        zvector = e(:,block_cells)
        call viz_file%write_cell_dataset(name, 'E_re', zvector%re, stat, errmsg)
        INSIST(stat == 0)

#ifdef GNU_PR117774
        call viz_file%write_cell_dataset(name, 'E_im', reshape([zvector%im],shape(zvector)), stat, errmsg)
#else
        call viz_file%write_cell_dataset(name, 'E_im', zvector%im, stat, errmsg)
#endif
        INSIST(stat == 0)

        call viz_file%write_cell_dataset(name, '|E|', abs(zvector), stat, errmsg)
        INSIST(stat == 0)

        zvector = h(:,block_cells)
        call viz_file%write_cell_dataset(name, 'H_re', zvector%re, stat, errmsg)
        INSIST(stat == 0)

#ifdef GNU_PR117774
        call viz_file%write_cell_dataset(name, 'H_im', reshape([zvector%im],shape(zvector)), stat, errmsg)
#else
        call viz_file%write_cell_dataset(name, 'H_im', zvector%im, stat, errmsg)
#endif
        INSIST(stat == 0)

        call viz_file%write_cell_dataset(name, '|H|', abs(zvector), stat, errmsg)
        INSIST(stat == 0)

        !! Output the mesh partition
        scalar = spread(real(this_PE,kind=r8), dim=1, ncopies=ncell)
        call viz_file%write_cell_dataset(name, 'MPI rank', scalar, stat, errmsg)
        INSIST(stat == 0)

        zscalar = dd(block_nodes)
        call viz_file%write_point_dataset(name, 'div_D_re', zscalar%re, stat, errmsg)
        INSIST(stat == 0)

#ifdef GNU_PR117774
        call viz_file%write_point_dataset(name, 'div_D_im', [zscalar%im], stat, errmsg)
#else
        call viz_file%write_point_dataset(name, 'div_D_im', zscalar%im, stat, errmsg)
#endif
        INSIST(stat == 0)

        call viz_file%write_point_dataset(name, '|div_D|', abs(zscalar), stat, errmsg)
        INSIST(stat == 0)
      end associate
    end do

    call viz_file%close

  end subroutine fdme_vtk_graphics

end module fdme_vtk_graphics_proc
