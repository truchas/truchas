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
  use vtkhdf_mb_file_type
  use vtkhdf_vtk_cell_types, only: VTK_TETRA
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

    type(vtkhdf_mb_file) :: viz_file
    type(vtkhdf_block_handle) :: blk
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
        blk = viz_file%add_block(name, mode=UG_STATIC)
        if (stat /= 0) return
        bitmask = ibset(0,pos=j)
        call solver%mesh%get_mesh_block(bitmask, x, xcnode, cnode, block_cells, block_nodes)
        types = spread(VTK_TETRA, dim=1, ncopies=size(xcnode)-1)
        call viz_file%write_mesh(blk, x, cnode, xcnode, types)

        ncell = size(block_cells) ! number of local cells in block

        scalar = q(block_cells)
        call viz_file%write_cell_data(blk, 'Q_EM', scalar)

        zvector = e(:,block_cells)
        call viz_file%write_cell_data(blk, 'E_re', zvector%re)
        INSIST(stat == 0)

#ifdef GNU_PR117774
        call viz_file%write_cell_data(blk, 'E_im', reshape([zvector%im],shape(zvector)))
#else
        call viz_file%write_cell_data(blk, 'E_im', zvector%im)
#endif

        call viz_file%write_cell_data(blk, '|E|', abs(zvector))

        zvector = h(:,block_cells)
        call viz_file%write_cell_data(blk, 'H_re', zvector%re)

#ifdef GNU_PR117774
        call viz_file%write_cell_data(blk, 'H_im', reshape([zvector%im],shape(zvector)))
#else
        call viz_file%write_cell_data(blk, 'H_im', zvector%im)
#endif
        INSIST(stat == 0)

        call viz_file%write_cell_data(blk, '|H|', abs(zvector))

        !! Output the mesh partition FIXME: output integer; name "ProcessIds" ???
        scalar = spread(real(this_PE,kind=r8), dim=1, ncopies=ncell)
        call viz_file%write_cell_data(blk, 'MPI rank', scalar)

        zscalar = dd(block_nodes)
        call viz_file%write_point_data(blk, 'div_D_re', zscalar%re)

#ifdef GNU_PR117774
        call viz_file%write_point_data(blk, 'div_D_im', [zscalar%im])
#else
        call viz_file%write_point_data(blk, 'div_D_im', zscalar%im)
#endif

        call viz_file%write_point_data(blk, '|div_D|', abs(zscalar))
        INSIST(stat == 0)

        !SPECULATIVE!
        !call viz_file%write_point_data(blk, 'vtkGlobalPointIds', solver%mesh%xnode(block_nodes))
      end associate
    end do

    call viz_file%close

  end subroutine fdme_vtk_graphics

end module fdme_vtk_graphics_proc
