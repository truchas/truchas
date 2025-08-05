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

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use fdme_solver_type
  use vtkhdf_file_type
  use parallel_communication
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
    real(r8), allocatable :: g_scalar(:), l_scalar(:)
    complex(r8), allocatable :: g_vector(:,:), l_vector(:,:), g_zscalar(:), l_zscalar(:)

    if (is_IOP) call viz_file%create(filename, stat, errmsg)
    call broadcast(stat)
    if (stat /= 0) then
      call broadcast(errmsg)
      return
    end if

    call export_mesh

    allocate(l_scalar(solver%mesh%ncell_onp))
    call solver%get_heat_source(l_scalar)
    allocate(g_scalar(merge(solver%mesh%cell_imap%global_size, 0, is_IOP)))
    call gather(l_scalar, g_scalar)
    if (is_IOP) call viz_file%write_cell_dataset('Q_EM', g_scalar, stat, errmsg)
    call broadcast(stat)
    INSIST(stat == 0)

    allocate(l_vector(3,solver%mesh%ncell_onp))
    allocate(g_vector(3,merge(solver%mesh%cell_imap%global_size, 0, is_IOP)))

    call solver%get_cell_efield(l_vector)
    call gather(l_vector, g_vector)
    if (is_IOP) call viz_file%write_cell_dataset('E_re', g_vector%re, stat, errmsg)
    call broadcast(stat)
    INSIST(stat == 0)

#ifdef GNU_PR117774
    if (is_IOP) call viz_file%write_cell_dataset('E_im', reshape([g_vector%im],shape(g_vector)), stat, errmsg)
#else
    if (is_IOP) call viz_file%write_cell_dataset('E_im', g_vector%im, stat, errmsg)
#endif
    call broadcast(stat)
    INSIST(stat == 0)

    if (is_IOP) call viz_file%write_cell_dataset('|E|', abs(g_vector), stat, errmsg)

    call solver%get_cell_hfield(l_vector)
    call gather(l_vector, g_vector)
    if (is_IOP) call viz_file%write_cell_dataset('H_re', g_vector%re, stat, errmsg)
    call broadcast(stat)
    INSIST(stat == 0)

#ifdef GNU_PR117774
    if (is_IOP) call viz_file%write_cell_dataset('H_im', reshape([g_vector%im],shape(g_vector)), stat, errmsg)
#else
    if (is_IOP) call viz_file%write_cell_dataset('H_im', g_vector%im, stat, errmsg)
#endif
    call broadcast(stat)
    INSIST(stat == 0)

    if (is_IOP) call viz_file%write_cell_dataset('|H|', abs(g_vector), stat, errmsg)

    !! Output the mesh partition
    call gather(spread(real(this_PE,kind=r8), dim=1, ncopies=solver%mesh%ncell_onP), g_scalar)
    if (is_IOP) call viz_file%write_cell_dataset('MPI rank', g_scalar, stat, errmsg)
    call broadcast(stat)
    INSIST(stat == 0)

    !! Divergence of the electric flux
    allocate(l_zscalar(solver%mesh%nnode_onp))
    allocate(g_zscalar(merge(solver%mesh%node_imap%global_size, 0, is_IOP)))
    call solver%get_div_dfield(l_zscalar)
    call gather(l_zscalar, g_zscalar)
    if (is_IOP) call viz_file%write_point_dataset('div_D_re', g_zscalar%re, stat, errmsg)
    call broadcast(stat)
#ifdef GNU_PR117774
    if (is_IOP) call viz_file%write_point_dataset('div_D_im', [g_zscalar%im], stat, errmsg)
#else
    if (is_IOP) call viz_file%write_point_dataset('div_D_im', g_zscalar%im, stat, errmsg)
#endif
    call broadcast(stat)
    if (is_IOP) call viz_file%write_point_dataset('|div_D|', abs(g_zscalar), stat, errmsg)
    call broadcast(stat)
    INSIST(stat == 0)

    if (is_IOP) call viz_file%close

  contains

    subroutine export_mesh

      use,intrinsic :: iso_fortran_env, only: int8

      integer, allocatable, target :: cnode(:,:)
      integer, allocatable :: xcnode(:)
      integer(int8), allocatable :: types(:)
      real(r8), allocatable :: x(:,:)
      integer, pointer :: connectivity(:)
      integer :: j, stat
      character(:), allocatable :: errmsg

      !! Collate the mesh data structure onto the IO process
      call solver%mesh%get_global_cnode_array(cnode)
      call solver%mesh%get_global_x_array(x)

      if (is_IOP) then
        xcnode = [(1+4*j, j=0, size(cnode,dim=2))]
        connectivity(1:size(cnode)) => cnode ! flattened view
        types = spread(VTK_TETRA, dim=1, ncopies=size(cnode,dim=2))
        call viz_file%write_mesh(x, connectivity, xcnode, types, stat, errmsg)
      end if
      call broadcast(stat)
      INSIST(stat == 0)

    end subroutine export_mesh

  end subroutine fdme_vtk_graphics

end module fdme_vtk_graphics_proc
