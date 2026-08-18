!!
!! VECTOR_FUNC_PROJECTION
!!
!! General low-order projections of vector functions onto mesh entities.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause

#include "f90_assert.fpp"

module vector_func_projection

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_base_mesh_class
  use vector_func_class
  implicit none
  private

  public :: project_vector_func_to_cell_centers

contains

  subroutine project_vector_func_to_cell_centers(mesh, f, values)
    class(unstr_base_mesh), intent(inout) :: mesh
    class(vector_func), intent(in) :: f
    real(r8), intent(out) :: values(:,:)

    real(r8) :: args(size(mesh%x, dim=1))
    integer :: j

    ASSERT(size(values,1) == f%dim)
    ASSERT(size(values,2) == mesh%ncell_onP)
    call mesh%init_cell_centroid
    do j = 1, mesh%ncell_onP
      args = mesh%cell_centroid(:,j)
      values(:,j) = f%eval(args)
    end do
  end subroutine project_vector_func_to_cell_centers

end module vector_func_projection
