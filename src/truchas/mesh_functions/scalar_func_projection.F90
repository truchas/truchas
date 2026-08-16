!! SCALAR_FUNC_PROJECTION
!!
!! General low-order projections of scalar functions onto mesh entities.

#include "f90_assert.fpp"

module scalar_func_projection

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_base_mesh_class
  use scalar_func_class
  implicit none
  private

  public :: project_scalar_func_to_cell_centers

contains

  subroutine project_scalar_func_to_cell_centers(mesh, f, values)
    class(unstr_base_mesh), intent(inout) :: mesh
    class(scalar_func), intent(in) :: f
    real(r8), intent(out) :: values(:)

    real(r8) :: args(size(mesh%x, dim=1))
    integer :: j, dim

    dim = size(mesh%x, dim=1)
    ASSERT(size(values) == mesh%ncell_onP)
    call mesh%init_cell_centroid
    do j = 1, mesh%ncell_onP
      args = mesh%cell_centroid(:,j)
      values(j) = f%eval(args)
    end do
  end subroutine project_scalar_func_to_cell_centers

end module scalar_func_projection
