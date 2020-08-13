!!
!! SCALAR_MESH_FUNC_CLASS
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! September 2014
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module scalar_mesh_func_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: scalar_mesh_func
    real(r8), allocatable :: value(:)
  contains
    procedure(compute), deferred :: compute
  end type scalar_mesh_func

  abstract interface
    subroutine compute(this, t)
      import r8, scalar_mesh_func
      class(scalar_mesh_func), intent(inout) :: this
      real(r8), intent(in) :: t
    end subroutine
  end interface

end module scalar_mesh_func_class
