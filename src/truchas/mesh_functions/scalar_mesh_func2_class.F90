!!
!! SCALAR_MESH_FUNC2_CLASS
!!
!! Zach Jibben <zjibben@lanl.gov>
!! January 2024
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module scalar_mesh_func2_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: scalar_mesh_func2
    real(r8), allocatable :: value(:)
  contains
    procedure(compute), deferred :: compute
  end type scalar_mesh_func2

  abstract interface
    subroutine compute(this, t, v)
      import r8, scalar_mesh_func2
      class(scalar_mesh_func2), intent(inout) :: this
      real(r8), intent(in) :: t, v(:)
    end subroutine
  end interface

end module scalar_mesh_func2_class
