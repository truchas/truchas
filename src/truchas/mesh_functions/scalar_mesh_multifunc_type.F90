!!
!! SCALAR_MESH_MULTIFUNC
!!
!! Zach Jibben <zjibben@lanl.gov>
!! January 2024
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module scalar_mesh_multifunc_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use scalar_mesh_func_class
  use scalar_mesh_func2_class
  implicit none
  private

  type, public :: scalar_mesh_multifunc
    real(r8), allocatable :: value(:)
    class(scalar_mesh_func), allocatable :: f1
    class(scalar_mesh_func2), allocatable :: f2
  contains
    procedure :: compute
  end type scalar_mesh_multifunc

contains

  subroutine compute(this, t, v)

    class(scalar_mesh_multifunc), intent(inout) :: this
    real(r8), intent(in) :: t, v(:)

    integer :: i

    this%value = 0

    if (allocated(this%f1)) then
      call this%f1%compute(t)
      this%value = this%value + this%f1%value
    end if

    if (allocated(this%f2)) then
      call this%f2%compute(t, v)
      this%value = this%value + this%f2%value
    end if

  end subroutine compute

end module scalar_mesh_multifunc_type
