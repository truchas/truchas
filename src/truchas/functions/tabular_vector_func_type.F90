!!
!! TABULAR_VECTOR_FUNC_TYPE
!!
!! A concrete implementation of the abstract base class VECTOR_FUNC.  This
!! implementation defines a tabular function given by user-specified data
!! points with intervening linear interpolation.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! April 2014
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module tabular_vector_func_type

  use kinds, only: r8
  use vector_func_class
  implicit none
  private

  type, extends(vector_func), public :: tabular_vector_func
    private
    real(r8), allocatable :: x(:), y(:,:)
  contains
    procedure :: eval
    procedure :: eval_comp
  end type

  interface tabular_vector_func
    procedure tabular_vector_func_value
  end interface

contains

  function tabular_vector_func_value (x, y) result (f)
    real(r8), intent(in) :: x(:), y(:,:)
    type(tabular_vector_func) :: f
    ASSERT(size(x) > 1)
    ASSERT(size(y,dim=1) >= 1)
    ASSERT(size(x) == size(y,dim=2))
    f%dim = size(y,dim=1)
    f%x = x
    f%y = y
  end function tabular_vector_func_value

  function eval (this, x) result (fx)
    class(tabular_vector_func), intent(in) :: this
    real(r8), intent(in) :: x(:)  ! only x(1) is used
    real(r8) :: fx(this%dim)
    integer :: n, j, j1, j2
    n = size(this%x)
    if (x(1) <= this%x(1)) then
      fx = this%y(:,1)
    else if (x(1) >= this%x(n)) then
      fx = this%y(:,n)
    else
      !! Binary search to find the interval x(j1) < x <= x(j2), j2=j1+1.
      j1 = 1; j2 = n
      do while (j2 - j1 > 1)
        j = (j1 + j2) / 2
        if (x(1) > this%x(j)) then
          j1 = j
        else
          j2 = j
        end if
      end do
      !! Linearly interpolate over the interval
      fx = ((this%x(j2)-x(1))*this%y(:,j1) + (x(1)-this%x(j1))*this%y(:,j2))/(this%x(j2)-this%x(j1))
    end if
  end function eval

  function eval_comp (this, i, x) result (fx)
    class(tabular_vector_func), intent(in) :: this
    integer,  intent(in) :: i     ! desired function component
    real(r8), intent(in) :: x(:)  ! only x(1) is used
    real(r8) :: fx
    integer :: n, j, j1, j2
    n = size(this%x)
    if (x(1) <= this%x(1)) then
      fx = this%y(i,1)
    else if (x(1) >= this%x(n)) then
      fx = this%y(i,n)
    else
      !! Binary search to find the interval x(j1) < x <= x(j2), j2=j1+1.
      j1 = 1; j2 = n
      do while (j2 - j1 > 1)
        j = (j1 + j2) / 2
        if (x(1) > this%x(j)) then
          j1 = j
        else
          j2 = j
        end if
      end do
      !! Linearly interpolate over the interval
      fx = ((this%x(j2)-x(1))*this%y(i,j1) + (x(1)-this%x(j1))*this%y(i,j2))/(this%x(j2)-this%x(j1))
    end if
  end function eval_comp

end module tabular_vector_func_type
