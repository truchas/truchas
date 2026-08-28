!!
!! MPOLY_VECTOR_FUNC_TYPE
!!
!! A concrete implementation of the abstract base class VECTOR_FUNC.  This
!! implementation defines a vector-valued multivariable polynomial with
!! user-specified coefficients and integral exponents.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! August 2026
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module mpoly_vector_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use vector_func_class
  implicit none
  private

  type, extends(vector_func), public :: mpoly_vector_func
    real(r8), allocatable :: x0(:)       ! reference point
    integer, allocatable :: expon(:,:)   ! exponents, coordinate x term
    real(r8), allocatable :: coef(:,:)   ! coefficients, component x term
  contains
    procedure :: eval
    procedure :: eval_comp
  end type mpoly_vector_func

  interface mpoly_vector_func
    procedure mpoly_vector_func_value
  end interface

contains

  function mpoly_vector_func_value(c, e, x0) result(f)
    real(r8), intent(in) :: c(:,:), x0(:)
    integer, intent(in) :: e(:,:)
    type(mpoly_vector_func) :: f

    INSIST(size(c,dim=2) > 0)
    INSIST(size(c,dim=2) == size(e,dim=2))
    INSIST(size(e,dim=1) == size(x0))

    f%dim = size(c,dim=1)
    f%x0 = x0
    f%expon = e
    f%coef = c
  end function mpoly_vector_func_value

  function eval(this, x) result(fx)
    class(mpoly_vector_func), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: fx(this%dim)

    integer :: i, j
    real(r8) :: t

    ASSERT(size(x) >= size(this%x0))

    fx = 0.0_r8
    do j = 1, size(this%coef,dim=2)
      t = 1.0_r8
      do i = 1, size(this%x0)
        if (this%expon(i,j) /= 0) t = t * (x(i)-this%x0(i))**this%expon(i,j)
      end do
      fx = fx + this%coef(:,j)*t
    end do
  end function eval

  function eval_comp(this, i, x) result(fx)
    class(mpoly_vector_func), intent(in) :: this
    integer, intent(in) :: i
    real(r8), intent(in) :: x(:)
    real(r8) :: fx

    integer :: j, k
    real(r8) :: t

    ASSERT(size(x) >= size(this%x0))

    fx = 0.0_r8
    do j = 1, size(this%coef,dim=2)
      t = this%coef(i,j)
      do k = 1, size(this%x0)
        if (this%expon(k,j) /= 0) t = t * (x(k)-this%x0(k))**this%expon(k,j)
      end do
      fx = fx + t
    end do
  end function eval_comp

end module mpoly_vector_func_type
