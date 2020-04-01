!!
!! LIN_COMB_FUNC_CLASS
!!
!! A commonly encountered computational pattern is the evaluation of a linear
!! combination of functions of the form $f(x,w) = \sum_i w_i f_i(x)$. This
!! module defines the abstract base class LIN_COMB_FUNC that defines a common
!! interface for various concrete implementations.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! December 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module lin_comb_func_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: lin_comb_func
  contains
    generic :: compute_value => compute_value_1, compute_value_2
    generic :: compute_deriv => compute_deriv_1, compute_deriv_2
    procedure(compute_value_1), deferred :: compute_value_1
    procedure(compute_value_2), deferred :: compute_value_2
    procedure(compute_deriv_1), deferred :: compute_deriv_1
    procedure(compute_deriv_2), deferred :: compute_deriv_2
  end type

  abstract interface
    subroutine compute_value_1(this, w, state, value)
      import lin_comb_func, r8
      class(lin_comb_func), intent(in) :: this
      real(r8), intent(in) :: w(:), state(:)
      real(r8), intent(out) :: value
    end subroutine
    subroutine compute_value_2(this, w, state, value)
      import lin_comb_func, r8
      class(lin_comb_func), intent(in) :: this
      real(r8), intent(in) :: w(:,:), state(:,:)
      real(r8), intent(out) :: value(:)
    end subroutine
  end interface

  abstract interface
    subroutine compute_deriv_1(this, w, state, n, deriv)
      import lin_comb_func, r8
      class(lin_comb_func), intent(in) :: this
      real(r8), intent(in) :: w(:), state(:)
      integer, intent(in) :: n
      real(r8), intent(out) :: deriv
    end subroutine
    subroutine compute_deriv_2(this, w, state, n, deriv)
      import lin_comb_func, r8
      class(lin_comb_func), intent(in) :: this
      real(r8), intent(in) :: w(:,:), state(:,:)
      integer, intent(in) :: n
      real(r8), intent(out) :: deriv(:)
    end subroutine
  end interface

end module lin_comb_func_class
