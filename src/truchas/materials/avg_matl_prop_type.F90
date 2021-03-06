!!
!! AVG_MATL_PROP_TYPE
!!
!! A freqently encountered pattern is the evaluation of a linear combination
!! function of the form $f(x, w) = \sum_{i=1}^m w_i f_i(x)$ for a fixed set
!  of functions $f_i$. This module defines the type AVG_MATL_PROP that
!! encapsulates this linear combination for the case of MATL_PROP class
!! functions.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! December 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module avg_matl_prop_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use matl_prop_class
  implicit none
  private

  type :: matl_prop_box
    class(matl_prop), allocatable :: prop
  end type

  type, public :: avg_matl_prop
    !private ! public only for matl_model:alloc_avg_matl_prop
    type(matl_prop_box), allocatable :: matl(:)
  contains
    generic   :: compute_value => compute_value_1, compute_value_2
    generic   :: compute_deriv => compute_deriv_1, compute_deriv_2
    procedure, private :: compute_value_1, compute_value_2
    procedure, private :: compute_deriv_1, compute_deriv_2
  end type

contains

  subroutine compute_value_1(this, w, state, value)
    class(avg_matl_prop), intent(in) :: this
    real(r8), intent(in) :: w(:), state(:)
    real(r8), intent(out) :: value
    integer :: n
    real(r8) :: tmp1, tmp2
    ASSERT(size(w) == size(this%matl))
    tmp2 = 0.0_r8
    do n = 1, size(w)
      if (w(n) > 0) then
        call this%matl(n)%prop%compute_value(state, tmp1)
        tmp2 = tmp2 + w(n)*tmp1
      end if
    end do
    value = tmp2
  end subroutine compute_value_1

  subroutine compute_value_2(this, w, state, value)
    class(avg_matl_prop), intent(in) :: this
    real(r8), intent(in) :: w(:,:), state(:,:)
    real(r8), intent(out) :: value(:)
    integer :: n, j
    real(r8) :: tmp1
    ASSERT(size(w,2) == size(this%matl))
    ASSERT(size(w,1) == size(state,2))
    ASSERT(size(value) == size(state,2))
    value = 0.0_r8
    do n = 1, size(w,2)
      do j = 1, size(w,1)
        if (w(j,n) > 0) then
          call this%matl(n)%prop%compute_value(state(:,j), tmp1)
          value(j) = value(j) + w(j,n)*tmp1
        end if
      end do
    end do
  end subroutine compute_value_2

  subroutine compute_deriv_1(this, w, state, n, deriv)
    class(avg_matl_prop), intent(in) :: this
    real(r8), intent(in) :: w(:), state(:)
    integer, intent(in) :: n
    real(r8), intent(out) :: deriv
    integer :: k
    real(r8) :: tmp1, tmp2
    ASSERT(size(w) == size(this%matl))
    tmp2 = 0.0_r8
    do k = 1, size(w)
      if (w(k) > 0) then
        call this%matl(k)%prop%compute_deriv(state, n, tmp1)
        tmp2 = tmp2 + w(k)*tmp1
      end if
    end do
    deriv = tmp2
  end subroutine compute_deriv_1

  subroutine compute_deriv_2(this, w, state, n, deriv)
    class(avg_matl_prop), intent(in) :: this
    real(r8), intent(in) :: w(:,:), state(:,:)
    integer, intent(in) :: n
    real(r8), intent(out) :: deriv(:)
    integer :: k, j
    real(r8) :: tmp1
    ASSERT(size(w,2) == size(this%matl))
    ASSERT(size(w,1) == size(state,2))
    ASSERT(size(deriv) == size(state,2))
    deriv = 0.0_r8
    do k = 1, size(w,2)
      do j = 1, size(w,1)
        if (w(j,k) > 0) then
          call this%matl(k)%prop%compute_deriv(state(:,j), n, tmp1)
          deriv(j) = deriv(j) + w(j,k)*tmp1
        end if
      end do
    end do
  end subroutine compute_deriv_2

end module avg_matl_prop_type
