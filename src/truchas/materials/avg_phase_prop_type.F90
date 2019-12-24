!!
!! AVG_PHASE_PROP_TYPE
!!
!! A freqently encountered pattern is the evaluation of a linear combination
!! function of the form $f(x, w) = \sum_{i=1}^m w_i f_i(x)$ for a fixed set
!  of functions $f_i$. This module defines the type AVG_PHASE_PROP that
!! encapsulates this linear combination for the case of SCALAR_FUNC class
!! phase property functions
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

module avg_phase_prop_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use scalar_func_class
  implicit none
  private

  type :: scalar_func_box
    class(scalar_func), allocatable :: func
  end type

  type, public :: avg_phase_prop
    !private ! public only for phase_model:alloc_avg_phase_prop
    type(scalar_func_box), allocatable :: phase(:)
  contains
    generic   :: compute_value => compute_value_1
    procedure, private :: compute_value_1
  end type

contains

  subroutine compute_value_1(this, w, state, value)
    class(avg_phase_prop), intent(in) :: this
    real(r8), intent(in) :: w(:), state(:)
    real(r8), intent(out) :: value
    integer :: n
    real(r8) :: tmp2
    ASSERT(size(w) >= size(this%phase))
    tmp2 = 0.0_r8
    do n = 1, size(this%phase)
      if (w(n) /= 0) tmp2 = tmp2 + w(n)*this%phase(n)%func%eval(state)
    end do
    value = tmp2
  end subroutine

end module avg_phase_prop_type
