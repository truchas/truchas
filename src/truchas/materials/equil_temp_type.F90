!!
!! EQUIL_TEMP_TYPE
!!
!! This class provides a method for computing the equilibrium temperature of
!! a body formed from a collection of material phases at possibly different
!! temperatures.
!!
!! See MATERIAL_MODEL%ALLOC_EQUIL_TEMP for a factory method for this type.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! May 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Given a set of weights $w_j$, phase IDs $p_j$, and temperatures $T_j$,
!! $j=1,\ldots n$, this solves
!! \begin{equation*}
!!   \sum_{j=1}^n w_j H_{p_j}(T) = \sum_{j=1}^n w_j H_{p_j}(T_j)$$
!! \end{equation*}
!! for the equilbrium temperature $T$, where $H_p$ is the enthalpy density of
!! phase $p$ as a function of temperature. The weights can be either volumes
!! or volume fractions.
!!

#include "f90_assert.fpp"

module equil_temp_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use inverse_func_class
  use avg_phase_prop_type
  implicit none
  private

  type, extends(inverse_func) :: invf
    type(avg_phase_prop), pointer :: H_of_T => null()
    real(r8), pointer :: w(:) => null()
  contains
    procedure :: g
  end type

  type, public :: equil_temp
    type(avg_phase_prop), allocatable :: H_of_T
    type(invf) :: T_of_H
  contains
    procedure :: compute
  end type

contains

  subroutine compute(this, w, temps, temp)

    class(equil_temp), intent(inout) :: this
    real(r8), intent(in), target :: w(:)
    real(r8), intent(in) :: temps(:)
    real(r8), intent(out) :: temp

    integer :: n, stat
     real(r8) :: tmin, tmax, H

    ASSERT(allocated(this%H_of_T))
    ASSERT(size(w) == size(this%H_of_T%phase))
    ASSERT(size(w) == size(temps))

    tmin = minval(temps, mask=(w > 0))
    tmax = maxval(temps, mask=(w > 0))
    INSIST(tmin <= tmax)

    if (tmin == tmax) then
      temp = tmin
      return
    end if

    !! Total enthalpy
    H = 0
    do n = 1, size(w)
      if (w(n) > 0) H = H + w(n)*this%H_of_T%phase(n)%func%eval([temps(n)])
    end do

    !! Equilibrium temperature yielding the same total enthalpy
    this%T_of_H%w => w
    call this%T_of_H%compute(H, tmin, tmax, temp, stat)
    INSIST(stat == 0)

  end subroutine

  function g(this, x) result(H)
    class(invf), intent(in) :: this
    real(r8), intent(in) :: x
    real(r8) :: H
    call this%H_of_T%compute_value(this%w, [x], H)
  end function

end module equil_temp_type
