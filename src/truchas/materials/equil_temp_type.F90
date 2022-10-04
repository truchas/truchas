!!
!! EQUIL_TEMP_TYPE
!!
!! This class provides a method for computing the equilibrium temperature of
!! a body formed from a collection of material phases at possibly different
!! temperatures.
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
    type(avg_phase_prop) :: H_of_T
    real(r8), pointer :: w(:) => null()
  contains
    procedure :: g
  end type

  type, public :: equil_temp
    private
    type(invf) :: T_of_H
    real(r8), allocatable :: v(:) ! persistent workspace for compute
  contains
    procedure :: init
    procedure :: compute
  end type

contains

  subroutine init(this, pids, model, stat, errmsg)
    use material_model_type
    class(equil_temp), intent(out) :: this
    integer, intent(in) :: pids(:)
    type(material_model), intent(in) :: model
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call this%T_of_H%H_of_T%init('enthalpy', pids, model, stat, errmsg)
    if (stat /= 0) return
    allocate(this%v(size(pids)))
    call this%T_of_H%init(eps=0.0_r8)
  end subroutine

  subroutine compute(this, w, temps, temp)

    class(equil_temp), intent(inout) :: this
    real(r8), intent(in), target :: w(:)
    real(r8), intent(in) :: temps(:)
    real(r8), intent(out) :: temp

    integer :: n, stat
    real(r8) :: tmin, tmax, H, Hn

    ASSERT(size(w) == size(this%v))
    ASSERT(size(w) == size(temps))

    tmin = minval(temps, mask=(w > 0))
    tmax = maxval(temps, mask=(w > 0))
    INSIST(tmin <= tmax)

    if (tmin == tmax) then
      temp = tmin
      return
    end if

    !! Total enthalpy. This is a bit of a hack. We need to evaluate the
    !! enthalpy of the individual materials at *different* temperatures,
    !! but only have access to the average enthalpy function.
    H = 0
    do n = 1, size(w)
      if (w(n) > 0) then
        this%v = 0
        this%v(n) = w(n)
        call this%T_of_H%H_of_T%compute_value(this%v, [temps(n)], Hn)
        H = H + Hn
      end if
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
