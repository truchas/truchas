!!
!! TABULAR_PHASE_CHANGE
!!
!! This extension of the PHASE_CHANGE abstract base class implements the solid
!! fraction function using a user-specified temperature-solid fraction table.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! December 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Note 1: This expression converts the real floating point value 0 <= fs <= 1
!! to a fixed point real value relative to 64-bit real 1. The idea here is that
!! we would like phase fraction addition/subtraction to be exact, and this can
!! be accomplished by working in fixed point arithmetic (or equivalently in an
!! integer space).
!!

#include "f90_assert.fpp"

module tabular_phase_change_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use phase_change_class
  use tabular_scalar_func_type
  implicit none
  private

  public :: alloc_tabular_phase_change

  type, extends(phase_change), public :: tabular_phase_change
    private
    real(r8) :: tlo, thi
    type(tabular_scalar_func) :: fs
  contains
    procedure :: init
    procedure :: solidus_temp
    procedure :: liquidus_temp
    procedure :: solid_frac
  end type tabular_phase_change

contains

  subroutine alloc_tabular_phase_change(this, temp, fsol, errmsg)
    class(phase_change), allocatable, intent(out) :: this
    real(r8), intent(in) :: temp(:), fsol(:)
    character(:), allocatable, intent(out) :: errmsg
    type(tabular_phase_change), allocatable :: pc
    allocate(pc)
    call pc%init(temp, fsol, errmsg)
    if (.not.allocated(errmsg)) call move_alloc(pc, this)
  end subroutine alloc_tabular_phase_change

  pure function solidus_temp(this)
    class(tabular_phase_change), intent(in) :: this
    real(r8) :: solidus_temp
    solidus_temp = this%tlo
  end function

  pure function liquidus_temp(this)
    class(tabular_phase_change), intent(in) :: this
    real(r8) :: liquidus_temp
    liquidus_temp = this%thi
  end function

  function solid_frac(this, temp) result(fs)
    use,intrinsic :: iso_fortran_env, only: i8 => int64
    class(tabular_phase_change), intent(in) :: this
    real(r8), intent(in) :: temp
    real(r8) :: fs
    integer, parameter :: D = digits(fs)
    if (temp <= this%tlo) then
      fs = 1
    else if (temp >= this%thi) then
      fs = 0
    else
      fs = this%fs%eval([temp])
      fs = min(1.0_r8, max(0.0_r8, fs))
      fs = scale(real(int(scale(fs,D),i8),r8),-D) ! See Note 1
    end if
  end function solid_frac

  subroutine init(this, temp, fsol, errmsg)

    class(tabular_phase_change), intent(out) :: this
    real(r8), intent(in) :: temp(:), fsol(:)
    character(:), allocatable, intent(out) :: errmsg

    integer :: n, j
    logical :: decreasing
    real(r8) :: dx
    real(r8), allocatable :: x(:), y(:)

    if (size(temp) /= size(fsol) .or. size(temp) < 2) then
      errmsg = 'malformed table'
      return
    end if

    n = size(temp)

    !! We can accomodate a table with either increasing or decreasing temps
    decreasing = (temp(2) < temp(1))

    !! Check for ordered temperature values.
    if (decreasing) then
      do j = 2, n
        if (temp(j) >= temp(j-1)) exit
      end do
    else
      do j = 2, n
        if (temp(j) <= temp(j-1)) exit
      end do
    end if
    if (j <= n) then
      errmsg = 'unordered temperature values'
      return
    end if

    !! Check for strictly monotonic solid fraction values from 0 to 1.
    if (decreasing) then
      do j = 2, n
        if (fsol(j) <= fsol(j-1)) exit
      end do
      if (j <= n) then
        errmsg = 'solid fraction not monotonically increasing with decreasing temperature'
      else if (fsol(1) /= 0 .or. fsol(n) /= 1) then
        errmsg = 'solid fraction must be 0 and 1 at the extreme temperatures'
      end if
    else
      do j = 2, n
        if (fsol(j) >= fsol(j-1)) exit
      end do
      if (j <= n) then
        errmsg = 'solid fraction not monotonically decreasing with increasing temperature'
      else if (fsol(1) /= 1 .or. fsol(n) /= 0) then
        errmsg = 'solid fraction must be 0 and 1 at the extreme temperatures'
      end if
    end if
    if (allocated(errmsg)) return

    if (decreasing) then
      this%tlo = temp(n)
      this%thi = temp(1)
    else
      this%tlo = temp(1)
      this%thi = temp(n)
    end if

    !! To prepare for the possiblilty of using Akima smoothing, we pad the
    !! table on either end, continuing the end point solid fraction. This
    !! will ensure that the smoothed solid fraction comes into the input
    !! end points with 0 slope, and is constant (0 or 1) outside the phase
    !! change temperature interval.
    !!
    !! It is not obvious that the tabular_scalar_func type allows x values
    !! in decreasing order and I'm too lazy to figure it out, so I'll just
    !! reverse the order here as needed.

    allocate(x(-1:n+2), y(-1:n+2))
    if (decreasing) then
      x(n:1:-1) = temp
      y(n:1:-1) = fsol
    else
      x(1:n) = temp
      y(1:n) = fsol
    end if

    dx = x(2) - x(1)
    x(0)  = x(1) - dx; y(0)  = y(1)
    x(-1) = x(0) - dx; y(-1) = y(1)

    dx = x(n) - x(n-1)
    x(n+1) = x(n)   + dx; y(n+1) = y(n)
    x(n+2) = x(n+1) + dx; y(n+2) = y(n)

#ifdef INTEL_BUG20191229
    ! this avoids an incredibly nasty wrong-results bug with Intel
    this%fs = tabular_scalar_func(x(-1:n+2), y(-1:n+2), smooth=.true.)
#else
    this%fs = tabular_scalar_func(x, y, smooth=.true.)
#endif

  end subroutine init

end module tabular_phase_change_type
