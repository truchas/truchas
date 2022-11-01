!!
!! SMOOTH_PHASE_CHANGE_TYPE
!!
!! This extension of the PHASE_CHANGE abstract class gives the solid fraction
!! as a smooth C1/C2 temperature polynomial. This has no basis as a physics
!! model, and is best seen either as an approximation to an isothermal phase
!! change, which the heat conduction discretization is not able to treat, or
!! as a purely mathematical model.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! December 2019
!!
!! NB: This could have been implemented as a table with Akima smoothing.
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

module smooth_phase_change_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use phase_change_class
  implicit none
  private

  public :: alloc_smooth_phase_change

  type, extends(phase_change), public :: smooth_phase_change
    private
    real(r8) :: tlo, thi
  contains
    procedure :: init
    procedure :: solid_frac
    procedure :: solidus_temp
    procedure :: liquidus_temp
  end type

contains

  subroutine alloc_smooth_phase_change(this, tlo, thi)
    class(phase_change), allocatable, intent(out) :: this
    real(r8), intent(in) :: tlo, thi
    type(smooth_phase_change), allocatable :: pc
    allocate(pc)
    call pc%init(tlo, thi)
    call move_alloc(pc, this)
  end subroutine

  subroutine init(this, tlo, thi)
    class(smooth_phase_change), intent(out) :: this
    real(r8), intent(in) :: tlo, thi
    this%tlo = tlo
    this%thi = thi
  end subroutine

  pure function solidus_temp(this)
    class(smooth_phase_change), intent(in) :: this
    real(r8) :: solidus_temp
    solidus_temp = this%tlo
  end function

  pure function liquidus_temp(this)
    class(smooth_phase_change), intent(in) :: this
    real(r8) :: liquidus_temp
    liquidus_temp = this%thi
  end function

  function solid_frac(this, temp) result(fs)
    use,intrinsic :: iso_fortran_env, only: i8 => int64
    class(smooth_phase_change), intent(in) :: this
    real(r8), intent(in) :: temp
    real(r8) :: fs
    integer, parameter :: D = digits(fs)
    real(r8) :: z
    if (temp <= this%tlo) then
      fs = 1
    else if (temp >= this%thi) then
      fs = 0
    else
      z = (this%thi - temp) / (this%thi - this%tlo)
      fs = (z*z)*(3 - 2*z) ! C1 Hermite cubic
      !fs = (z*z*z)*(10 - z*(15 - 6*z))   ! C2
      fs = scale(real(int(scale(fs,D),i8),r8),-D) ! See Note 1
    end if
  end function

end module smooth_phase_change_type
