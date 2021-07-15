!!
!! PHASE_CHANGE_CLASS
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

module phase_change_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: phase_change
    real(r8), allocatable :: latent_heat
  contains
    procedure(solid_frac), deferred :: solid_frac
    procedure(solidus_temp), deferred :: solidus_temp
    procedure(liquidus_temp), deferred :: liquidus_temp
    procedure :: ref_liquidus_temp
    procedure :: write_solid_frac_plotfile
  end type

  abstract interface
    function solid_frac(this, temp) result(fs)
    import r8, phase_change
    class(phase_change), intent(in) :: this
    real(r8), intent(in) :: temp
    real(r8) :: fs
    end function

    !! Solid fraction is 1 for T <= solidus_temp and < 1 for T > solidus_temp
    pure function solidus_temp(this)
    import r8, phase_change
    class(phase_change), intent(in) :: this
    real(r8) :: solidus_temp
    end function

    !! Solid fraction is 0 for T >= liquidus_temp and > 0 for T < liquidus_temp
    pure function liquidus_temp(this)
    import r8, phase_change
    class(phase_change), intent(in) :: this
    real(r8) :: liquidus_temp
    end function
  end interface

contains

  !! The latent heat is the difference in specific enthalpy between liquid and
  !! solid phases at this temperature. Normally this would equal liquidus_temp(),
  !! but in some cases the user-specified phase change interval gets spread due
  !! to internal numerical smoothing. Extensions should override this if needed.
  pure function ref_liquidus_temp(this)
    class(phase_change), intent(in) :: this
    real(r8) :: ref_liquidus_temp
    ref_liquidus_temp = this%liquidus_temp()
  end function

  subroutine write_solid_frac_plotfile(this, filename, digits, npoints, iostat)

    class(phase_change), intent(in) :: this
    character(*), intent(in) :: filename
    integer, intent(in) :: digits, npoints
    integer, intent(out) :: iostat

    integer :: lun, j
    character(12) :: fmt
    real(r8) :: Tsol, Tliq, dT, T

    ASSERT(npoints > 0)
    ASSERT(digits > 1)

    Tsol = this%solidus_temp()
    Tliq = this%liquidus_temp()

    open(newunit=lun,file=filename,status='replace',action='write',iostat=iostat)
    if (iostat /= 0) return
    write(fmt,'("(*(es",i0,".",i0,"))")') digits+7, digits-1

    dT = (Tliq - Tsol) / npoints

    do j = 0, npoints
      T = Tsol + j*dT
      write(lun,fmt) T, this%solid_frac(T)
    end do

    close(lun)

  end subroutine write_solid_frac_plotfile

end module phase_change_class

