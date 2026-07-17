module enthalpy_advector_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: enthalpy_advector
  contains
    procedure(p1), deferred :: get_advected_enthalpy1
    procedure(p2), deferred :: get_advected_enthalpy2
    generic :: get_advected_enthalpy => get_advected_enthalpy1, get_advected_enthalpy2
  end type enthalpy_advector

  abstract interface
    subroutine p1(this, tcell, dq)
      import
      class(enthalpy_advector), intent(in) :: this  !TODO: inout?
      real(r8), intent(in) :: tcell(:)
      real(r8), intent(out) :: dq(:)
    end subroutine
    subroutine p2(this, tcell, dq, tmin, tmax)
      import
      class(enthalpy_advector), intent(in) :: this  !TODO: inout?
      real(r8), intent(in) :: tcell(:)
      real(r8), intent(out) :: dq(:), tmin(:), tmax(:)
    end subroutine
  end interface

end module enthalpy_advector_class
