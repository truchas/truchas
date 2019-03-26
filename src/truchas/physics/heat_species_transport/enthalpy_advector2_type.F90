module enthalpy_advector2_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use enthalpy_advector_class
  implicit none
  private

  type, extends(enthalpy_advector), public :: enthalpy_advector2
  contains
    procedure :: get_advected_enthalpy1
    procedure :: get_advected_enthalpy2
  end type

contains

  !! Input/output arrays are on process cells only
  subroutine get_advected_enthalpy1(this, tcell, dq)
    use advection_module, only: compute_advected_enthalpy
    class(enthalpy_advector2), intent(in) :: this
    real(r8), intent(in) :: tcell(:)
    real(r8), intent(out) :: dq(:)
    call compute_advected_enthalpy(tcell, dq)
  end subroutine get_advected_enthalpy1

  !! Input/output arrays are on process cells only
  subroutine get_advected_enthalpy2(this, tcell, dq, tmin, tmax)
    use advection_module, only: compute_advected_enthalpy
    class(enthalpy_advector2), intent(in) :: this
    real(r8), intent(in) :: tcell(:)
    real(r8), intent(out) :: dq(:), tmin(:), tmax(:)
    call compute_advected_enthalpy(tcell, dq, tmin, tmax)
  end subroutine get_advected_enthalpy2

end module enthalpy_advector2_type
