subroutine vortex_velocity(x, p, r) bind(c)

  use,intrinsic :: iso_fortran_env, only: r8 => real64

  real(r8), intent(in) :: x(*), p(*)
  real(r8), intent(out) :: r(*)

  real(r8), parameter :: pi = acos(-1.0_r8)
  real(r8) :: factor

  factor = cos(pi * x(1) / 8.0_r8)
  r(1) = -2.0_r8 * sin(pi*x(2))**2 * sin(pi*x(3)) * cos(pi*x(3)) * factor
  r(2) =  2.0_r8 * sin(pi*x(3))**2 * sin(pi*x(2)) * cos(pi*x(2)) * factor

end subroutine vortex_velocity
