subroutine axisymmetric_velocity(x, p, r) bind(c)

  use,intrinsic :: iso_fortran_env, only: r8 => real64

  real(r8), intent(in) :: x(*), p(*)
  real(r8), intent(out) :: r(*)

  real(r8), parameter :: pi = acos(-1.0_r8)

  r(1) = 0.5_r8 * x(2) * cos(pi * x(1))
  r(2) = -x(3) * cos(pi * x(1))

end subroutine axisymmetric_velocity
