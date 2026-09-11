function rt_interface(x, p) result(f) bind(c)
  use,intrinsic :: iso_fortran_env, only: r8 => real64
  real(r8), intent(in) :: x(*), p(*)
  real(r8) :: f
  real(r8), parameter :: TWOPI = 8*atan(1.0_r8)
  f = x(2) - (2 + p(1)*cos(TWOPI*x(1)))
end function
