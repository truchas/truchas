subroutine te01(x, p, fx) bind(c)
  use,intrinsic :: iso_fortran_env, only: r8 => real64
  real(r8) :: x(0:4) ! t, x, y, z
  real(r8) :: p(3) ! a, h0, E0
  complex(r8) :: fx(3)
  real(r8), parameter :: PI = 3.1415926535897932385_r8
  fx = 0.0_r8
  fx(2)%im = 2*p(2)*p(3)*cos(PI*x(1)/p(1))
end subroutine
