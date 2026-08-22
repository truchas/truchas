subroutine test_vector_func(x, p, r) bind(c)

  use,intrinsic :: iso_c_binding, only: c_double

  real(c_double), intent(in) :: x(*), p(*)
  real(c_double), intent(out) :: r(*)

  r(1) = x(2) + p(1)
  r(2) = x(3) * p(2)

end subroutine test_vector_func
