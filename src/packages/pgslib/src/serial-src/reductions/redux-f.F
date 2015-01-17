!!
!! redux-f.F
!!
!! Serial-only dummy routines that correspond to the MPI parallel
!! C functions from par-src/reductions/redux-c.c.
!!
!! This is a modern rewrite of the original code by Robert Ferrell
!! using the C interoperability features of Fortran 2003.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! November 2014
!!

subroutine pgslib_global_min_int_c (MinA) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int
  integer(c_int), intent(inout):: MinA
end subroutine

subroutine pgslib_global_min_float_c (MinA) bind(c)
  use,intrinsic :: iso_c_binding, only: c_float
  real(c_float), intent(inout):: MinA
end subroutine

subroutine pgslib_global_min_double_c (MinA) bind(c)
  use,intrinsic :: iso_c_binding, only: c_double
  real(c_double), intent(inout):: MinA
end subroutine

subroutine pgslib_global_max_int_c (MaxA) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int
  integer(c_int), intent(inout):: MaxA
end subroutine

subroutine pgslib_global_max_float_c (MaxA) bind(c)
  use,intrinsic :: iso_c_binding, only: c_float
  real(c_float), intent(inout):: MaxA
end subroutine

subroutine pgslib_global_max_double_c (MaxA) bind(c)
  use,intrinsic :: iso_c_binding, only: c_double
  real(c_double), intent(inout):: MaxA
end subroutine

subroutine pgslib_global_sum_int_c (SumA) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int
  integer(c_int), intent(inout):: SumA
end subroutine

subroutine pgslib_global_sum_float_c (SumA) bind(c)
  use,intrinsic :: iso_c_binding, only: c_float
  real(c_float), intent(inout):: SumA
end subroutine

subroutine pgslib_global_sum_double_c (SumA) bind(c)
  use,intrinsic :: iso_c_binding, only: c_double
  real(c_double), intent(inout):: SumA
end subroutine

subroutine pgslib_global_all_log_c (AllA) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int
  integer(c_int), intent(inout) :: AllA
end subroutine

subroutine pgslib_global_any_log_c (AnyA) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int
  integer(c_int), intent(inout) :: AnyA
end subroutine

subroutine pgslib_global_minloc_int_c (MinV, GlobalIndex) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int
  integer(c_int), intent(inout) :: GlobalIndex
  integer(c_int), intent(inout) :: MinV
end subroutine

subroutine pgslib_global_minloc_float_c (MinV, GlobalIndex) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int, c_float
  integer(c_int), intent(inout) :: GlobalIndex
  real(c_float),  intent(inout) :: MinV
end subroutine

subroutine pgslib_global_minloc_double_c (MinV, GlobalIndex) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int, c_double
  integer(c_int), intent(inout) :: GlobalIndex
  real(c_double), intent(inout) :: MinV
end subroutine

subroutine pgslib_global_maxloc_int_c (MaxV, GlobalIndex) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int
  integer(c_int), intent(inout) :: GlobalIndex
  integer(c_int), intent(inout) :: MaxV
end subroutine

subroutine pgslib_global_maxloc_float_c (MaxV, GlobalIndex) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int, c_float
  integer(c_int), intent(inout) :: GlobalIndex
  real(c_float),  intent(inout) :: MaxV
end subroutine

subroutine pgslib_global_maxloc_double_c (MaxV, GlobalIndex) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int, c_double
  integer(c_int), intent(inout) :: GlobalIndex
  real(c_double), intent(inout) :: MaxV
end subroutine
