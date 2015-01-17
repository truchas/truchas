!!
!! input-f.F
!!
!! Serial-only dummy routines that correspond to the MPI parallel
!! C functions from par-src/io/io-c.c
!!
!! This is a modern rewrite of the original code by Robert Ferrell
!! using the C interoperability features of Fortran 2003.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! November 2014
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NOTES
!!
!! 1. The lower level C component (which this emulates for the serial case)
!! uses ints to represent logical values.  However, Fortran logical variables
!! are not interoperable C int variables, and so we must resort to old-school
!! hackery to interface Fortran logical arguments to C int arguments.  This
!! involves two key aspects: the data size and the representation of logical
!! values.  Here we make the assumption (potentially not portable) that integer
!! and logical types with the same type parameter are the same size; that is,
!! logical(c_int) and integer(c_int) are the same size.
!!
!! 2. To match the corresponding C code we should declare these variables as
!! integers.  However the NAG compiler with -C run-time checking will detect
!! that the corresponding actual arguments are logical(c_int).  Thus we
!! declare them as logical(c_int) to avoid those spurious errors.
!!

!!!! BROADCAST PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pgslib_bcast_int_scalar_c (scalar) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int
  integer(c_int) :: scalar
end subroutine

subroutine pgslib_bcast_float_scalar_c (scalar) bind(c)
  use,intrinsic :: iso_c_binding, only: c_float
  real(c_float) :: scalar
end subroutine

subroutine pgslib_bcast_double_scalar_c (scalar) bind(c)
  use,intrinsic :: iso_c_binding, only: c_double
  real(c_double) :: scalar
end subroutine

subroutine pgslib_bcast_log_scalar_c (scalar) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int
  logical(c_int) :: scalar  ! SEE NOTES 1 and 2
end subroutine

subroutine pgslib_bcast_int_vector_c (vector, vec_len) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int
  integer(c_int) :: vector(*)
  integer(c_int), intent(in) :: vec_len
end subroutine

subroutine pgslib_bcast_float_vector_c (vector, vec_len) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int, c_float
  real(c_float) :: vector(*)
  integer(c_int), intent(in) :: vec_len
end subroutine

subroutine pgslib_bcast_double_vector_c (vector, vec_len) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int, c_double
  real(c_double) :: vector(*)
  integer(c_int), intent(in) :: vec_len
end subroutine

subroutine pgslib_bcast_log_vector_c (vector, vec_len) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int
  logical(c_int) :: vector(*)  ! SEE NOTES 1 and 2
  integer(c_int), intent(in) :: vec_len
end subroutine

subroutine pgslib_bcast_char_vector_c (vector, vec_len) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int, c_char
  character(kind=c_char) :: vector(*)
  integer(c_int), intent(in) :: vec_len
end subroutine

!!!! DISTRIBUTE PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pgslib_dist_int_scalar_c (scalar_out, scalarv_in) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int
  integer(c_int), intent(out) :: scalar_out
  integer(c_int), intent(in)  :: scalarv_in(*)
  scalar_out = scalarv_in(1)
end subroutine

subroutine pgslib_dist_float_scalar_c (scalar_out, scalarv_in) bind(c)
  use,intrinsic :: iso_c_binding, only: c_float
  real(c_float), intent(out) :: scalar_out
  real(c_float), intent(in)  :: scalarv_in(*)
  scalar_out = scalarv_in(1)
end subroutine

subroutine pgslib_dist_double_scalar_c (scalar_out, scalarv_in) bind(c)
  use,intrinsic :: iso_c_binding, only: c_double
  real(c_double), intent(out) :: scalar_out
  real(c_double), intent(in)  :: scalarv_in(*)
  scalar_out = scalarv_in(1)
end subroutine

subroutine pgslib_dist_log_scalar_c (scalar_out, scalarv_in) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int
  logical(c_int), intent(out) :: scalar_out     ! SEE NOTES 1 and 2
  logical(c_int), intent(in)  :: scalarv_in(*)  ! SEE NOTES 1 and 2
  scalar_out = scalarv_in(1)
end subroutine

subroutine pgslib_dist_int_vector_c (vector_out, out_len, vector_in, lengths) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int
  integer(c_int), intent(out) :: vector_out(*)
  integer(c_int), intent(in)  :: vector_in(*)
  integer(c_int), intent(in)  :: out_len, lengths(*)
  vector_out(:out_len) = vector_in(:out_len)
end subroutine

subroutine pgslib_dist_float_vector_c (vector_out, out_len, vector_in, lengths) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int, c_float
  real(c_float), intent(out) :: vector_out(*)
  real(c_float), intent(in)  :: vector_in(*)
  integer(c_int), intent(in) :: out_len, lengths(*)
  vector_out(:out_len) = vector_in(:out_len)
end subroutine

subroutine pgslib_dist_double_vector_c (vector_out, out_len, vector_in, lengths) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int, c_double
  real(c_double), intent(out) :: vector_out(*)
  real(c_double), intent(in)  :: vector_in(*)
  integer(c_int), intent(in) :: out_len, lengths(*)
  vector_out(:out_len) = vector_in(:out_len)
end subroutine

subroutine pgslib_dist_log_vector_c (vector_out, out_len, vector_in, lengths) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int
  logical(c_int), intent(out) :: vector_out(*)  ! SEE NOTES 1 and 2
  logical(c_int), intent(in)  :: vector_in(*)   ! SEE NOTES 1 and 2
  integer(c_int), intent(in)  :: out_len, lengths(*)
  vector_out(:out_len) = vector_in(:out_len)
end subroutine

!!!! COLLATE PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pgslib_collate_int_scalar_c (scalarv_out, scalar_in) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int
  integer(c_int), intent(out) :: scalarv_out(*)
  integer(c_int), intent(in)  :: scalar_in
  scalarv_out(1) = scalar_in
end subroutine

subroutine pgslib_collate_float_scalar_c (scalarv_out, scalar_in) bind(c)
  use,intrinsic :: iso_c_binding, only: c_float
  real(c_float), intent(out) :: scalarv_out(*)
  real(c_float), intent(in)  :: scalar_in
  scalarv_out(1) = scalar_in
end subroutine

subroutine pgslib_collate_double_scalar_c (scalarv_out, scalar_in) bind(c)
  use,intrinsic :: iso_c_binding, only: c_double
  real(c_double), intent(out) :: scalarv_out(*)
  real(c_double), intent(in)  :: scalar_in
  scalarv_out(1) = scalar_in
end subroutine

subroutine pgslib_collate_log_scalar_c (scalarv_out, scalar_in) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int
  logical(c_int), intent(out) :: scalarv_out(*) ! SEE NOTES 1 and 2
  logical(c_int), intent(in)  :: scalar_in      ! SEE NOTES 1 and 2
  scalarv_out(1) = scalar_in
end subroutine

subroutine pgslib_collate_int_vector_c (vector_out, Lengths, vector_in, In_Length) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int
  integer(c_int), intent(out) :: vector_out(*)
  integer(c_int), intent(in)  :: vector_in(*)
  integer(c_int), intent(in)  :: Lengths(*), In_Length
  vector_out(:Lengths(1)) = vector_in(:Lengths(1))
end subroutine

subroutine pgslib_collate_float_vector_c (vector_out, Lengths, vector_in, In_Length) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int, c_float
  real(c_float),  intent(out) :: vector_out(*)
  real(c_float),  intent(in)  :: vector_in(*)
  integer(c_int), intent(in)  :: Lengths(*), In_Length
  vector_out(:Lengths(1)) = vector_in(:Lengths(1))
end subroutine

subroutine pgslib_collate_double_vector_c (vector_out, Lengths, vector_in, In_Length) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int, c_double
  real(c_double), intent(out) :: vector_out(*)
  real(c_double), intent(in)  :: vector_in(*)
  integer(c_int), intent(in)  :: Lengths(*), In_Length
  vector_out(:Lengths(1)) = vector_in(:Lengths(1))
end subroutine

subroutine pgslib_collate_log_vector_c (vector_out, Lengths, vector_in, In_Length) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int
  logical(c_int), intent(out) :: vector_out(*)  ! SEE NOTES 1 and 2
  logical(c_int), intent(in)  :: vector_in(*)   ! SEE NOTES 1 and 2
  integer(c_int), intent(in)  :: Lengths(*), In_Length
  vector_out(:Lengths(1)) = vector_in(:Lengths(1))
end subroutine

subroutine pgslib_collate_char_vector_c (vector_out, Lengths, vector_in, In_Length) bind(c)
  use,intrinsic :: iso_c_binding, only: c_int, c_char
  character(kind=c_char), intent(out) :: vector_out(*)
  character(kind=c_char), intent(in)  :: vector_in(*)
  integer(c_int), intent(in) :: Lengths(*), In_Length
  vector_out(:Lengths(1)) = vector_in(:Lengths(1))
end subroutine
