!!
!! gs-f.F
!!
!! Serial-only dummy routines that correspond to the MPI parallel
!! C functions from par-src/gath-scatt.
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

!!!! PROCEDURES FROM gs-util-c.c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pgslib_gs_init_trace_c(GS_Trace) bind(c)
  use,intrinsic :: iso_c_binding, only: c_ptr
  type(c_ptr) :: GS_Trace
end subroutine

subroutine pgslib_gs_release_trace_c(GS_Trace) bind(c)
  use,intrinsic :: iso_c_binding, only: c_ptr
  type(c_ptr) :: GS_Trace
end subroutine

subroutine pgslib_setup_n_duplicate_c (N_Duplicate, GS_Trace) bind(c)
  use,intrinsic :: iso_c_binding, only: c_ptr, c_int
  integer(c_int) :: N_Duplicate
  type(c_ptr) :: GS_Trace
  N_Duplicate = 0
end subroutine

subroutine pgslib_prep_supplement_c (N_Supplement, PEs, GS_Trace) bind(c)
  use,intrinsic :: iso_c_binding, only: c_ptr, c_int
  integer(c_int) :: N_Supplement, PEs(*)
  type(c_ptr) :: GS_Trace
end subroutine

subroutine pgslib_setup_duplicate_buffer_c (GS_Trace) bind(c)
  use,intrinsic :: iso_c_binding, only: c_ptr
  type(c_ptr) :: GS_Trace
end subroutine

subroutine pgslib_trace_degree_c(Scatter_Degree, Gather_Degree, GS_Trace) bind(c)
  use,intrinsic :: iso_c_binding, only: c_ptr
  integer :: Scatter_Degree, Gather_Degree
  type(c_ptr) :: GS_Trace
  ! For serial code there should not be any off-PE moves, so return 0
  Scatter_Degree = 0
  Gather_Degree = 0
end subroutine

!!!! PROCEDURES FROM gather-*-c.c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pgslib_gather_buf_int_c (Supplement, Duplicate, BlockSize, GS_Trace) bind(c)
  use,intrinsic :: iso_c_binding, only: c_ptr, c_int
  integer(c_int) :: Supplement(*), Duplicate(*)
  integer(c_int) :: BlockSize
  type(c_ptr) :: GS_Trace
  ! For the serial code, there shouldn't be any off-PE moves, so don't need to do anything.
end subroutine

subroutine pgslib_gather_buf_float_c (Supplement, Duplicate, BlockSize, GS_Trace) bind(c)
  use,intrinsic :: iso_c_binding, only: c_ptr, c_int, c_float
  real(c_float) :: Supplement(*), Duplicate(*)
  integer(c_int) :: BlockSize
  type(c_ptr) :: GS_Trace
  ! For the serial code, there shouldn't be any off-PE moves, so don't need to do anything.
end subroutine

subroutine pgslib_gather_buf_double_c (Supplement, Duplicate, BlockSize, GS_Trace) bind(c)
  use,intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double
  real(c_double) :: Supplement(*), Duplicate(*)
  integer(c_int) :: BlockSize
  type(c_ptr) :: GS_Trace
end subroutine

subroutine pgslib_gather_buf_log_c (Supplement, Duplicate, BlockSize, GS_Trace) bind(c)
  use,intrinsic :: iso_c_binding, only: c_ptr, c_int
  logical(c_int) :: Supplement(*), Duplicate(*) ! SEE NOTES 1 AND 2
  integer(c_int) :: BlockSize
  type(c_ptr) :: GS_Trace
end subroutine

!!!! PROCEDURES FROM scatter-*-c.c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pgslib_scatter_buf_int_c (Duplicate, Supplement, BlockSize, GS_Trace) bind(c)
  use,intrinsic :: iso_c_binding, only: c_ptr, c_int
  integer(c_int) :: Duplicate(*), Supplement(*)
  integer(c_int) :: BlockSize
  type(c_ptr) :: GS_Trace
end subroutine

subroutine pgslib_scatter_buf_float_c (Duplicate, Supplement, BlockSize, GS_Trace) bind(c)
  use,intrinsic :: iso_c_binding, only: c_ptr, c_int, c_float
  real(c_float) :: Duplicate(*), Supplement(*)
  integer(c_int) :: BlockSize
  type(c_ptr) :: GS_Trace
end subroutine

subroutine pgslib_scatter_buf_double_c (Duplicate, Supplement, BlockSize, GS_Trace) bind(c)
  use,intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double
  real(c_double) :: Duplicate(*), Supplement(*)
  integer(c_int) :: BlockSize
  type(c_ptr) :: GS_Trace
end subroutine

subroutine pgslib_scatter_buf_log_c (Duplicate, Supplement, BlockSize, GS_Trace) bind(c)
  use,intrinsic :: iso_c_binding, only: c_ptr, c_int
  logical(c_int) :: Duplicate(*), Supplement(*) ! SEE NOTES 1 AND 2
  integer(c_int) :: BlockSize
  type(c_ptr) :: GS_Trace
end subroutine
