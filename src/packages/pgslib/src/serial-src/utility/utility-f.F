!!
!! utility-f.F
!!
!! Serial-only Fortran replacements for the MPI parallel C functions from
!! par-src/utility/utility-c.c.
!!
!! This is a modification of the original code that uses the C
!! interoperability features of Fortran 2003.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! November 2014
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NOTES
!!
!! 1. Many of the C functions take C strings (null-terminated) as arguments.
!! These are interoperable with Fortran character(len=1,kind=c_char) rank-1
!! assumed-size arrays and that is how those arguments are declared in the
!! the bind(c) interfaces defined for the functions.  Fortran allows scalar
!! character variables of arbitrary length to be used as actual arguments
!! for such dummy arrays, and the Fortran client code uses such variables as
!! it is much more convenient than working with arrays of single-character
!! elements.  We would like to do the same thing here in these replacement
!! Fortran procedures (as it was originally), but such dummy arguments are
!! not  allowed in bind(c) procedures.  Thus we are forced to implement the
!! character dummy arguments exactly as declared in the interface blocks.
!! This has required a bit unpleasant new code (pgslib_c_string_to_f).
!!

module pgslib_c_string_to_f
  use,intrinsic :: iso_c_binding, only: c_char
  implicit none
contains
  subroutine c_string_to_f (c_string, f_string)
    character(kind=c_char), intent(in) :: c_string(*)
    character(:,kind=c_char), allocatable, intent(out) :: f_string
    integer :: n
    n = 0
    do while (c_string(n+1) /= achar(0))
      n = n + 1
    end do
    allocate(character(n)::f_string)
    do n = 1, len(f_string)
      f_string(n:n) = c_string(n)
    end do
  end subroutine
end module pgslib_c_string_to_f

subroutine pgslib_initialize_c (nPE, thisPE, IO_ROOT_PE, File_Per_PE, File_Prefix) bind(c)

  use pgslib_c_string_to_f
  use,intrinsic :: iso_c_binding, only: c_int, c_char
  implicit none

  integer(c_int), intent(out)   :: nPE, thisPE
  integer(c_int), intent(inout) :: IO_ROOT_PE
  integer(c_int), intent(in)    :: File_Per_PE
  character(kind=c_char), intent(in) :: File_Prefix(*)

  ! File names are setup and used by output procedures.
  character(2048) :: F_OUT_NAME, F_ERR_NAME
  integer :: ERR_UNIT, OUT_UNIT
  common /F_NAMES/ F_OUT_NAME, F_ERR_NAME
  common /F_UNITS/ ERR_UNIT, OUT_UNIT

  character(:,kind=c_char), allocatable :: prefix
  logical, save :: been_initialized = .false.

  ! This routine may get called more than once, but can do init only once.
  if (.not.been_initialized) then
    been_initialized = .true.
    call c_string_to_f (File_Prefix, prefix)
    F_OUT_NAME = trim(prefix) // '-out'
    OUT_UNIT   = 29
    F_ERR_NAME = trim(prefix) // '-err'
    ERR_UNIT   = 30
  end if

  ! NOTE: This routine must return valid values on all calls.
  nPE = 1
  thisPE = 0
  IO_ROOT_PE = 0

end subroutine

subroutine pgslib_finalize_c() bind(c)
end subroutine

subroutine pgslib_error_c (Estring) bind(c)

  use pgslib_c_string_to_f
  use,intrinsic :: iso_c_binding, only: c_char
  implicit none

  character(kind=c_char), intent(in) :: Estring(*)

  ! File names are setup and used by output procedures.
  character(2048) :: F_OUT_NAME, F_ERR_NAME
  integer :: ERR_UNIT, OUT_UNIT
  common /F_NAMES/ F_OUT_NAME, F_ERR_NAME
  common /F_UNITS/ ERR_UNIT, OUT_UNIT

  character(:,kind=c_char), allocatable :: f_string
  logical, save :: Err_Opened = .false.
  integer :: ioerror

  if (.not.Err_Opened) then
    open(unit=ERR_UNIT, file=trim(F_ERR_NAME), iostat=ioerror)
    if (ioerror /= 0) then
      write(*,*) 'Error opening error output file:', trim(F_ERR_NAME)
    end if
    Err_Opened = .true.
  end if
  call c_string_to_f (Estring, f_string)
  write(ERR_UNIT,*) 'ERROR: ', f_string

end subroutine

subroutine pgslib_output_c (message) bind(c)

  use pgslib_c_string_to_f
  use,intrinsic :: iso_c_binding, only: c_char
  implicit none

  character(kind=c_char), intent(in) :: message(*)

  ! File names are setup and used by output procedures.
  character(2048) :: F_OUT_NAME, F_ERR_NAME
  integer :: ERR_UNIT, OUT_UNIT
  common /F_NAMES/ F_OUT_NAME, F_ERR_NAME
  common /F_UNITS/ ERR_UNIT, OUT_UNIT

  character(:,kind=c_char), allocatable :: f_string
  logical, save :: Out_Opened = .false.
  integer :: ioerror

  if (.not.Out_Opened) then
    open(unit=OUT_UNIT, file=trim(F_OUT_NAME), iostat=ioerror)
    if (ioerror /= 0) then
      write(*,*) 'Error opening standard output file:',TRIM(F_OUT_NAME)
    end if
    Out_Opened = .true.
  end if
  call c_string_to_f (message, f_string)
  write(OUT_UNIT,*) 'Output: ', f_string

end subroutine

subroutine pgslib_flush_output_c() bind(c)
end subroutine

subroutine pgslib_abort_c() bind(c)
  stop
end subroutine

subroutine pgslib_barrier_c() bind(c)
end subroutine

subroutine pgslib_barrier_time_c(bt) bind(c)
  use,intrinsic :: iso_c_binding, only: c_float
  real(c_float) :: bt
  bt = 0.0_c_float
end subroutine

subroutine pgslib_sr_time_c(t) bind(c)
  use,intrinsic :: iso_c_binding, only: c_float
  real(c_float) :: t
  t = 0.0_c_float
end subroutine

! in the C version, these three routines are used to get the command
! line from C back to Fortran

subroutine pgslib_get_argc (argc) bind (c)
  use,intrinsic :: iso_c_binding, only: c_int
  integer(c_int) :: argc
  argc = 1 + command_argument_count()
end subroutine

subroutine pgslib_get_argv (i, string_l, string) bind (c)
  use,intrinsic :: iso_c_binding, only: c_int, c_char
  integer(c_int) :: i, string_l, n
  character(kind=c_char) :: string(string_l)
  character(len=string_l,kind=c_char) :: arg
  call get_command_argument (i, value=arg)
  string = ' '
  do n = 1, string_l
    string(n) = arg(n:n)
  end do
end subroutine

subroutine pgslib_cleanup_argv() bind(c)
end subroutine
