!!
!! SERIALIZATION_TOOLS
!!
!! Provides procedures for copying intrinsic variables to/from a sequence of bytes.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! October 2015
!!
!! PROGRAMMING INTERFACE
!!
!!  The module defines two generic subroutines:
!!
!!  COPY_TO_BYTES (VAR, BUFFER, OFFSET)
!!    TYPE(INT8), INTENT(INOUT) :: BUFFER(:)
!!    INTEGER, INTENT(INOUT) :: OFFSET
!!
!!  COPY_FROM_BYTES (BUFFER, OFFSET, VAR)
!!    TYPE(INT8), INTENT(IN) :: BUFFER(:)
!!    INTEGER, INTENT(INOUT) :: OFFSET
!!
!!  VAR is a scalar variable of integer type (kinds int8, int16, int32 or int64)
!!  or real type (kinds real32 or real64). The contents of the bytes occupied by
!!  VAR are copied to/from the byte array BUFFER starting at element OFFSET+1.
!!  and OFFSET is updated and points to the last element of BUFFER written/read.
!!  BUFFER must be of sufficient size to handle the write/read; no bounds
!!  checking is performed.  The intrinsic function STORAGE_SIZE returns the
!!  number of bits needed to store a scalar value of the same type as a given
!!  variable.  This can be used to determine how to size BUFFER (divide by 8
!!  to get bytes.)
!!
!! IMPLEMENTATION NOTES
!!
!!  This is a prototype motivated by a specific immediate need; not much
!!  thought has gone into it.
!!
!!  A combination of C_LOC and C_F_POINTER is used to essentially equivalence
!!  the dummy VAR with an array of bytes.  That array is then copied to the
!!  indicated section of BUFFER (assignment of INT8 kind integers).  The
!!  expectation and assumption is that the bits are simply copied as is with
!!  no interpretation of the type of data. (Note that this would definitely
!!  not be the case if the lhs and rhs were floating point variables.)
!!
!!  What is wanted is just a dumb C memcopy, and that may actually be a better
!!  way of doing it -- use iso_c_binding to access the C library function.
!!

module serialization_tools

  use,intrinsic :: iso_c_binding, only: c_loc, c_f_pointer
  use,intrinsic :: iso_fortran_env, only: int8, int16, int32, int64, real32, real64
  implicit none
  private

  public :: copy_to_bytes, copy_from_bytes

  interface copy_to_bytes
    procedure copy_to_bytes_int8, copy_to_bytes_int16, copy_to_bytes_int32, copy_to_bytes_int64
    procedure copy_to_bytes_real32, copy_to_bytes_real64, copy_to_bytes_real64_1
    procedure copy_to_bytes_log32
  end interface

  interface copy_from_bytes
    procedure copy_from_bytes_int8, copy_from_bytes_int16, copy_from_bytes_int32, copy_from_bytes_int64
    procedure copy_from_bytes_real32, copy_from_bytes_real64, copy_from_bytes_real64_1
    procedure copy_from_bytes_log32
  end interface

contains

  subroutine copy_to_bytes_int8 (var, buffer, offset)
    integer(int8), intent(in), target :: var
    integer, intent(inout) :: offset
    integer(int8), intent(inout) :: buffer(:)
    integer(int8), pointer :: ptr(:)
#ifdef INTEL_BUG20210619
    integer :: n
    n = storage_size(var)/storage_size(buffer)
#else
    integer, parameter :: n = storage_size(var)/storage_size(buffer)
#endif
    call c_f_pointer (c_loc(var), ptr, shape=[n])
    buffer(offset+1:offset+n) = ptr
    offset = offset + n
  end subroutine copy_to_bytes_int8

  subroutine copy_to_bytes_int16 (var, buffer, offset)
    integer(int16), intent(in), target :: var
    integer, intent(inout) :: offset
    integer(int8), intent(inout) :: buffer(:)
    integer(int8), pointer :: ptr(:)
#ifdef INTEL_BUG20210619
    integer :: n
    n = storage_size(var)/storage_size(buffer)
#else
    integer, parameter :: n = storage_size(var)/storage_size(buffer)
#endif
    call c_f_pointer (c_loc(var), ptr, shape=[n])
    buffer(offset+1:offset+n) = ptr
    offset = offset + n
  end subroutine copy_to_bytes_int16

  subroutine copy_to_bytes_int32 (var, buffer, offset)
    integer(int32), intent(in), target :: var
    integer, intent(inout) :: offset
    integer(int8), intent(inout) :: buffer(:)
    integer(int8), pointer :: ptr(:)
#ifdef INTEL_BUG20210619
    integer :: n
    n = storage_size(var)/storage_size(buffer)
#else
    integer, parameter :: n = storage_size(var)/storage_size(buffer)
#endif
    call c_f_pointer (c_loc(var), ptr, shape=[n])
    buffer(offset+1:offset+n) = ptr
    offset = offset + n
  end subroutine copy_to_bytes_int32

  subroutine copy_to_bytes_int64 (var, buffer, offset)
    integer(int64), intent(in), target :: var
    integer, intent(inout) :: offset
    integer(int8), intent(inout) :: buffer(:)
    integer(int8), pointer :: ptr(:)
#ifdef INTEL_BUG20210619
    integer :: n
    n = storage_size(var)/storage_size(buffer)
#else
    integer, parameter :: n = storage_size(var)/storage_size(buffer)
#endif
    call c_f_pointer (c_loc(var), ptr, shape=[n])
    buffer(offset+1:offset+n) = ptr
    offset = offset + n
  end subroutine copy_to_bytes_int64

  subroutine copy_to_bytes_real32 (var, buffer, offset)
    real(real32), intent(in), target :: var
    integer, intent(inout) :: offset
    integer(int8), intent(inout) :: buffer(:)
    integer(int8), pointer :: ptr(:)
#ifdef INTEL_BUG20210619
    integer :: n
    n = storage_size(var)/storage_size(buffer)
#else
    integer, parameter :: n = storage_size(var)/storage_size(buffer)
#endif
    call c_f_pointer (c_loc(var), ptr, shape=[n])
    buffer(offset+1:offset+n) = ptr
    offset = offset + n
  end subroutine copy_to_bytes_real32

  subroutine copy_to_bytes_real64 (var, buffer, offset)
    real(real64), intent(in), target :: var
    integer, intent(inout) :: offset
    integer(int8), intent(inout) :: buffer(:)
    integer(int8), pointer :: ptr(:)
#ifdef INTEL_BUG20210619
    integer :: n
    n = storage_size(var)/storage_size(buffer)
#else
    integer, parameter :: n = storage_size(var)/storage_size(buffer)
#endif
    call c_f_pointer (c_loc(var), ptr, shape=[n])
    buffer(offset+1:offset+n) = ptr
    offset = offset + n
  end subroutine copy_to_bytes_real64

  subroutine copy_to_bytes_real64_1 (var, buffer, offset)
    real(real64), intent(in), contiguous, target :: var(:)
    integer, intent(inout) :: offset
    integer(int8), intent(inout) :: buffer(:)
    integer(int8), pointer :: ptr(:)
    integer :: len
#ifdef INTEL_BUG20210619
    integer :: n
    n = storage_size(var)/storage_size(buffer)
#else
    integer, parameter :: n = storage_size(var)/storage_size(buffer)
#endif
    len = n*size(var)
    call c_f_pointer (c_loc(var(1)), ptr, shape=[len])
    buffer(offset+1:offset+len) = ptr
    offset = offset + len
  end subroutine copy_to_bytes_real64_1

  subroutine copy_to_bytes_log32 (var, buffer, offset)
    logical(int32), intent(in), target :: var
    integer, intent(inout) :: offset
    integer(int8), intent(inout) :: buffer(:)
    integer(int8), pointer :: ptr(:)
#ifdef INTEL_BUG20210619
    integer :: n
    n = storage_size(var)/storage_size(buffer)
#else
    integer, parameter :: n = storage_size(var)/storage_size(buffer)
#endif
    call c_f_pointer (c_loc(var), ptr, shape=[n])
    buffer(offset+1:offset+n) = ptr
    offset = offset + n
  end subroutine copy_to_bytes_log32

  subroutine copy_from_bytes_int8 (buffer, offset, var)
    integer(int8), intent(in) :: buffer(:)
    integer, intent(inout) :: offset
    integer(int8), intent(out), target :: var
    integer(int8), pointer :: ptr(:)
#ifdef INTEL_BUG20210619
    integer :: n
    n = storage_size(var)/storage_size(buffer)
#else
    integer, parameter :: n = storage_size(var)/storage_size(buffer)
#endif
    call c_f_pointer (c_loc(var), ptr, shape=[n])
    ptr = buffer(offset+1:offset+n)
    offset = offset + n
  end subroutine copy_from_bytes_int8

  subroutine copy_from_bytes_int16 (buffer, offset, var)
    integer(int8), intent(in) :: buffer(:)
    integer, intent(inout) :: offset
    integer(int16), intent(out), target :: var
    integer(int8), pointer :: ptr(:)
#ifdef INTEL_BUG20210619
    integer :: n
    n = storage_size(var)/storage_size(buffer)
#else
    integer, parameter :: n = storage_size(var)/storage_size(buffer)
#endif
    call c_f_pointer (c_loc(var), ptr, shape=[n])
    ptr = buffer(offset+1:offset+n)
    offset = offset + n
  end subroutine copy_from_bytes_int16

  subroutine copy_from_bytes_int32 (buffer, offset, var)
    integer(int8), intent(in) :: buffer(:)
    integer, intent(inout) :: offset
    integer(int32), intent(out), target :: var
    integer(int8), pointer :: ptr(:)
#ifdef INTEL_BUG20210619
    integer :: n
    n = storage_size(var)/storage_size(buffer)
#else
    integer, parameter :: n = storage_size(var)/storage_size(buffer)
#endif
    call c_f_pointer (c_loc(var), ptr, shape=[n])
    ptr = buffer(offset+1:offset+n)
    offset = offset + n
  end subroutine copy_from_bytes_int32

  subroutine copy_from_bytes_int64 (buffer, offset, var)
    integer(int8), intent(in) :: buffer(:)
    integer, intent(inout) :: offset
    integer(int64), intent(out), target :: var
    integer(int8), pointer :: ptr(:)
#ifdef INTEL_BUG20210619
    integer :: n
    n = storage_size(var)/storage_size(buffer)
#else
    integer, parameter :: n = storage_size(var)/storage_size(buffer)
#endif
    call c_f_pointer (c_loc(var), ptr, shape=[n])
    ptr = buffer(offset+1:offset+n)
    offset = offset + n
  end subroutine copy_from_bytes_int64

  subroutine copy_from_bytes_real32 (buffer, offset, var)
    integer(int8), intent(in) :: buffer(:)
    integer, intent(inout) :: offset
    real(real32), intent(out), target :: var
    integer(int8), pointer :: ptr(:)
#ifdef INTEL_BUG20210619
    integer :: n
    n = storage_size(var)/storage_size(buffer)
#else
    integer, parameter :: n = storage_size(var)/storage_size(buffer)
#endif
    call c_f_pointer (c_loc(var), ptr, shape=[n])
    ptr = buffer(offset+1:offset+n)
    offset = offset + n
  end subroutine copy_from_bytes_real32

  subroutine copy_from_bytes_real64 (buffer, offset, var)
    integer(int8), intent(in) :: buffer(:)
    integer, intent(inout) :: offset
    real(real64), intent(out), target :: var
    integer(int8), pointer :: ptr(:)
#ifdef INTEL_BUG20210619
    integer :: n
    n = storage_size(var)/storage_size(buffer)
#else
    integer, parameter :: n = storage_size(var)/storage_size(buffer)
#endif
    call c_f_pointer (c_loc(var), ptr, shape=[n])
    ptr = buffer(offset+1:offset+n)
    offset = offset + n
  end subroutine copy_from_bytes_real64

  subroutine copy_from_bytes_real64_1 (buffer, offset, var)
    integer(int8), intent(in) :: buffer(:)
    integer, intent(inout) :: offset
    real(real64), intent(out), contiguous, target :: var(:)
    integer(int8), pointer :: ptr(:)
    integer :: len
#ifdef INTEL_BUG20210619
    integer :: n
    n = storage_size(var)/storage_size(buffer)
#else
    integer, parameter :: n = storage_size(var)/storage_size(buffer)
#endif
    len = n*size(var)
    call c_f_pointer (c_loc(var(1)), ptr, shape=[len])
    ptr = buffer(offset+1:offset+len)
    offset = offset + len
  end subroutine copy_from_bytes_real64_1

  subroutine copy_from_bytes_log32 (buffer, offset, var)
    integer(int8), intent(in) :: buffer(:)
    integer, intent(inout) :: offset
    logical(int32), intent(out), target :: var
    integer(int8), pointer :: ptr(:)
#ifdef INTEL_BUG20210619
    integer :: n
    n = storage_size(var)/storage_size(buffer)
#else
    integer, parameter :: n = storage_size(var)/storage_size(buffer)
#endif
    call c_f_pointer (c_loc(var), ptr, shape=[n])
    ptr = buffer(offset+1:offset+n)
    offset = offset + n
  end subroutine copy_from_bytes_log32

end module serialization_tools
