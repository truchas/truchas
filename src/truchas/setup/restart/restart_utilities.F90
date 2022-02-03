!!
!! RESTART_UTILITIES
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! 26 Apr 2005, revised 13 May 2005
!!
!! This module provides a collection of utility procedures that simplify and
!! streamline the process of reading the restart file, and it should be used
!! only in that context.  If something here proves useful in a broader
!! setting, please move it to a different module---don't use this module.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PUBLIC PROCEDURES
!!
!!  First some general comments about the procedures:
!!
!!    o These are designed as collective procedures (termed 'global'), and
!!      thus must be called from all processors with arguments of the same
!!      type, kind, and rank.
!!
!!    o Some of the procedures have two forms distinguished by whether the
!!      argument STAT or ERRMSG is specified.  In the STAT form, the non-zero
!!      iostat value from a failed I/O operation is returned in the intent(out)
!!      integer STAT; otherwise 0 is returned.  As global procedures, the same
!!      STAT value is returned on all processors.  In the ERRMSG form, an error
!!      message is written using the intent(in) string ERRMSG and execution
!!      halted if an I/O error occurs.
!!
!!    o The subroutines with the UNIT argument assume that a file has been
!!      opened for unformatted reading on UNIT.
!!
!!  CALL INFO (MESSAGE) writes the string MESSAGE to the out and tty units.
!!
!!  CALL HALT (ERRMSG) writes the string MESSAGE prepended with 'FATAL: ' to
!!    the out and tty units, finalizes all I/O and parallel communication,
!!    and then halts execution of the program.  This is more or less
!!    equivalent to PUNT, but less hideous.
!!
!!  CALL MAKE_UNDEFINED (VALUE) formally causes the argument VALUE to become
!!    undefined.  With certain compilers and options (e.g., the NAG compiler
!!    with '-C') the argument actually becomes undefined so that a runtime
!!    error results if its value is ever used without being defined.  This
!!    is useful for catching programming errors.  This subroutine is
!!    elemental: the argument can be any intrinsic numeric type, scalar or
!!    array of any rank.
!!
!!  CALL SKIP_RECORDS (UNIT, N, {STAT|ERRMSG}) skips over N records in the
!!    file opened on UNIT.
!!
!!  CALL READ_VAR (UNIT, VAR, {STAT|ERRMSG}) reads VAR from the current record
!!    of the file opened on UNIT.  The intent(out) argument VAR may be a scalar
!!    or rank-1 array of integer, real, or character type.  The data read must
!!    belong to a single record.  Although the input occurs only on the I/O
!!    processor, the data is broadcast to all other processors so that VAR
!!    returns the same value on all processors.
!!
!!  CALL READ_DIST_ARRAY (UNIT, ARRAY, [PERM,] {STAT|ERRMSG}) reads the
!!    distributed ARRAY from the file opened on UNIT.  ARRAY is an intent(out)
!!    array of integer or real type.  It is interpreted as a distributed array
!!    with the final index being the distributed index.  The extent for each
!!    dimension except the last must be the same across all processors.
!!    If ARRAY is rank-1, the I/O processor reads the data for all arrays from
!!    a single record, and then distributes to each processor its piece, which
!!    is returned in ARRAY.  The optional intent(in) argument PERM is a
!!    rank-1 integer array that specifies a permutation (internal-to-external)
!!    of the distributed index.  If PERM is present, the data is reordered
!!    accordingly.  PERM may be a distributed array, distributed as ARRAY;
!!    otherwise PERM on the I/O processor is used for the permutation, and
!!    PERM should have size 0 on all other processors.  When ARRAY is multi-
!!    dimensional, the preceding operations are done for each rank-1 slice
!!    through the last dimension, in rray element order.  For example,
!!    ARRAY(1,1,:) is read from the current record, ARRAY(2,1,:) is read from
!!    the next record, and so on.  Currently, ARRAY may be rank-1, 2, or 3.
!!

#include "f90_assert.fpp"

module restart_utilities

  use kinds, only: r8
  use,intrinsic :: iso_fortran_env, only: int8

  implicit none
  private

  public :: read_var, read_dist_array, skip_records, make_undefined, info, halt

  interface read_var
    module procedure read_var_stat_int8_0, read_var_halt_int8_0
    module procedure read_var_stat_int8_1, read_var_halt_int8_1
    module procedure read_var_stat_I0, read_var_halt_I0
    module procedure read_var_stat_I1, read_var_halt_I1
    module procedure read_var_stat_R0, read_var_halt_R0
    module procedure read_var_stat_R1, read_var_halt_R1
    module procedure read_var_stat_S0, read_var_halt_S0
    module procedure read_var_stat_S1, read_var_halt_S1
  end interface

  interface read_dist_array
    module procedure read_array_stat_int8_1, read_array_halt_int8_1
    module procedure read_array_stat_int8_2, read_array_halt_int8_2
    module procedure read_array_stat_int8_3, read_array_halt_int8_3
    module procedure read_array_stat_I1, read_array_halt_I1
    module procedure read_array_stat_I2, read_array_halt_I2
    module procedure read_array_stat_I3, read_array_halt_I3
    module procedure read_array_stat_R1, read_array_halt_R1
    module procedure read_array_stat_R2, read_array_halt_R2
    module procedure read_array_stat_R3, read_array_halt_R3
  end interface

  interface make_undefined
    module procedure make_undefined_integer, make_undefined_single, make_undefined_double
  end interface

  interface skip_records
    module procedure skip_records_stat, skip_records_halt
  end interface

contains

  subroutine info (message)
    use truchas_logging_services
    character(len=*), intent(in) :: message
    call TLS_info (message)
  end subroutine info

  subroutine halt (errmsg)
    use truchas_logging_services
    character(len=*), intent(in) :: errmsg
    call TLS_fatal (errmsg)
  end subroutine halt

  !!
  !! Specific procedures for the generic subroutine MAKE_UNDEFINED.
  !!

  elemental subroutine make_undefined_integer (value)
    integer, intent(out) :: value
  end subroutine make_undefined_integer

  elemental subroutine make_undefined_single (value)
    real, intent(out) :: value
  end subroutine make_undefined_single

  elemental subroutine make_undefined_double (value)
    real(kind(1.0d0)), intent(out) :: value
  end subroutine make_undefined_double

  !!
  !! Specific procedures for the generic subroutine SKIP_RECORDS.
  !!

  subroutine skip_records_stat (unit, n, stat)
    use parallel_communication, only: is_IOP
    use pgslib_module, only:  pgslib_bcast
    integer, intent(in)  :: unit, n
    integer, intent(out) :: stat
    integer :: j
    if (is_IOP) then
      do j = 1, n
        read(unit,iostat=stat)
        if (stat /= 0) exit
      end do
    end if
    call pgslib_bcast (stat)
  end subroutine skip_records_stat

  subroutine skip_records_halt (unit, n, errmsg)
    integer, intent(in) :: unit, n
    character(len=*), intent(in) :: errmsg
    integer :: stat
    call skip_records_stat (unit, n, stat)
    if (stat /= 0) call halt (errmsg)
  end subroutine skip_records_halt

  !!
  !! Specific procedures for the generic subroutines READ_VAR and READ_DIST_ARRAY.
  !!

#define _INT8_DATA_
#include "restart_utilities.fpp"

#define _INTEGER_DATA_
#include "restart_utilities.fpp"

#define _REAL_DATA_
#include "restart_utilities.fpp"

#define _STRING_DATA_
#include "restart_utilities.fpp"

end module restart_utilities
