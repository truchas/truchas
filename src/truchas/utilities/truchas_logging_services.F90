!!
!! TRUCHAS_LOGGING_SERVICES
!!
!! This module provides a common set of procedures for writing log/trace
!! messages -- the types of messages typically written to the terminal or
!! simulation log file.  All Truchas code is expected to use the facilities
!! provided here, and eschew use of raw Fortran writes for this purpose.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! April 2012
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  There are three categories of log messages. Informational messages are the
!!  typical generic messages, written as-is without any adornment.  These are
!!  logged using TLS_INFO.  Warning messages are intended to warn the user of
!!  a possible non-fatal error or exceptional condition.  These are written
!!  prefixed with "Warning:" to be readily distinguished from other messages.
!!  These are logged using TLS_WARN.  The final category is fatal error
!!  messages, and these are logged using either TLS_FATAL or TLS_ERROR.  The
!!  former prefixes the message with "FATAL:" and gracefully terminates
!!  execution.  The latter prefixes the message with "ERROR:" but does not
!!  terminate execution.  However in this case it is required that execution
!!  ultimately be terminated using TLS_FATAL.  The intent is to permit the
!!  logging of multiple fatal errors before finally terminating execution.
!!
!!  Informational messages are subdivided into several verbosity levels.
!!  Only those messages not exceeding the system verbosity level are logged.
!!  The verbosity levels are provided by the named module constants
!!  TLS_VERB_SILENT, TLS_VERB_NORMAL, and TLS_VERB_NOISY.  The module variable
!!  TLS_VERBOSITY (read-only!) gives the system verbosity level.  These values
!!  are of a private derived-type, and thus are only useful for passing to the
!!  procedures described below (the VERBOSITY argument), or in a comparison
!!  expression like "if (TLS_VERBOSITY >= TLS_VERB_NOISY) ..."  All of the
!!  comparison operators are defined for the type.
!!
!!  The procedures TLS_INFO, TLS_WARN, TLS_ERROR, TLS_FATAL should be regarded
!!  as collective procedures, needing to be invoked by every parallel process.
!!  In reality only TLS_FATAL is truly collective (to terminate execution.)
!!  Logging is done by the designated I/O process only, using the message on
!!  that process; other processes do nothing, except when terminating execution
!!  with TLS_FATAL.
!!
!!  Because logging a fatal error message is a collective operation, it is
!!  often required to communicate the presence of an error condition on one
!!  process across the entire parallel machine so that it can be properly
!!  logged.  The convenience procedures TLS_FATAL_IF_ANY, TLS_FATAL_IF_IOP
!!  which wrap TLS_FATAL, handles this communication.
!!
!!  Finally, there are occasions when it is not possible or convenient to
!!  gracefully terminate execution in a collective manner.  The procedure
!!  TLS_PANIC exists for this purpose.  It can be called by one or more
!!  processes, and each writes its message prefixed with "PANIC[N]:" to the
!!  terminal (only), where N is the process rank, and then execution is
!!  terminated via MPI_ABORT, which is supposed to terminate all processes,
!!  including those not calling TLS_PANIC.  There may be issues with losing
!!  some output.  Use of this procedure should be kept to a minimum.
!!
!!  As currently implemented, messages are logged to both the terminal (C's
!!  stdout or Fortran's output_unit) and to the *.log file in the run's
!!  output directory.  TLS_PANIC messages are written only to the terminal.
!!
!!  What follows is a more detailed description of the procedures.  Note
!!  that some variants allow an array MESSAGE which is handled as a multiline
!!  message.  This is mostly for compatability with existing Truchas usage;
!!  new use of multiline messages is discouraged.
!!
!!  CALL TLS_INFO (MESSAGE[, VERBOSITY]) writes the message passed in the
!!    character argument MESSAGE, which may be either a scalar or a rank-1
!!    array.  In the latter case each array element is written starting on
!!    a new line.  VERBOSITY specifies the verbosity level of the message,
!!    and defaults to TLS_VERB_NORMAL.  If the message verbosity level exceeds
!!    the system verbosity level (TLS_VERBOSITY) the message is not written.
!!    For clarity this should be called collectively, but only the designated
!!    I/O process logs its message.
!!
!!  CALL TLS_INFO (MESSAGE, ADVANCE[, VERBOSITY]) writes the message passed
!!    in the scalar character argument MESSAGE.  If the logical argument
!!    ADVANCE is false, the usual trailing newline is not written, otherwise
!!    the behavior is exactly as for the first variant of TLS_INFO.  Note
!!    that this variant is not defined for array MESSAGE.  VERBOSITY specifies
!!    the verbosity level of the message and defaults to TLS_VERB_NORMAL.
!!    If the message verbosity level exceeds the system verbosity level
!!    (TLS_VERBOSITY) the message is not written.
!!    For clarity this should be called collectively, but only the designated
!!    I/O process logs its message.
!!
!!  CALL TLS_WARN (MESSAGE) writes the warning message passed in the character
!!    argument MESSAGE, which may be either scalar or a rank-1 array.  The
!!    written message is prefixed with the string "Warning: ".  When MESSAGE
!!    is an array, the second and subsequent lines are indented so that the
!!    message appears as a block to the right of the "Warning: " prefix.
!!    For clarity this should be called collectively, but only the designated
!!    I/O process logs its message.
!!
!!  CALL TLS_ERROR (MESSAGE) writes the error message passed in the character
!!    argument MESSAGE, which may be either scalar or a rank-1 array.  The
!!    written message is prefixed with the string "ERROR: ".  When MESSAGE
!!    is an array, the second and subsequent lines are indented so that the
!!    message appears as a block to the right of the "ERROR: " prefix.  This
!!    procedure does not initiate termination of the simulation, however,
!!    the caller is expected to ultimately do so using TLS_FATAL.
!!    For clarity this should be called collectively, but only the designated
!!    I/O process logs its message.
!!
!!  CALL TLS_FATAL (MESSAGE) writes the error message passed in the scalar
!!    character argument MESSAGE.  The written message is prefixed with the
!!    string "FATAL: ".  After logging the message the parallel machine is
!!    brought down (MPI_FINALIZE) and execution gracefully terminated.  This
!!    is a collective procedure.
!!
!!  CALL TLS_FATAL_IF_ANY (ERRC, MESSAGE) calls TLS_FATAL with the passed
!!    message if the logical argument ERRC has value true on any of the
!!    processes.  This is a collective procedure.
!!
!!  CALL TLS_FATAL_IF_IOP (ERRC, MESSAGE) calls TLS_FATAL with the passed
!!    message if the logical argument ERRC has value true on the IO process;
!!    ERRC is ignored on other processes.  This is a collective procedure.
!!
!!  CALL TLS_PANIC (MESSAGE) writes the message passed in the scalar character
!!    argument MESSAGE to the terminal, prefixed with the string "PANIC[X]: ",
!!    where X is the rank of the process.  This is not a collective procedure;
!!    any number of processes may call it and each process writes its own
!!    message.  After writing the message execution is ungracefully terminated
!!    via MPI_ABORT.
!!
!!  CALL TLS_EXIT logs the fixed normal termination message (no passed message)
!!    and the parallel machine is brought down (MPI_FINALIZE) and execution
!!    gracefully terminated.  This is a collective procedure.
!!
!!  Packages that handle their own logging (e.g. Ubik) often accept a unit
!!  number to use for that purpose.  The function TLS_LOGGING_UNIT() returns
!!  a unit number used by this module for logging (on the I/O process, -1 on
!!  all other processes.  Currently this is the unit for the .log file.
!!  [This use case argues for logging to a single device, and not several.]
!!
!!  The following procedures initialize and finalize the logging system.  They
!!  are called by the Truchas driver, and should never be invoked by normal
!!  application code.
!!
!!    CALL TLS_INITIALIZE
!!    CALL TLS_SET_VERBOSITY (VERBOSITY) sets the system verbosity level.
!!    CALL TLS_FINALIZE
!!
!!  The following procedure is provided for use in ad hoc debugging only.
!!  Use when needed, but then remove it -- DO NOT LEAVE IT IN THE CODE.
!!
!!  CALL TLS_DEBUG (MESSAGE) writes the input scalar character MESSAGE to a
!!    processor-specific file.  The file is named debug.<N>, where <N> is the
!!    process rank, and it is created in the directory from which the
!!    executable was run.  After each write the output buffer is flushed.
!!    The file is created only if needed.  This is NOT a collective procedure.
!!
!!  CALL TLS_GET_DEBUG_UNIT (UNIT) returns the logical unit that the
!!    processor-specific debug file is opened on.  If the file is not already
!!    open, this call opens the file as well.  This is NOT a collective procedure.
!!
!! IMPLEMENTATION NOTES
!!
!! * I've continued the tradition of logging to both the terminal (stdout) and
!!   the .log file. We should revisit this choice of two output devices at some
!!   point. It would be simpler to limit to a single device, terminal or disk
!!   file.
!!
!! * I would like the routines that terminate execution to return a status
!!   code to the OS (via the system exit command, for example).  This would
!!   make it trivial for shell scripts to know whether a run was successful
!!   or not.  This isn't fully possible for the following reasons:
!!   = Serial.  PGSLIB_Abort invokes STOP.  PGSLIB_Finalize is a no-op and
!!     returns control, so there at least exit could be called.
!!   = Open MPI (1.5).  MPI_Finalize returns control and exit can be called,
!!     but it is mpirun that receives the status code, not the OS, and mpirun
!!     does not itself return the status code.  All it does is write useless
!!     info to stderr.  The situation is effectively the same for MPI_Abort,
!!     though it's not clear that it returns control to the caller.
!!   As a consequence, I'm not bothering with calling exit at this time.
!!
!! POSSIBLE TO-DOS
!! * Create debug files in output file directory along with everything else.
!! * option to suppress warning messages
!! * pass the is_IOP flag into initialization and eliminate dependence on
!!   parallel_communication module.
!! * flushing of the output buffers.
!!

#include "f90_assert.fpp"

module truchas_logging_services

  use parallel_communication, only: is_IOP
  implicit none
  private

  public :: TLS_initialize, TLS_finalize, TLS_set_verbosity, TLS_logging_unit
  public :: TLS_info, TLS_warn, TLS_error, TLS_fatal, TLS_exit
  public :: TLS_fatal_if_any, TLS_fatal_if_IOP
  public :: TLS_panic, TLS_debug, TLS_debug_unit

  interface TLS_info
    module procedure TLS_info_scalar, TLS_info_array, TLS_info_advance
  end interface TLS_info

  interface TLS_warn
    module procedure TLS_warn_scalar, TLS_warn_array
  end interface TLS_warn

  interface TLS_error
    module procedure TLS_error_scalar, TLS_error_array
  end interface

  !! Log messages are written to these units.
  integer, allocatable, save :: log_unit(:)

  !! Debug messages are written to this unit.
  integer, save :: dbg_unit = -1

  !! Private type to describe verbosity levels.
  type :: verb_level
    private
    integer :: level
  end type verb_level

  !! Named verbosity level constants; normal and noisy will do for now.
  type(verb_level), parameter, public :: TLS_VERB_SILENT = verb_level(0)
  type(verb_level), parameter, public :: TLS_VERB_NORMAL = verb_level(1)
  type(verb_level), parameter, public :: TLS_VERB_NOISY  = verb_level(2)

  !! The system verbosity level; message levels are compared against this.
  type(verb_level), save, public :: TLS_verbosity = TLS_VERB_NORMAL ! PROTECTED!

  interface operator(==)
    module procedure verb_level_eq
  end interface

  interface operator(/=)
    module procedure verb_level_ne
  end interface

  interface operator(<)
    module procedure verb_level_lt
  end interface

  interface operator(<=)
    module procedure verb_level_le
  end interface

  interface operator(>)
    module procedure verb_level_gt
  end interface

  interface operator(>=)
    module procedure verb_level_ge
  end interface

  public :: operator(==), operator(/=), operator(<), operator(<=), operator(>), operator(>=)

contains

  subroutine TLS_initialize
    use truchas_env, only: output_file_name
    use,intrinsic :: iso_fortran_env, only: output_unit
    !! Initialize the devices the log messages will be written to.
    if (is_IOP) then
      allocate(log_unit(2))
      log_unit(1) = output_unit ! pre-connected output unit (stdout)
      open(newunit=log_unit(2), file=output_file_name('log'))
    end if
  end subroutine TLS_initialize

  subroutine TLS_set_verbosity (verbosity)
    type(verb_level), intent(in) :: verbosity
    TLS_verbosity = verbosity
  end subroutine TLS_set_verbosity

  integer function TLS_logging_unit ()
    if (is_IOP .and. allocated(log_unit)) then
      TLS_logging_unit = log_unit(2)  ! we have to choose one
    else
      TLS_logging_unit = -1
    end if
  end function TLS_logging_unit

  subroutine TLS_finalize
    !! If/when TLS_initialize opens a file, it should be closed here.
    if (is_IOP .and. allocated(log_unit)) then
      close(log_unit(2))
      deallocate(log_unit)
    end if
    if (dbg_unit /= -1) then
      close(dbg_unit)
      dbg_unit = -1
    end if
  end subroutine TLS_finalize

 !!
 !! SPECIFIC TLS_INFO PROCEDURES
 !!

  subroutine TLS_info_scalar (message, verbosity)
    character(*), intent(in) :: message
    type(verb_level), intent(in), optional :: verbosity
    integer :: n
    type(verb_level) :: level
    if (is_IOP) then
      if (present(verbosity)) then
        ASSERT(verbosity > TLS_VERB_SILENT)
        level = verbosity
      else
        level = TLS_VERB_NORMAL
      end if
      if (level <= TLS_verbosity) then
        do n = 1, size(log_unit)
          write(log_unit(n),'(a)') message(:len_trim(message))
        end do
      end if
    end if
  end subroutine TLS_info_scalar

  subroutine TLS_info_array (message, verbosity)
    character(*), intent(in) :: message(:)
    type(verb_level), intent(in), optional :: verbosity
    integer :: n, m
    type(verb_level) :: level
    if (is_IOP) then
      if (present(verbosity)) then
        ASSERT(verbosity > TLS_VERB_SILENT)
        level = verbosity
      else
        level = TLS_VERB_NORMAL
      end if
      if (level <= TLS_verbosity) then
        do m = 1, size(message)
          do n = 1, size(log_unit)
            write(log_unit(n),'(a)') message(m)(:len_trim(message(m)))
          end do
        end do
      end if
    end if
  end subroutine TLS_info_array

  subroutine TLS_info_advance (message, advance, verbosity)
    character(*), intent(in) :: message
    logical, intent(in) :: advance
    type(verb_level), intent(in), optional :: verbosity
    integer :: n
    type(verb_level) :: level
    if (is_IOP) then
      if (present(verbosity)) then
        ASSERT(verbosity > TLS_VERB_SILENT)
        level = verbosity
      else
        level = TLS_VERB_NORMAL
      end if
      if (level <= TLS_verbosity) then
        if (advance) then
          do n = 1, size(log_unit)
            write(log_unit(n),'(a)') message(:len_trim(message))
          end do
        else
          do n = 1, size(log_unit)
            write(log_unit(n),'(a)',advance='no') message ! intentionally not trimmed
          end do
        end if
      end if
    end if
  end subroutine TLS_info_advance

 !!
 !! SPECIFIC TLS_WARN AND TLS_ERROR PROCEDURES
 !!

  subroutine TLS_warn_scalar (message)
    character(*), intent(in) :: message
    call labeled_message_scalar ('Warning: ', message)
  end subroutine TLS_warn_scalar

  subroutine TLS_warn_array (message)
    character(*), intent(in) :: message(:)
    call labeled_message_array ('Warning: ', message)
  end subroutine TLS_warn_array

  subroutine TLS_error_scalar (message)
    character(*), intent(in) :: message
    call labeled_message_scalar ('ERROR: ', message)
  end subroutine TLS_error_scalar

  subroutine TLS_error_array (message)
    character(*), intent(in) :: message(:)
    call labeled_message_array ('ERROR: ', message)
  end subroutine TLS_error_array

  subroutine labeled_message_scalar (label, message)
    character(*), intent(in) :: label, message
    integer :: n
    if (is_IOP) then
      do n = 1, size(log_unit)
        write(log_unit(n),'(2a)') label, message(:len_trim(message))
      end do
    end if
  end subroutine labeled_message_scalar

  subroutine labeled_message_array (label, message)
    character(*), intent(in) :: label, message(:)
    integer :: n, m
    ASSERT(size(message) >= 1)
    if (is_IOP) then
      do n = 1, size(log_unit)
        write(log_unit(n),'(2a)') label, message(1)(:len_trim(message(1)))
      end do
      do m = 2, size(message)
        do n = 1, size(log_unit)
          write(log_unit(n),'(2a)') repeat(' ',len(label)), message(m)(:len_trim(message(m)))
        end do
      end do
    end if
  end subroutine labeled_message_array

  subroutine TLS_fatal (message)
#ifdef NAG_COMPILER
    use,intrinsic :: f90_unix, only: exit
#endif
    use truchas_danu_output_data, only: outfile
    use pgslib_module, only: pgslib_finalize
    use utilities_module, only: timestamp
    character(*), intent(in) :: message
    character(32) :: date_time
    call labeled_message_scalar ('FATAL: ', message)
    call outfile%close ()
    call timestamp (date_time)
    call TLS_info ('truchas terminated abnormally on '//date_time(5:13)//' at '//date_time(15:22))
    call TLS_finalize
    call pgslib_finalize
    call exit (1)
  end subroutine TLS_fatal

  subroutine TLS_fatal_if_any (errc, message)
    use parallel_communication, only: global_any
    logical, intent(in) :: errc
    character(*), intent(in) :: message
    if (global_any(errc)) call TLS_fatal (message)
  end subroutine TLS_fatal_if_any

  subroutine TLS_fatal_if_IOP (errc, message)
    use parallel_communication, only: is_IOP, broadcast
    logical, intent(in) :: errc
    character(*), intent(in) :: message
    logical :: global_errc
    if (is_IOP) global_errc = errc
    call broadcast (global_errc)
    if (global_errc) call TLS_fatal (message)
  end subroutine TLS_fatal_if_IOP

  subroutine TLS_exit
#ifdef NAG_COMPILER
    use,intrinsic :: f90_unix, only: exit
#endif
    use pgslib_module, only: pgslib_finalize
    use utilities_module, only: timestamp
    character(32) :: date_time
    call timestamp (date_time)
    call TLS_info ('')
    call TLS_info ('truchas terminated normally on '//date_time(5:13)//' at '//date_time(15:22))
    call TLS_finalize
    call pgslib_finalize
    call exit (0)
  end subroutine TLS_exit

 !!
 !! TLS_PANIC
 !!
 !! Write message and ungracefully abort execution.  This is not collective and
 !! intended only for situations where a graceful collective termination using
 !! TLS_FATAL is not possible or convenient.  Because only some processes may
 !! call this, the usual log file (known only to the designated I/O process) is
 !! not available.  Thus each process writes its own message (they need not be
 !! the same) to the terminal, prefixed with the process rank, and we rely on
 !! the OS to do the collation.  Termination of *all* processes is accomplished
 !! by calling (indirectly) MPI_ABORT.  This is somewhat nasty and can result in
 !! some lost ouput, and perhaps even hung processes.  This procedure should not
 !! be used lightly or routinely.
 !!

  subroutine TLS_panic (message)
#ifdef NAG_COMPILER
    use,intrinsic :: f90_unix, only: exit
#endif
    use,intrinsic :: iso_fortran_env, only: output_unit
    use parallel_communication, only: this_PE
    use pgslib_module, only: pgslib_abort
    use utilities_module, only: timestamp
    character(*), intent(in) :: message
    character(32) :: date_time
    write(output_unit,'(a,i0,a)') 'PANIC[', this_PE, ']: ' // message(:len_trim(message))
    call timestamp (date_time)
    write(output_unit,'(a)') 'Truchas terminated abnormally on '//date_time(5:13)//' at '//date_time(15:22)
    call pgslib_abort
    call exit (1)
  end subroutine TLS_panic

 !!
 !! TLS_DEBUG -- Write debug message to processor-specific file
 !!

  subroutine TLS_debug (message)
    character(*), intent(in) :: message
    if (dbg_unit == -1) call open_debug_file
    write(dbg_unit,'(a)') message
  end subroutine TLS_debug

  subroutine TLS_get_debug_unit (unit)
    integer, intent(out) :: unit
    if (dbg_unit == -1) call open_debug_file
    unit = dbg_unit
  end subroutine TLS_get_debug_unit

  subroutine open_debug_file
    use parallel_communication, only: nPE, this_PE
    use string_utilities, only: i_to_c
    character(16) :: string
    character(128) :: dbg_file
    integer :: n, ios
    if (dbg_unit == -1) then  ! open the debug file
      write(string,'(i0)') nPE
      n = len_trim(string)
      write(string,'("(a,i",i0,".",i0,")")') n, n
      write(dbg_file,string) 'debug.', this_PE
      open(newunit=dbg_unit,file=trim(dbg_file),action='write',status='replace',iostat=ios)
      if (ios /= 0) call TLS_panic ('TLS_DEBUG: error opening file "' // trim(dbg_file) // &
                                    '": iostat=' // i_to_c(ios))
    end if
  end subroutine open_debug_file
    

  integer function TLS_debug_unit () result (unit)
    use parallel_communication, only: nPE, this_PE
    use string_utilities, only: i_to_c
    character(16) :: string
    character(128) :: dbg_file
    integer :: n, ios
    if (dbg_unit == -1) then  ! first time -- open the debug file
      write(string,'(i0)') nPE
      n = len_trim(string)
      write(string,'("(a,i",i0,".",i0,")")') n, n
      write(dbg_file,string) 'debug.', this_PE
      open(newunit=dbg_unit,file=trim(dbg_file),action='write',status='replace',iostat=ios)
      if (ios /= 0) call TLS_panic ('TLS_DEBUG: error opening file "' // trim(dbg_file) // &
                                    '": iostat=' // i_to_c(ios))
    end if
    unit = dbg_unit
  end function TLS_debug_unit

 !!
 !! COMPARISON OPERATORS FOR TYPE(VERB_LEVEL) OBJECTS
 !!

  pure logical function verb_level_eq (a, b)
    type(verb_level), intent(in) :: a, b
    verb_level_eq = (a%level == b%level)
  end function verb_level_eq

  pure logical function verb_level_ne (a, b)
    type(verb_level), intent(in) :: a, b
    verb_level_ne = (a%level /= b%level)
  end function verb_level_ne

  pure logical function verb_level_lt (a, b)
    type(verb_level), intent(in) :: a, b
    verb_level_lt = (a%level < b%level)
  end function verb_level_lt

  pure logical function verb_level_le (a, b)
    type(verb_level), intent(in) :: a, b
    verb_level_le = (a%level <= b%level)
  end function verb_level_le

  pure logical function verb_level_gt (a, b)
    type(verb_level), intent(in) :: a, b
    verb_level_gt = (a%level > b%level)
  end function verb_level_gt

  pure logical function verb_level_ge (a, b)
    type(verb_level), intent(in) :: a, b
    verb_level_ge = (a%level >= b%level)
  end function verb_level_ge

end module truchas_logging_services
