!!
!! SIMULATION_LOG_TYPE
!!
!! This module defines SIMULATION_LOG, which writes messages for one
!! simulation to a disk log and, optionally, to standard output.  The log is
!! written by the I/O process only.  Standard output is always a subset of
!! the disk log: callers may suppress the standard-output copy of an
!! individual message, but cannot write to standard output alone.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

module simulation_log_type

  use,intrinsic :: iso_fortran_env, only: output_unit
  use parallel_communication, only: is_IOP, broadcast, broadcast_alloc_char
  implicit none
  private

  integer, parameter, public :: LOG_NORMAL = 0
  integer, parameter, public :: LOG_DETAIL = 1

  type, public :: simulation_log
    private
    integer :: log_unit = 0
    integer :: verbosity = LOG_NORMAL
    logical :: terminal_output = .true.
  contains
    procedure :: init
    procedure :: close
    procedure :: info
    procedure :: warn
    procedure :: error
    final :: finalize
  end type simulation_log

contains

  subroutine init(this, filename, stat, errmsg, verbosity, terminal_output)

    class(simulation_log), intent(out) :: this
    character(*), intent(in) :: filename
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, intent(in), optional :: verbosity
    logical, intent(in), optional :: terminal_output

    character(256) :: iomsg

    stat = 0
    if (present(verbosity)) then
      this%verbosity = verbosity
    end if
    if (present(terminal_output)) then
      this%terminal_output = terminal_output
    end if
    if (this%verbosity < LOG_NORMAL .or. this%verbosity > LOG_DETAIL) then
      stat = 1
      errmsg = 'invalid log verbosity'
    else if (len_trim(filename) == 0) then
      stat = 1
      errmsg = 'log filename must not be empty'
    else if (is_IOP) then
      open(newunit=this%log_unit, file=filename, status='replace', action='write', iostat=stat, iomsg=iomsg)
      if (stat /= 0) errmsg = trim(iomsg)
    end if
    call broadcast(stat)
    if (stat /= 0) then
      call broadcast_alloc_char(errmsg)
      return
    end if

  end subroutine init


  subroutine close(this)
    class(simulation_log), intent(inout) :: this
    if (is_IOP .and. this%log_unit /= 0) close(this%log_unit)
    this%log_unit = 0
  end subroutine

  subroutine info(this, message, level, terminal)
    class(simulation_log), intent(in) :: this
    character(*), intent(in) :: message
    integer, intent(in), optional :: level
    logical, intent(in), optional :: terminal
    integer :: message_level
    logical :: write_terminal
    message_level = LOG_NORMAL
    if (present(level)) message_level = level
    if (message_level > this%verbosity .or. .not.is_IOP) return
    write(this%log_unit, '(a)') trim(message)
    write_terminal = .true.
    if (present(terminal)) write_terminal = terminal
    if (this%terminal_output .and. write_terminal) write(output_unit,'(a)') trim(message)
  end subroutine

  subroutine warn(this, message, terminal)
    class(simulation_log), intent(in) :: this
    character(*), intent(in) :: message
    logical, intent(in), optional :: terminal
    logical :: write_terminal
    if (.not.is_IOP) return
    write(this%log_unit, '(2a)') 'Warning: ', trim(message)
    write_terminal = .true.
    if (present(terminal)) write_terminal = terminal
    if (this%terminal_output .and. write_terminal) write(output_unit,'(2a)') 'Warning: ', trim(message)
  end subroutine

  subroutine error(this, message, terminal)
    class(simulation_log), intent(in) :: this
    character(*), intent(in) :: message
    logical, intent(in), optional :: terminal
    logical :: write_terminal
    if (.not.is_IOP) return
    write(this%log_unit, '(2a)') 'ERROR: ', trim(message)
    write_terminal = .true.
    if (present(terminal)) write_terminal = terminal
    if (this%terminal_output .and. write_terminal) write(output_unit,'(2a)') 'ERROR: ', trim(message)
  end subroutine

  subroutine finalize(this)
    type(simulation_log), intent(inout) :: this
    call this%close()
  end subroutine

end module simulation_log_type
