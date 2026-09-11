!!
!! SIMULATION_LOG_TYPE
!!
!! This module defines SIMULATION_LOG, which writes messages for one
!! simulation to a disk log and, optionally, to standard output.  The log is
!! written by the I/O process only. Standard output is always a subset of the
!! disk log: callers may suppress the terminal copy of an individual message
!! or throttle a group of messages, but cannot write to the terminal alone.
!!
!! Terminal throttling is controlled explicitly with begin/end commands. One
!! visibility decision applies to normal INFO messages in the group; all
!! records are still written to the disk log. These commands do not affect
!! section nesting or indentation.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "t2d_assert.inc"

module simulation_log_type

  use,intrinsic :: iso_fortran_env, only: output_unit, real64
  use mpi_f08
  implicit none
  private

  integer, parameter, public :: LOG_NORMAL = 0
  integer, parameter, public :: LOG_DETAIL = 1
  real(real64), parameter :: DEFAULT_TERMINAL_STEP_INTERVAL = 3.0_real64

  type, public :: simulation_log
    private
    integer :: log_unit = 0
    integer :: verbosity = LOG_NORMAL
    integer :: indentation = 0
    real(real64) :: terminal_step_interval = DEFAULT_TERMINAL_STEP_INTERVAL
    real(real64) :: last_terminal_step_time = 0.0_real64
    logical :: terminal_output = .true.
    logical :: io_process = .false.
    logical :: terminal_throttle_active = .false.
    logical :: terminal_group_visible = .true.
    logical :: have_terminal_step_time = .false.
  contains
    procedure :: init
    procedure :: close
    procedure :: unit
    procedure :: terminal_output_enabled
    procedure :: is_enabled
    procedure :: info
    procedure :: begin_section
    procedure :: end_section
    procedure :: begin_terminal_throttle_group
    procedure :: end_terminal_throttle_group
    procedure :: warn
    procedure :: error
    final :: finalize
  end type simulation_log

contains

  subroutine init(this, comm, filename, stat, errmsg, verbosity, terminal_output, terminal_step_interval)

    use,intrinsic :: ieee_arithmetic, only: ieee_is_finite

    class(simulation_log), intent(out) :: this
    type(MPI_Comm), intent(in) :: comm
    character(*), intent(in) :: filename
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, intent(in), optional :: verbosity
    logical, intent(in), optional :: terminal_output
    !! Global terminal switch (default true) and throttle-group interval
    !! in seconds (default three; zero disables throttling).
    real(real64), intent(in), optional :: terminal_step_interval

    character(256) :: iomsg
    integer :: ierr, rank

    call MPI_Comm_rank(comm, rank, ierr)
    INSIST(ierr == MPI_SUCCESS)

    this%io_process = rank == 0
    if (present(verbosity)) this%verbosity = verbosity
    if (present(terminal_output)) this%terminal_output = terminal_output
    if (present(terminal_step_interval)) this%terminal_step_interval = terminal_step_interval
    stat = 0
    if (this%io_process) then
      if (this%verbosity < LOG_NORMAL .or. this%verbosity > LOG_DETAIL) then
        stat = 1
        errmsg = 'invalid log verbosity'
      else if (.not.ieee_is_finite(this%terminal_step_interval) .or. &
          this%terminal_step_interval < 0.0_real64) then
        stat = 1
        errmsg = 'terminal step interval must be nonnegative'
      else if (len_trim(filename) == 0) then
        stat = 1
        errmsg = 'log filename must not be empty'
      else
        open(newunit=this%log_unit, file=filename, status='replace', action='write', iostat=stat, iomsg=iomsg)
        if (stat /= 0) errmsg = trim(iomsg)
      end if
    end if
    call MPI_Bcast(stat, 1, MPI_INTEGER, 0, comm, ierr)
    INSIST(ierr == MPI_SUCCESS)
    if (stat /= 0) then
      call broadcast_errmsg(comm, this%io_process, errmsg)
      return
    end if

  end subroutine init

  subroutine close(this)
    class(simulation_log), intent(inout) :: this
    if (this%io_process .and. this%log_unit /= 0) close(this%log_unit)
    this%log_unit = 0
  end subroutine

  integer function unit(this)
    class(simulation_log), intent(in) :: this
    INSIST(this%io_process .and. this%log_unit /= 0)
    unit = this%log_unit
  end function

  logical function terminal_output_enabled(this)
    class(simulation_log), intent(in) :: this

    terminal_output_enabled = this%io_process .and. this%terminal_output
  end function

  logical function is_enabled(this, level)
    class(simulation_log), intent(in) :: this
    integer, intent(in) :: level
    is_enabled = level <= this%verbosity
  end function

  subroutine info(this, message, level, terminal)

    class(simulation_log), intent(in) :: this
    character(*), intent(in) :: message
    integer, intent(in), optional :: level
    logical, intent(in), optional :: terminal

    call write_info(this, message, level, terminal)

  end subroutine


  subroutine begin_section(this, message, level, terminal)

    class(simulation_log), intent(inout) :: this
    character(*), intent(in) :: message
    integer, intent(in), optional :: level
    logical, intent(in), optional :: terminal

    call write_info(this, message, level, terminal)
    this%indentation = this%indentation + 1

  end subroutine


  subroutine end_section(this, message, level, terminal)

    class(simulation_log), intent(inout) :: this
    character(*), intent(in) :: message
    integer, intent(in), optional :: level
    logical, intent(in), optional :: terminal

    INSIST(this%indentation > 0)
    this%indentation = this%indentation - 1
    call write_info(this, message, level, terminal)

  end subroutine


  !! Begin/end a silent terminal-output scope. The I/O process decides whether
  !! the group is visible based on the configured interval; all ranks must
  !! call these methods in matching order. Sections inside the group retain
  !! their ordinary output and indentation behavior.
  subroutine begin_terminal_throttle_group(this)

    class(simulation_log), intent(inout) :: this
    real(real64) :: now

    INSIST(.not.this%terminal_throttle_active)
    this%terminal_throttle_active = .true.
    this%terminal_group_visible = .true.
    if (.not.this%io_process) return

    if (.not.this%terminal_output) then
      this%terminal_group_visible = .false.
    else if (this%terminal_step_interval > 0.0_real64) then
      now = MPI_Wtime()
      this%terminal_group_visible = .not.this%have_terminal_step_time .or. &
          now - this%last_terminal_step_time >= this%terminal_step_interval
      if (this%terminal_group_visible) then
        this%last_terminal_step_time = now
        this%have_terminal_step_time = .true.
      end if
    end if

  end subroutine


  subroutine end_terminal_throttle_group(this)

    class(simulation_log), intent(inout) :: this

    INSIST(this%terminal_throttle_active)
    this%terminal_throttle_active = .false.
    this%terminal_group_visible = .true.

  end subroutine


  subroutine write_info(this, message, level, terminal)

    class(simulation_log), intent(in) :: this
    character(*), intent(in) :: message
    integer, intent(in), optional :: level
    logical, intent(in), optional :: terminal

    integer :: message_level
    logical :: write_terminal

    message_level = LOG_NORMAL
    if (present(level)) message_level = level
    if (message_level > this%verbosity .or. .not.this%io_process) return
    write(this%log_unit, '(2a)') indentation(this), trim(message)
    write_terminal = .true.
    if (present(terminal)) write_terminal = terminal
    if (this%terminal_output .and. write_terminal .and. terminal_section_visible(this)) &
      write(output_unit,'(2a)') indentation(this), trim(message)

  end subroutine

  logical function terminal_section_visible(this)
    class(simulation_log), intent(in) :: this
    terminal_section_visible = .not.this%terminal_throttle_active .or. this%terminal_group_visible
  end function

  subroutine warn(this, message, terminal)
    class(simulation_log), intent(in) :: this
    character(*), intent(in) :: message
    logical, intent(in), optional :: terminal
    logical :: write_terminal
    if (.not.this%io_process) return
    write(this%log_unit, '(3a)') indentation(this), 'Warning: ', trim(message)
    write_terminal = .true.
    if (present(terminal)) write_terminal = terminal
    if (this%terminal_output .and. write_terminal) write(output_unit,'(3a)') indentation(this), 'Warning: ', trim(message)
  end subroutine

  subroutine error(this, message, terminal)
    class(simulation_log), intent(in) :: this
    character(*), intent(in) :: message
    logical, intent(in), optional :: terminal
    logical :: write_terminal
    if (.not.this%io_process) return
    write(this%log_unit, '(3a)') indentation(this), 'ERROR: ', trim(message)
    write_terminal = .true.
    if (present(terminal)) write_terminal = terminal
    if (this%terminal_output .and. write_terminal) write(output_unit,'(3a)') indentation(this), 'ERROR: ', trim(message)
  end subroutine

  subroutine finalize(this)
    type(simulation_log), intent(inout) :: this
    call this%close
  end subroutine


  function indentation(this) result(prefix)

    class(simulation_log), intent(in) :: this
    character(:), allocatable :: prefix

    prefix = repeat('  ', this%indentation)

  end function

  subroutine broadcast_errmsg(comm, io_process, errmsg)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(in) :: io_process
    character(:), allocatable, intent(inout) :: errmsg
    integer :: ierr, length
    length = 0
    if (io_process) length = len(errmsg)
    call MPI_Bcast(length, 1, MPI_INTEGER, 0, comm, ierr)
    INSIST(ierr == MPI_SUCCESS)
    if (.not.io_process) allocate(character(length) :: errmsg)
    call MPI_Bcast(errmsg, length, MPI_CHARACTER, 0, comm, ierr)
    INSIST(ierr == MPI_SUCCESS)
  end subroutine

end module simulation_log_type
