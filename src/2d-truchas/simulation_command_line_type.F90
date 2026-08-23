!!
!! SIMULATION_COMMAND_LINE_TYPE
!!
!! This module defines the common command-line interface used by the JSON
!! driven simulation programs.  It handles the input file, output directory,
!! overwrite policy, and help request.  Parsing is MPI-independent; the
!! drivers initialize MPI before reporting help or command-line errors so
!! those messages are emitted only by the I/O process.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

module simulation_command_line_type

  use,intrinsic :: iso_fortran_env, only: output_unit
  use utilities_module, only: make_directory_hierarchy
  implicit none
  private

  type, public :: simulation_command_line
    character(:), allocatable :: program
    character(:), allocatable :: input_file
    character(:), allocatable :: input_dir
    character(:), allocatable :: output_dir
    logical :: force = .false.
    logical :: help = .false.
  contains
    procedure :: parse
    procedure :: prepare_output_dir
    procedure :: write_help
  end type simulation_command_line

contains

  subroutine parse(this, stat, errmsg)

    class(simulation_command_line), intent(out) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: i, narg
    character(4096) :: arg

    stat = 0
    call get_command_argument(0, arg)
    this%program = basename(trim(arg))
    if (len_trim(this%program) == 0) this%program = 'simulation'

    narg = command_argument_count()
    i = 0
    do while (i < narg)
      i = i + 1
      call get_command_argument(i, arg)
      select case (trim(arg))
      case ('-h', '--help')
        this%help = .true.
      case ('-f', '--force')
        this%force = .true.
      case ('-o', '--output-dir')
        if (i == narg) then
          stat = 1
          errmsg = trim(arg) // ' requires an argument'
          return
        end if
        i = i + 1
        call get_command_argument(i, arg)
        if (len_trim(arg) == 0) then
          stat = 1
          errmsg = 'output directory must not be empty'
          return
        end if
        this%output_dir = trim(arg)
      case default
        if (len_trim(arg) > 0 .and. arg(1:1) == '-') then
          stat = 1
          errmsg = 'invalid option: ' // trim(arg)
          return
        end if
        if (allocated(this%input_file)) then
          stat = 1
          errmsg = 'more than one input file was specified'
          return
        end if
        this%input_file = trim(arg)
      end select
    end do

    if (this%help) return
    if (.not.allocated(this%input_file)) then
      stat = 1
      errmsg = 'an input file is required'
      return
    end if
    this%input_dir = dirname(this%input_file)
    if (.not.allocated(this%output_dir)) then
      this%output_dir = file_stem(this%input_file)
    end if

  end subroutine parse


  subroutine prepare_output_dir(this, stat, errmsg)

    class(simulation_command_line), intent(in) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    logical :: exists

    stat = 0
    inquire(file=trim(this%output_dir), exist=exists)
    if (exists .and. .not.this%force) then
      stat = 1
      errmsg = 'output directory already exists; use --force: ' // &
          trim(this%output_dir)
    else
      stat = make_directory_hierarchy(trim(this%output_dir))
      if (stat /= 0) then
        errmsg = 'unable to create output directory: ' // &
            trim(this%output_dir)
      end if
    end if

  end subroutine prepare_output_dir


  subroutine write_help(this, description)

    class(simulation_command_line), intent(in) :: this
    character(*), intent(in) :: description

    write(output_unit,'(a)') 'Usage: ' // trim(this%program) // &
        ' [options] INFILE'
    write(output_unit,'(a)') ''
    write(output_unit,'(a)') trim(description)
    write(output_unit,'(a)') ''
    write(output_unit,'(a)') 'Options:'
    write(output_unit,'(a)') '  -o DIR, --output-dir DIR'
    write(output_unit,'(a)') '                   Write run products in DIR.'
    write(output_unit,'(a)') '  -f, --force       Permit use of an existing output directory.'
    write(output_unit,'(a)') '  -h, --help        Display this help and exit.'

  end subroutine write_help


  pure function basename(path) result(name)

    character(*), intent(in) :: path
    character(:), allocatable :: name
    integer :: i

    i = scan(path, '/', back=.true.)
    if (i == 0) then
      name = trim(path)
    else if (i == len_trim(path)) then
      name = ''
    else
      name = trim(path(i+1:))
    end if

  end function basename


  pure function dirname(path) result(name)

    character(*), intent(in) :: path
    character(:), allocatable :: name
    integer :: i

    i = scan(trim(path), '/', back=.true.)
    if (i == 0) then
      name = './'
    else
      name = trim(path(:i))
    end if

  end function dirname


  pure function file_stem(path) result(stem)

    character(*), intent(in) :: path
    character(:), allocatable :: stem
    character(:), allocatable :: name
    integer :: i

    name = basename(path)
    i = scan(name, '.', back=.true.)
    if (i > 1) then
      stem = name(:i-1)
    else
      stem = name
    end if
    if (len_trim(stem) == 0) stem = 'output'

  end function file_stem

end module simulation_command_line_type
