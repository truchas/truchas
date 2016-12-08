!!
!! Neil N. Carlson <nnc@lanl.gov>
!! December 2016
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_toolpath_factory

  use toolpath_factory
  use parameter_list_type
  use,intrinsic :: iso_fortran_env, only: r8 => real64
#ifdef NAGFOR
  use,intrinsic :: f90_unix, only: exit
#endif
  implicit none

  integer :: status = 0
  character(:), allocatable :: indir

  call process_command_line(indir, status)
  if (status /= 0) call exit(status)

  call test1
  call test2

  call exit(status)

contains

  !! Test getting the toolpath commands from a string.  Uses (and tests)
  !! the plotfile method to generate data for comparison to reference data

  subroutine test1
    type(parameter_list) :: params
    call params%set('command-string', &
        '[["setflag",1],["dwell",2],&
         &["moverel",[3,0,-4],1,0.5,0.25],&
         &["setflag",3,4],["moverel",[-3,4,0],1,0.5],&
         &["clrflag",1,4],["dwell",1],&
         &["moverel",[3,0,4],1],&
         &["clrflag",3]]')
    call test_aux('test1: ', 'tpfac-test1.dat', indir//'/tpfac-ref1.dat', params)
  end subroutine test1

  !! Test getting the toolpath commands from a file.  Uses (and tests)
  !! the plotfile method to generate data for comparison to reference data

  subroutine test2
    type(parameter_list) :: params
    call params%set('command-file', indir//'/tpfac.json')
    call test_aux('test2: ', 'tpfac-test2.dat', indir//'/tpfac-ref1.dat', params)
  end subroutine test2

  subroutine test_aux(label, tmpfile, reffile, params)

    character(*), intent(in) :: label, tmpfile, reffile
    type(parameter_list) :: params

    type(toolpath), allocatable :: tp
    integer :: stat, lun
    character(:), allocatable :: errmsg, command

    call alloc_toolpath(tp, params, stat, errmsg)
    if (stat /= 0) then
      call write_fail(label//errmsg)
      return
    end if

    open(newunit=lun,file=tmpfile,action='write')
    call tp%write_plotfile(lun, dt=1.0_r8)
    close(lun)

    command = 'diff -q ' // reffile // ' ' // tmpfile
    call execute_command_line(command, exitstat=stat)
    if (stat /= 0) call write_fail(label//'"'//command//'" failed')

  end subroutine test_aux

  subroutine write_fail(errmsg)
    use,intrinsic :: iso_fortran_env, only: error_unit
    character(*), intent(in) :: errmsg
    status = 1
    write(error_unit,'(a)') errmsg
  end subroutine

  subroutine process_command_line(indir, status)

    use,intrinsic :: iso_fortran_env, only: output_unit

    character(:), allocatable, intent(out) :: indir
    integer, intent(out) :: status

    character(:), allocatable :: prog
    character(256) :: arg
    integer :: n

    call get_command_argument(0, arg)
    n = scan(arg, '/', back=.true.) ! remove the leading path component if any
    prog = trim(arg(n+1:))

    status = 0
    select case (command_argument_count())
    case (0)
      indir = '.'
    case (1)
      call get_command_argument(1, arg)
      select case (arg)
      case ('-h', '--help')
        status = 1
      case default
        indir = trim(arg)
      end select
    case default
      status = 1
    end select

    if (status /= 0) then
      write(output_unit,'(a)') 'Usage: ' // prog // ' [INDIR]'
      write(output_unit,'(a)') 'INDIR is the directory containing the input file (default ".")'
    end if

  end subroutine process_command_line

end program test_toolpath_factory
