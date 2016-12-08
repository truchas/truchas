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

program test_dwell_xyz_motion_type

  use dwell_xyz_motion_type
  use,intrinsic :: iso_fortran_env, only: r8 => real64
#ifdef NAGFOR
  use,intrinsic :: f90_unix, only: exit
#endif
  implicit none

  integer :: status = 0

  call test_usual
  call test_start
  call test_final

  call exit(status)

contains

  subroutine test_usual

    type(dwell_xyz_motion) :: move
    real(r8) :: t0 = 10, r(3)=[1,2,3], dt=2
    real(r8), allocatable :: times(:)

    move = dwell_xyz_motion(r, t0=t0, dt=dt)

    if (move%start_time() /= t0) call write_fail('test_usual: wrong start time')
    if (move%final_time() /= t0+dt) call write_fail('test_usual: wrong final time')
    if (any(move%start_coord() /= r)) call write_fail('test_usual: wrong start coord')
    if (any(move%final_coord() /= r)) call write_fail('test_usual: wrong final coord')
    if (any(move%coord(t0) /= r)) call write_fail('test_usual: wrong coord left')
    if (any(move%coord(t0+dt/2) /= r)) call write_fail('test_usual: wrong coord mid')
    if (any(move%coord(t0+dt) /= r)) call write_fail('test_usual: wrong coord right')

    times = move%partition(0.1_r8)
    if (size(times) == 2) then
      if (any(times /= [t0,t0+dt])) call write_fail('test_usual: wrong partition')
    else
      call write_fail('test_usual: wrong partition size')
    end if

  end subroutine test_usual

  subroutine test_start

    type(dwell_xyz_motion) :: move
    real(r8) :: t1 = 10, r(3)=[1,2,3]
    real(r8), allocatable :: times(:)

    move = dwell_xyz_motion(r, t1=t1)

    if (move%start_time() /= -huge(t1)) call write_fail('test_start: wrong start time')
    if (move%final_time() /= t1) call write_fail('test_start: wrong final time')
    if (any(move%start_coord() /= r)) call write_fail('test_start: wrong start coord')
    if (any(move%final_coord() /= r)) call write_fail('test_start: wrong final coord')
    if (any(move%coord(t1-1) /= r)) call write_fail('test_start: wrong coord mid')
    if (any(move%coord(t1) /= r)) call write_fail('test_start: wrong coord right')

    times = move%partition(0.1_r8)
    if (size(times) == 1) then
      if (any(times /= [t1])) call write_fail('test_start: wrong partition')
    else
      call write_fail('test_start: wrong partition size')
    end if

  end subroutine test_start

  subroutine test_final

    type(dwell_xyz_motion) :: move
    real(r8) :: t0 = 10, r(3)=[1,2,3]
    real(r8), allocatable :: times(:)

    move = dwell_xyz_motion(r, t0=t0)

    if (move%start_time() /= t0) call write_fail('test_final: wrong start time')
    if (move%final_time() /= huge(t0)) call write_fail('test_final: wrong final time')
    if (any(move%start_coord() /= r)) call write_fail('test_final: wrong start coord')
    if (any(move%final_coord() /= r)) call write_fail('test_final: wrong final coord')
    if (any(move%coord(t0) /= r)) call write_fail('test_final: wrong coord left')
    if (any(move%coord(t0+1) /= r)) call write_fail('test_final: wrong coord mid')

    times = move%partition(0.1_r8)
    if (size(times) == 1) then
      if (any(times /= [t0])) call write_fail('test_final: wrong partition')
    else
      call write_fail('test_final: wrong partition size')
    end if

  end subroutine test_final

  subroutine write_fail(errmsg)
    use,intrinsic :: iso_fortran_env, only: error_unit
    character(*), intent(in) :: errmsg
    status = 1
    write(error_unit,'(a)') errmsg
  end subroutine

end program test_dwell_xyz_motion_type
