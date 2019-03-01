!!
!! Neil N. Carlson <nnc@lanl.gov>
!! December 2016
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_linear_xyz_motion_type

  use linear_xyz_motion_type
  use,intrinsic :: iso_fortran_env, only: r8 => real64
#ifdef NAGFOR
  use,intrinsic :: f90_unix, only: exit
#endif
  implicit none

  integer :: status = 0

  call test_sad
  call test_sa
  call test_s
  call test_sad_crit
  call test_sa_crit

  call exit(status)

contains

  subroutine test_sad

    type(linear_xyz_motion) :: move
    real(r8) :: t0 = 10, r0(3)=[1,2,3], dr(3)=[2,3,6], s=1, a=0.5_r8, d=0.25_r8, t1 = 20
    real(r8), allocatable :: times(:)

    move = linear_xyz_motion(t0, r0, dr, s, a, d)

    if (move%start_time() /= t0) call write_fail('test_sad: wrong start time')
    if (move%final_time() /= t1) call write_fail('test_sad: wrong final time')
    if (any(move%start_coord() /= r0)) call write_fail('test_sad: wrong start coord')
    if (any(move%final_coord() /= r0+dr)) call write_fail('test_sad: wrong final coord')
    if (any(move%coord(t0) /= r0)) call write_fail('test_sad: wrong coord left')
    if (any(move%coord(t1) /= r0+dr)) call write_fail('test_sad: wrong coord right')

    if (any(move%coord(t0+2) /= r0+dr/7)) call write_fail('test_sad: wrong coord at speed')
    if (any(move%coord(t0+6) /= r0+5*dr/7)) call write_fail('test_sad: wrong coord at decel')

    times = move%partition(6.99_r8)
    if (size(times) == 3) then
      if (any(times /= [t0,t0+4.5,t1])) call write_fail('test_sad: wrong partition 1')
    else
      call write_fail('test_sad: wrong partition 1 size')
    end if

    times = move%partition(7.01_r8)
    if (size(times) == 2) then
      if (any(times /= [t0,t1])) call write_fail('test_sad: wrong partition 2')
    else
      call write_fail('test_sad: wrong partition 2 size')
    end if

  end subroutine test_sad

  subroutine test_sa

    type(linear_xyz_motion) :: move
    real(r8) :: t0 = 10, r0(3)=[1,2,3], dr(3)=[2,3,6], s=1, a=0.25_r8, t1 = 21
    real(r8), allocatable :: times(:)

    move = linear_xyz_motion(t0, r0, dr, s, a)

    if (move%start_time() /= t0) call write_fail('test_sa: wrong start time')
    if (move%final_time() /= t1) call write_fail('test_sa: wrong final time')
    if (any(move%start_coord() /= r0)) call write_fail('test_sa: wrong start coord')
    if (any(move%final_coord() /= r0+dr)) call write_fail('test_sa: wrong final coord')
    if (any(move%coord(t0) /= r0)) call write_fail('test_sa: wrong coord left')
    if (any(move%coord(t1) /= r0+dr)) call write_fail('test_sa: wrong coord right')

    if (any(move%coord(t0+4) /= r0+2*dr/7)) call write_fail('test_sa: wrong coord at speed')
    if (any(move%coord(t0+7) /= r0+5*dr/7)) call write_fail('test_sa: wrong coord at decel')

    times = move%partition(6.99_r8)
    if (size(times) == 3) then
      if (any(times /= [t0,t0+5.5,t1])) call write_fail('test_sa: wrong partition 1')
    else
      call write_fail('test_sa: wrong partition 1 size')
    end if

    times = move%partition(7.01_r8)
    if (size(times) == 2) then
      if (any(times /= [t0,t1])) call write_fail('test_sa: wrong partition 2')
    else
      call write_fail('test_sa: wrong partition 2 size')
    end if

  end subroutine test_sa

  subroutine test_s

    type(linear_xyz_motion) :: move
    real(r8) :: t0 = 10, r0(3)=[1,2,3], dr(3)=[2,3,6], s=1, t1 = 17
    real(r8), allocatable :: times(:)

    move = linear_xyz_motion(t0, r0, dr, s)

    if (move%start_time() /= t0) call write_fail('test_s: wrong start time')
    if (move%final_time() /= t1) call write_fail('test_s: wrong final time')
    if (any(move%start_coord() /= r0)) call write_fail('test_s: wrong start coord')
    if (any(move%final_coord() /= r0+dr)) call write_fail('test_s: wrong final coord')
    if (any(move%coord(t0) /= r0)) call write_fail('test_s: wrong coord left')
    if (any(move%coord(t1) /= r0+dr)) call write_fail('test_s: wrong coord right')

    if (any(move%coord(t0+1) /= r0+dr/7)) call write_fail('test_s: wrong coord at 11')
    if (any(move%coord(t0+6) /= r0+6*dr/7)) call write_fail('test_s: wrong coord at 16')

    times = move%partition(6.99_r8)
    if (size(times) == 3) then
      if (any(times /= [t0,t0+3.5,t1])) call write_fail('test_s: wrong partition 1')
    else
      call write_fail('test_s: wrong partition 1 size')
    end if

    times = move%partition(7.01_r8)
    if (size(times) == 2) then
      if (any(times /= [t0,t1])) call write_fail('test_s: wrong partition 2')
    else
      call write_fail('test_s: wrong partition 2 size')
    end if

  end subroutine test_s

  subroutine test_sad_crit

    type(linear_xyz_motion) :: move
    real(r8) :: t0 = 10, r0(3)=[0,0,1], dr(3)=[0,0,3], s=2+epsilon(2.0_r8), a=2, d=1, t1=13

    move = linear_xyz_motion(t0, r0, dr, s, a, d)

    if (move%start_time() /= t0) call write_fail('test_sad_crit: wrong start time')
    if (move%final_time() /= t1) call write_fail('test_sad_crit: wrong final time')
    if (any(move%start_coord() /= r0)) call write_fail('test_sad_crit: wrong start coord')
    if (any(move%final_coord() /= r0+dr)) call write_fail('test_sad_crit: wrong final coord')
    if (any(move%coord(t0) /= r0)) call write_fail('test_sad_crit: wrong coord left')
    if (any(move%coord(t1) /= r0+dr)) call write_fail('test_sad_crit: wrong coord right')

    if (any(move%coord(10.5_r8) /= r0+dr/12)) call write_fail('test_sad_crit: wrong coord at 10.5')
    if (any(move%coord(11.0_r8) /= r0+dr/3)) call write_fail('test_sad_crit: wrong coord at 11.0')
    if (any(move%coord(12.0_r8) /= r0+5*dr/6)) call write_fail('test_sad_crit: wrong coord at 12.0')

  end subroutine test_sad_crit

  subroutine test_sa_crit

    type(linear_xyz_motion) :: move
    real(r8) :: t0 = 10, r0(3)=[0,0,1], dr(3)=[0,0,2], s=2+epsilon(2.0_r8), a=2, t1=12

    move = linear_xyz_motion(t0, r0, dr, s, a)

    if (move%start_time() /= t0) call write_fail('test_sa_crit: wrong start time')
    if (move%final_time() /= t1) call write_fail('test_sa_crit: wrong final time')
    if (any(move%start_coord() /= r0)) call write_fail('test_sa_crit: wrong start coord')
    if (any(move%final_coord() /= r0+dr)) call write_fail('test_sa_crit: wrong final coord')
    if (any(move%coord(t0) /= r0)) call write_fail('test_sa_crit: wrong coord left')
    if (any(move%coord(t1) /= r0+dr)) call write_fail('test_sa_crit: wrong coord right')

    if (any(move%coord(10.5_r8) /= r0+dr/8)) call write_fail('test_sa_crit: wrong coord at 10.5')
    if (any(move%coord(11.0_r8) /= r0+dr/2)) call write_fail('test_sa_crit: wrong coord at 11.0')
    if (any(move%coord(11.5_r8) /= r0+7*dr/8)) call write_fail('test_sa_crit: wrong coord at 11.5')

  end subroutine test_sa_crit

  subroutine write_fail(errmsg)
    use,intrinsic :: iso_fortran_env, only: error_unit
    character(*), intent(in) :: errmsg
    status = 1
    write(error_unit,'(a)') errmsg
  end subroutine

end program test_linear_xyz_motion_type
