!!
!! Neil N. Carlson <nnc@lanl.gov>
!! December 2016
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_toolpath_type

  use toolpath_type
  use,intrinsic :: iso_fortran_env, only: r8 => real64
#ifdef NAGFOR
  use,intrinsic :: f90_unix, only: exit
#endif
  implicit none

  integer :: status = 0

  call test_start_seg
  call test_final_seg
  call test_minimal
  call test_path
  call test_partition

  call exit(status)

contains

  !! Test appending a start segment.  This does not produce a valid toolpath
  !! but we can still position the tool path to the segment and query data
  !! from the segment.

  subroutine test_start_seg

    type(toolpath) :: tp
    real(r8) :: t1 = 1, r1(3) = [1,2,3], r(3)
    integer :: flags = 2

    call tp%append_path_segment(new_start_segment(t1, r1, flags))
    call tp%set_segment(-1.0_r8)

    call tp%get_position(-1.0_r8, r)
    if (any(r /= r1)) call write_fail('test_start_seg: wrong position')

    if (tp%is_flag_set(0)) call write_fail('test_start_seg: flag 0 wrong')
    if (.not.tp%is_flag_set(1)) call write_fail('test_start_seg: flag 1 wrong')

    if (tp%is_valid()) call write_fail('test_start_seg: validity check failed')

  end subroutine test_start_seg

  !! Test appending a final segment.  This does not produce a valid toolpath
  !! but we can still position the tool path to the segment and query data
  !! from the segment.

  subroutine test_final_seg

    type(toolpath) :: tp
    real(r8) :: t0 = 1, r0(3) = [1,2,3], r(3)
    integer :: flags = 2

    call tp%append_path_segment(new_final_segment(t0, r0, flags))
    call tp%set_segment(2.0_r8)

    call tp%get_position(2.0_r8, r)
    if (any(r /= r0)) call write_fail('test_final_seg: wrong position')

    if (tp%is_flag_set(0)) call write_fail('test_final_seg: flag 0 wrong')
    if (.not.tp%is_flag_set(1)) call write_fail('test_final_seg: flag 1 wrong')

    if (tp%is_valid()) call write_fail('test_final_seg: validity check failed')

  end subroutine test_final_seg

  !! Test a minimally valid toolpath (start and final segments)

  subroutine test_minimal

    type(toolpath) :: tp
    real(r8) :: t1 = 1, r1(3) = [1,2,3], r0(3) = [3,2,1], r(3)

    !! A discontinuous path for checking correct toolpath positioning
    call tp%append_path_segment(new_start_segment(t1, r1, flags=2))
    call tp%append_path_segment(new_final_segment(t1, r0, flags=1))

    if (.not.tp%is_valid()) call write_fail('test_minimal: validity check failed')

    !! Check final segment
    call tp%set_segment(2.0_r8)
    call tp%get_position(t1, r)
    if (any(r /= r0)) call write_fail('test_minimal: final wrong position')
    if (.not.tp%is_flag_set(0)) call write_fail('test_minimal: final flag 0 wrong')
    if (tp%is_flag_set(1)) call write_fail('test_minimal: final flag 1 wrong')

    !! Check start segment
    call tp%set_segment(0.0_r8)
    call tp%get_position(t1, r)
    if (any(r /= r1)) call write_fail('test_minimal: start wrong position')
    if (tp%is_flag_set(0)) call write_fail('test_minimal: start flag 0 wrong')
    if (.not.tp%is_flag_set(1)) call write_fail('test_minimal: start flag 1 wrong')

    !! Move to final segment and check again
    call tp%next_segment
    call tp%get_position(t1, r)
    if (any(r /= r0)) call write_fail('test_minimal: next wrong position')
    if (.not.tp%is_flag_set(0)) call write_fail('test_minimal: next flag 0 wrong')
    if (tp%is_flag_set(1)) call write_fail('test_minimal: next flag 1 wrong')

  end subroutine test_minimal

  !! Check a hand-built toolpath.  More checks of toolpath segment positioning
  !! and adds check of new_path_segment and get_segment starts.

  subroutine test_path

    use dwell_xyz_motion_type
    use linear_xyz_motion_type
    use xyz_motion_class

    type(toolpath) :: tp
    class(xyz_motion), allocatable :: move

    real(r8) :: t = 1, r(3) = [1,2,3]
    real(r8), allocatable :: times(:)
    logical, allocatable :: discont(:)

    call tp%append_path_segment(new_start_segment(t, r, flags=0))
#ifdef NO_2008_LHS_POLY_REALLOC
    allocate(move, source=dwell_xyz_motion(r, t0=t, dt=1.0_r8))
#else
    move = dwell_xyz_motion(r, t0=t, dt=1.0_r8)
#endif
    t = move%final_time(); r = move%final_coord()
    call tp%append_path_segment(new_path_segment(move, flags=1))
#ifdef NO_2008_LHS_POLY_REALLOC
    allocate(move, source=linear_xyz_motion(t, r, [1.0_r8, 0.0_r8, 0.0_r8], 1.0_r8))
#else
    move = linear_xyz_motion(t, r, [1.0_r8, 0.0_r8, 0.0_r8], 1.0_r8)
#endif
    t = move%final_time(); r = move%final_coord()
    call tp%append_path_segment(new_path_segment(move, flags=1))
#ifdef NO_2008_LHS_POLY_REALLOC
    allocate(move, source=dwell_xyz_motion(r, t0=t, dt=1.0_r8))
#else
    move = dwell_xyz_motion(r, t0=t, dt=1.0_r8)
#endif
    t = move%final_time(); r = move%final_coord()
    call tp%append_path_segment(new_path_segment(move, flags=0))
    call tp%append_path_segment(new_final_segment(t, r, flags=0))

    if (.not.tp%is_valid()) call write_fail('test_path: validity check failed')

    call tp%set_segment(1.0_r8) ! the dwell, not start segment
    if (.not.tp%is_flag_set(0)) call write_fail('test_path: A: flag is wrong')
    call tp%get_position(1.0_r8,r)
    if (any(r /= [1,2,3])) call write_fail('test_path: A: start position is wrong')
    call tp%get_position(2.0_r8,r)
    if (any(r /= [1,2,3])) call write_fail('test_path: A: final position is wrong')

    call tp%next_segment  ! the linear move
    if (.not.tp%is_flag_set(0)) call write_fail('test_path: B: flag is wrong')
    call tp%get_position(2.0_r8,r)
    if (any(r /= [1,2,3])) call write_fail('test_path: B: start position is wrong')
    call tp%get_position(3.0_r8,r)
    if (any(r /= [2,2,3])) call write_fail('test_path: B: final position is wrong')

    call tp%next_segment  ! the dwell
    if (tp%is_flag_set(0)) call write_fail('test_path: C: flag is wrong')
    call tp%get_position(3.0_r8,r)
    if (any(r /= [2,2,3])) call write_fail('test_path: C: start position is wrong')
    call tp%get_position(4.0_r8,r)
    if (any(r /= [2,2,3])) call write_fail('test_path: C: final position is wrong')

    call tp%next_segment  ! the final segment
    if (tp%is_flag_set(0)) call write_fail('test_path: D: flag is wrong')
    call tp%get_position(4.0_r8,r)
    if (any(r /= [2,2,3])) call write_fail('test_path: D: start position is wrong')

    call tp%get_segment_starts(times, discont)
    if (size(times) /= 4) then
      call write_fail('test_path: times wrong size')
    else if (any(times /= [1,2,3,4])) then
      call write_fail('test_path: times is wrong')
    end if
    if (size(discont) /= 4) then
      call write_fail('test_path: discont wrong size')
    else if (any(discont .neqv. [.true.,.false.,.true.,.false.])) then
      call write_fail('test_path: discont is wrong')
    end if

  end subroutine test_path

  !! This checks the special partitioning feature using another hand-built toolpath.

  subroutine test_partition

    use dwell_xyz_motion_type
    use linear_xyz_motion_type
    use xyz_motion_class

    type(toolpath) :: tp
    class(xyz_motion), allocatable :: move

    real(r8) :: t = 0, r(3) = [0,0,0]
    real(r8), allocatable :: time(:), coord(:,:)
    character(:), allocatable :: hash(:)

    call tp%append_path_segment(new_start_segment(t, r, flags=0))
#ifdef NO_2008_LHS_POLY_REALLOC
    allocate(move, source=dwell_xyz_motion(r, t0=t, dt=1.0_r8))
#else
    move = dwell_xyz_motion(r, t0=t, dt=1.0_r8)
#endif
    t = move%final_time(); r = move%final_coord()
    call tp%append_path_segment(new_path_segment(move, flags=0))
#ifdef NO_2008_LHS_POLY_REALLOC
    allocate(move, source=linear_xyz_motion(t, r, [0.0_r8, 1.0_r8, 0.0_r8], 1.0_r8))
#else
    move = linear_xyz_motion(t, r, [0.0_r8, 1.0_r8, 0.0_r8], 1.0_r8)
#endif
    t = move%final_time(); r = move%final_coord()
    call tp%append_path_segment(new_path_segment(move, flags=0))
#ifdef NO_2008_LHS_POLY_REALLOC
    allocate(move, source=linear_xyz_motion(t, r, [2.0_r8, 0.0_r8, 0.0_r8], 2.0_r8))
#else
    move = linear_xyz_motion(t, r, [2.0_r8, 0.0_r8, 0.0_r8], 2.0_r8)
#endif
    t = move%final_time(); r = move%final_coord()
    call tp%append_path_segment(new_path_segment(move, flags=0))
#ifdef NO_2008_LHS_POLY_REALLOC
    allocate(move, source=linear_xyz_motion(t, r, [0.0_r8, 0.0_r8, 4.0_r8], 2.0_r8))
#else
    move = linear_xyz_motion(t, r, [0.0_r8, 0.0_r8, 4.0_r8], 2.0_r8)
#endif
    t = move%final_time(); r = move%final_coord()
    call tp%append_path_segment(new_path_segment(move, flags=0))
    call tp%append_path_segment(new_final_segment(t, r, flags=0))

    if (.not.tp%is_valid()) call write_fail('test_partition: validity check failed')

    call tp%set_partition(2.0_r8)
    call tp%get_partition(time=time, coord=coord, hash=hash)
    if (size(time) /= 6) then
      call write_fail('test_partition: A: wrong time size')
    else if (any(time /= [0,1,2,3,4,5])) then
      call write_fail('test_partition: A: wrong time values')
    end if
    call tp%get_partition(time=time, coord=coord, hash=hash)
    if (any(shape(coord) /= [3,6])) then
      call write_fail('test_partition: A: wrong coord shape')
    else if (any(coord /= reshape([0,0,0, 0,0,0, 0,1,0, 2,1,0, 2,1,2, 2,1,4],[3,6]))) then
      call write_fail('test_partition: A: wrong coord values')
    end if
    if (size(hash) /= 6) then
      call write_fail('test_partition: A: wrong hash size')
    else if (any(hash /= ['d3399b7','d3399b7','e28bd59','ecb21b7','81e2bcf','e2b04c2'])) then
      call write_fail('test_partition: A: wrong hash values')
    end if

  end subroutine test_partition

  subroutine write_fail(errmsg)
    use,intrinsic :: iso_fortran_env, only: error_unit
    character(*), intent(in) :: errmsg
    status = 1
    write(error_unit,'(a)') errmsg
  end subroutine

end program test_toolpath_type
