!!
!! Neil N. Carlson <nnc@lanl.gov>
!! July 2017
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_ded_head_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use,intrinsic :: iso_fortran_env, only: output_unit
#ifdef NAGFOR
  use,intrinsic :: f90_unix, only: exit
#endif
  use parameter_list_type
  use toolpath_type
  use toolpath_factory
  use toolpath_table
  use ded_head_type
  implicit none

  real(r8), parameter :: PI = 3.141592653589793_r8
  integer :: status = 0

  call make_test_toolpaths
  call test_no_fade
  call test_fade
  call test_fade_init

  call exit(status)

contains

  subroutine make_test_toolpaths

    type(parameter_list) :: params
    type(toolpath), allocatable :: tp
    integer :: stat
    character(:), allocatable :: errmsg

    call params%set('command-string', &
        '[["setflag",0],["dwell",2],["clrflag",0],["dwell",2],["setflag",0],["dwell",2]]')
    call alloc_toolpath(tp, params, stat, errmsg)
    if (stat /= 0) then
      call write_fail(errmsg)
      call exit(status)
    end if
    call insert_toolpath('tp1', tp)

  end subroutine make_test_toolpaths

  !! This tests the exponential modulation of the laser power when switching
  !! between on/off states.  The laser parameters conspire to give a peak
  !! amplitude of 1 so that evaluating its irradiance at its center (origin)
  !! effectively gives the modulation factor.

  subroutine test_fade

    integer :: n, stat
    type(parameter_list) :: params
    type(parameter_list), pointer :: laser_params
    type(ded_head) :: head
    real(r8) :: t, r(3) = 0.0_r8, ref(0:6), error

    !! Parameters for DED_HEAD object initialization
    call params%set('toolpath', 'tp1')
    call params%set('laser-absorp', 1.0_r8)
    call params%set('laser-time-constant', 0.5_r8)
    laser_params => params%sublist('laser')
    call laser_params%set('type', 'gaussian')
    call laser_params%set('power', 2.0_r8*PI)
    call laser_params%set('sigma', 1.0_r8)

    do n = 0, 2
      t = real(n,r8)
      ref(n) = 1 - exp(-t/0.5_r8)
    end do

    do n = 3, 4
      t = real(n,r8)
      ref(n) = ref(2)*exp((2-t)/0.5_r8)
    end do

    do n = 5, 6
      t = real(n,r8)
      ref(n) = 1 - (1-ref(4))*exp((4-t)/0.5_r8)
    end do

    t = 0.0_r8
    call head%init(params, t)

    stat = 0
    do n = 0, 2
      t = real(n,r8)
      error = abs(head%laser_irrad(t, r) - ref(n))
      write(output_unit,'(a,f5.1,a,es10.2)') 't=', t, ', error=', error
      if (error > epsilon(1.0_r8)) stat = 1
    end do

    call head%next_tp_segment(t)
    do n = 2, 4
      t = real(n,r8)
      error = abs(head%laser_irrad(t, r) - ref(n))
      write(output_unit,'(a,f5.1,a,es10.2)') 't=', t, ', error=', error
      if (error > epsilon(1.0_r8)) stat = 1
    end do

    call head%next_tp_segment(t)
    do n = 4, 6
      t = real(n,r8)
      error = abs(head%laser_irrad(t, r) - ref(n))
      write(output_unit,'(a,f5.1,a,es10.2)') 't=', t, ', error=', error
      if (error > epsilon(1.0_r8)) stat = 1
    end do

    if (stat /= 0) then
      call write_fail('test_fade: some errors exceed tolerance')
    else
      write(output_unit,'(a)') 'test_fade: OKAY'
    end if

  end subroutine test_fade

  !! This tests the discontinuous switching of the laser power between on/off
  !! states.  The laser parameters conspire to give a peak amplitude of 1 so
  !! that evaluating its irradiance at its center (origin) should give either
  !! 0 or 1 depending on whether the laser is off or on, respectively.

  subroutine test_no_fade

    integer :: n, stat
    type(parameter_list) :: params
    type(parameter_list), pointer :: laser_params
    type(ded_head) :: head
    real(r8) :: t, r(3) = 0.0_r8, error

    !! Parameters for DED_HEAD object initialization
    call params%set('toolpath', 'tp1')
    call params%set('laser-absorp', 1.0_r8)
    call params%set('laser-time-constant', 0.0_r8)
    laser_params => params%sublist('laser')
    call laser_params%set('type', 'gaussian')
    call laser_params%set('power', 2.0_r8*PI)
    call laser_params%set('sigma', 1.0_r8)

    t = 0.0_r8
    call head%init(params, t)

    stat = 0
    do n = 0, 2
      t = real(n,r8)
      error = abs(head%laser_irrad(t, r) - 1.0_r8)
      write(output_unit,'(a,f5.1,a,es10.2)') 't=', t, ', error=', error
      if (error > epsilon(1.0_r8)) stat = 1
    end do

    call head%next_tp_segment(t)
    do n = 2, 4
      t = real(n,r8)
      error = abs(head%laser_irrad(t, r))
      write(output_unit,'(a,f5.1,a,es10.2)') 't=', t, ', error=', error
      if (error > epsilon(1.0_r8)) stat = 1
    end do

    call head%next_tp_segment(t)
    do n = 4, 6
      t = real(n,r8)
      error = abs(head%laser_irrad(t, r) - 1.0_r8)
      write(output_unit,'(a,f5.1,a,es10.2)') 't=', t, ', error=', error
      if (error > epsilon(1.0_r8)) stat = 1
    end do

    if (stat /= 0) then
      call write_fail('test_no_fade: some errors exceed tolerance')
    else
      write(output_unit,'(a)') 'test_no_fade: OKAY'
    end if

  end subroutine test_no_fade

  !! This tests the initialization of the exponential modulation of laser
  !! when the initial time does not coincide with the start of a path segment
  !! where the laser switched state between on/off. The laser parameters
  !! conspire to give a peak amplitude of 1 so that evaluating its irradiance
  !! at its center (origin) effectively gives the modulation factor.

  subroutine test_fade_init

    type(parameter_list) :: params
    type(parameter_list), pointer :: laser_params
    type(ded_head) :: head
    real(r8) :: t, r(3) = 0.0_r8, ref, error
    integer :: stat

    !! Parameters for DED_HEAD object initialization
    call params%set('toolpath', 'tp1')
    call params%set('laser-absorp', 1.0_r8)
    call params%set('laser-time-constant', 0.5_r8)
    laser_params => params%sublist('laser')
    call laser_params%set('type', 'gaussian')
    call laser_params%set('power', 2.0_r8*PI)
    call laser_params%set('sigma', 1.0_r8)

    stat = 0

    !! Off-to-on at t=0
    call head%init(params, t=1.0_r8)
    t = 1.0_r8
    ref = 1 - exp(-t/0.5_r8)
    error = abs(head%laser_irrad(t, r) - ref)
    write(output_unit,'(a,f5.1,a,es10.2)') 't=', t, ', error=', error
    if (error > epsilon(1.0_r8)) then
      call write_fail('test_fade_init: off-to-on: error exceeds tolerance')
      stat = 1
    end if

    !! On-to-off at t=2
    call head%init(params, t=3.0_r8)
    t = 3.0_r8
    ref = exp((2-t)/0.5_r8)
    error = abs(head%laser_irrad(t, r) - ref)
    write(output_unit,'(a,f5.1,a,es10.2)') 't=', t, ', error=', error
    if (error > epsilon(1.0_r8)) then
      call write_fail('test_fade_init: on-to-off: error exceeds tolerance')
      stat = 1
    end if

    if (stat == 0) then
      write(output_unit,'(a)') 'test_fade_init: OKAY'
    else
      write(output_unit,'(a)') 'test_fade_init: FAIL'
    end if

  end subroutine test_fade_init


  subroutine write_fail(errmsg)
    use,intrinsic :: iso_fortran_env, only: output_unit
    character(*), intent(in) :: errmsg
    status = 1
    write(output_unit,'(a)') errmsg
  end subroutine

end program test_ded_head_type
