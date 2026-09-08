!!
!! SIMULATION_OUTPUT_SCHEDULE
!!
!! This module expands the simulation output schedule from its parameter-list
!! representation into a strictly increasing array of output times.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, September 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

module simulation_output_schedule

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  implicit none
  private

  public :: get_output_times

contains

  subroutine get_output_times(params, initial_time, output_times, stat, errmsg)

    type(parameter_list), intent(inout) :: params
    real(r8), intent(in) :: initial_time
    real(r8), allocatable, intent(out) :: output_times(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: schedule
    real(r8), allocatable :: times(:), expanded(:)
    real(r8) :: final_time, delta
    integer, allocatable :: subintervals(:)
    integer :: i, j, k, ntime, noutput

    stat = 0

    if (.not.params%is_sublist('output-times')) then
      stat = 1
      errmsg = 'processing ' // params%path() // ': missing "output-times" sublist parameter'
      return
    end if
    schedule => params%sublist('output-times')

    call schedule%get('times', times, stat, errmsg)
    if (stat /= 0) then
      errmsg = 'processing ' // schedule%path() // ': ' // errmsg
      return
    end if
    ntime = size(times)
    if (ntime == 0) then
      stat = 1
      errmsg = 'processing ' // schedule%path() // ': require a nonempty "times" array'
      return
    end if

    if (schedule%is_parameter('subintervals')) then
      call schedule%get('subintervals', subintervals, stat, errmsg)
      if (stat /= 0) then
        errmsg = 'processing ' // schedule%path() // ': ' // errmsg
        return
      end if
    else
      allocate(subintervals(ntime - 1))
      subintervals = 1
    end if

    if (size(subintervals) /= ntime - 1) then
      stat = 1
      errmsg = 'processing ' // schedule%path() // &
          ': require subintervals to have one fewer entry than times'
      return
    end if
    if (ntime > 1) then
      if (any(times(2:) <= times(:ntime-1))) then
        stat = 1
        errmsg = 'processing ' // schedule%path() // ': require times to be strictly increasing'
        return
      end if
    end if
    if (size(subintervals) > 0) then
      if (any(subintervals <= 0)) then
        stat = 1
        errmsg = 'processing ' // schedule%path() // ': require subintervals to be positive'
        return
      end if
    end if

    call params%get('final-time', final_time, default=times(ntime), stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      errmsg = 'processing ' // params%path() // ': ' // errmsg
      return
    end if
    if (final_time <= initial_time) then
      stat = 1
      errmsg = 'processing ' // params%path() // ': require final-time > initial-time'
      return
    end if
    if (final_time < times(ntime)) then
      stat = 1
      errmsg = 'processing ' // params%path() // &
          ': require final-time to be >= the last output schedule time'
      return
    end if

    allocate(expanded(1+sum(subintervals)))
    expanded(1) = times(1)
    k = 1
    do i = 1, ntime - 1
      delta = (times(i+1) - times(i)) / real(subintervals(i), r8)
      do j = 1, subintervals(i)
        k = k + 1
        if (j == subintervals(i)) then
          expanded(k) = times(i+1)
        else
          expanded(k) = times(i) + real(j, r8) * delta
        end if
      end do
    end do

    noutput = count(expanded > initial_time)
    if (final_time > times(ntime)) noutput = noutput + 1
    if (noutput == 0) then
      stat = 1
      errmsg = 'processing ' // params%path() // ': require an output time after initial-time'
      return
    end if
    allocate(output_times(noutput))
    k = 0
    do i = 1, size(expanded)
      if (expanded(i) > initial_time) then
        k = k + 1
        output_times(k) = expanded(i)
      end if
    end do
    if (final_time > times(ntime)) then
      k = k + 1
      output_times(k) = final_time
    end if

  end subroutine get_output_times

end module simulation_output_schedule
