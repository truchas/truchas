!!
!! TRUCHAS_TIMERS
!!
!! This module provides the subroutines START_TIMER and STOP_TIMER for timing
!! nested sections of code.  Use of the timers requires no more than calling
!! START_TIMER and STOP_TIMER in matching pairs with an arbitrary string
!! argument to name the timer.  The only requirement is that the timers be
!! nested; the timer most recently started (and still running) is the only one
!! eligible to be stopped.
!!
!! The timers come from Petaca's TIMER_TREE_TYPE module.
!!
!! The subroutine WRITE_TIMER_DATA (no arguments) writes the accumulated timer
!! data to the same destinations as TLS_info.  Note that this procedure only
!! writes data from one process under the assumption that each process generates
!! an identical timer tree with essentially the same times (wall clock).  This
!! will be the case if the timers are always started and stopped collectively.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module truchas_timers

  use timer_tree_type, only: start_timer, stop_timer, serialize_timer_tree
  use truchas_logging_services
  implicit none
  private

  public :: start_timer, stop_timer
  public :: write_timer_data

contains

  !! For simplicity each process goes through the motion of writing its
  !! timer tree data, however TLS_info is a no-op on all processes other
  !! than the IO process.

  subroutine write_timer_data

    integer :: j, n, m, indent
    integer, allocatable :: tree(:)
    character(:), allocatable :: name(:)
    real, allocatable :: time(:)
    character(80) :: line

    call TLS_info('')
    call TLS_info('TIMING SUMMARY')
    call TLS_info('--------------')

    call serialize_timer_tree(tree, name, time)

    !! Write the timer tree data using indentation to describe the tree hierarchy.
    !! See the timer_tree_type documentation for how to use the serialized data.
    m = 0 ! high-water mark of node indices encountered
    indent = 0
    do j = 1, size(tree)
      n = tree(j)
      if (n > m) then ! first encounter of this node index
        m = n
        indent = indent + 2
        write(line,'(a,es9.3," --- ",a)') repeat(' ',indent), time(n), trim(name(n))
        call TLS_info(trim(line))
      else
        indent = indent - 2
      end if
    end do
    call TLS_info('')

  end subroutine write_timer_data

end module truchas_timers
