!!===================================================================================
!!TRUCHAS TIMING
!! 
!!PURPOSE(S):
!! 
!!    Manage all the timing output for Truchas. Use TBrook to print XML information,
!!    and print more detailed output to the console and to the .out file. Uses 
!!    Timing_tree to get timing data.
!!
!! CONTAINS:
!!    report_tree_timing()  !Print timing data.
!!
!! AUTHOR(S):
!!    Brandon S. Runnels
!!    August 10, 2007
!!===================================================================================
#include "f90_assert.fpp"
module truchas_timing

use parallel_info_module, only : p_info
use timing_tree
use truchas_logging_services

implicit none

private
public :: report_tree_timing

contains

!!===================================================================================
!!REPORT TREE TIMING
!! 
!!PURPOSE(S):
!! 
!!    Manage the output printing for truchas. Prints output in 2 ways:
!!    1. Print XML output using TBrook. If using multiple processors, all of the
!!       processors send their data to the IO processor, and it prints a timing tree
!!       for each one.
!!    2. Print the average times for each timer and each processor's times to the
!!       command line. If only using one processor, only one appears.
!!    3. Do the same as #2, but print to a file specified by unit.
!!
!!===================================================================================
subroutine report_tree_timing()
  use pgslib_module

  integer             :: temp
  integer             :: i
  integer             :: j
  integer             :: depth
  integer             :: num
  integer             :: istat
  character(len=70)   :: attributes
  character(len=75)   :: formatstring
  character(len=175)   :: output_string

  integer,                    pointer     :: order(:)
  character(len=MAX_KEY_LEN), pointer     :: keys(:)
  real,                       pointer     :: walls(:)
  real,                       pointer     :: cpus(:)

  integer,                    allocatable :: order_grid(:, :)
  character(len=MAX_KEY_LEN), allocatable :: key_grid(:,:)
  real,                       allocatable :: wall_grid(:, :)
  real,                       allocatable :: cpu_grid(:, :)

  character(len=MAX_KEY_LEN) :: key
  real              :: wall
  real              :: cpu
  real              :: cpu_max
  real              :: cpu_min

  temp = 0
  num = 0
  istat = 0
  
  call flatten_timing_tree(order, keys, walls, cpus)

  num = size(keys)
  call pgslib_bcast(num)
  ASSERT(pgslib_global_all(num==size(keys)))

  allocate(order_grid(p_info%npe, num*2))
  allocate(key_grid(p_info%npe, num))
  allocate(wall_grid(p_info%npe, num))
  allocate(cpu_grid(p_info%npe, num))
  
  !!Store all the timing data for all processors into four 2 dimensional arrays

  do i = 1, size(order)
     call pgslib_collate(order_grid(:, i), order(i))
     if (p_info%iop) then
        ASSERT(all(order_grid(:, i) == order(i)))
     end if
  end do

  do i = 1, size(keys)
     call pgslib_collate(key_grid(:, i), keys(i))
     call pgslib_collate(wall_grid(:, i), walls(i))
     call pgslib_collate(cpu_grid(:, i), cpus(i))
     if (p_info%iop) then 
        ASSERT(all(key_grid(:, i) == keys(i)))
     end if
  end do

  !!Code to be executed by IO processor
  if (p_info%iop) then

     !!PRINT HEADERS
     call TLS_info ('')
     call TLS_info ('  TIMING SUMMARY                       AVERAGE        MIN          MAX   ')
     call TLS_info ('  --------------                      ---------    ---------    ---------')

     depth = 0
     temp = 0
     do i = 1, num*2

        key = key_grid(1, order_grid(1, i))
        cpu = 0
        cpu_max = cpu_grid(1, order_grid(1, i))
        cpu_min = cpu_grid(1, order_grid(1, i))

        do j = 1, p_info%npe
           cpu = cpu + cpu_grid(j, order_grid(1, i))
           if (cpu_max < cpu_grid(j, order_grid(1, i))) cpu_max = cpu_grid(j, order_grid(1, i))
           if (cpu_min > cpu_grid(j, order_grid(1, i))) cpu_min = cpu_grid(j, order_grid(1, i))
        end do
     
        cpu = cpu/p_info%npe

        if (temp < order_grid(1, i)) then
           depth = depth + 1

           if (.not.(depth == 0)) then
              write(formatString, "(a, i2, a)") "(", depth*2, "x, A, A, ES13.3, ES13.3, ES13.3)"
              write(output_string, fmt = formatString) key, " ", cpu, cpu_min, cpu_max
              call TLS_info (output_string)

           else
              write(formatString, "(a, i2, a)") "(A, A, ES13.3, ES13.3, ES13.3)"
              write(output_string, fmt = formatString) key, " ", cpu, cpu_min, cpu_max
              call TLS_info (output_string)
           end if

        else
           depth = depth - 1
        end if
 
        temp = order_grid(1, i)

     end do

  end if

  call TLS_info ('')

  deallocate(order_grid, key_grid, wall_grid, cpu_grid)

  call destroy_timing_tree()

end subroutine report_tree_timing

end module truchas_timing
