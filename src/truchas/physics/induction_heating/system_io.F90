!  Rudimentary (evolving) system-dependent I/O module, for IBM XLF Version 3.2
!
!  Neil N. Carlson, June 1996

module system_io

  public :: new_unit

  ! IOSTAT values... (there are a great many more!)
  integer, parameter, public :: END_OF_FILE = -1
  integer, parameter, public :: END_OF_INTERNAL_FILE = -2
  integer, parameter, public :: END_OF_RECORD = -4

  ! Default preconnected units...
  integer, parameter, public :: STDERR = 0
  integer, parameter, public :: STDIN  = 5
  integer, parameter, public :: STDOUT = 6
  integer, dimension(3), parameter, public :: PRECONNECTED_UNITS = (/ STDERR, STDIN, STDOUT /)

  ! Largest possible unit number.
  integer, parameter, public :: MAX_UNIT_NUMBER = 100

  contains

   !!!
   !!!  Find an unconnected unit number.

    subroutine new_unit (unit)

      integer, intent(out) :: unit
      
      logical :: exists, opened
      integer :: ios, j

      do j = 1, MAX_UNIT_NUMBER

        if (any (j == PRECONNECTED_UNITS)) then
          cycle
        end if

        inquire (unit=j, exist=exists, opened=opened, iostat=ios)
        if (exists .and. .not. opened .and. ios == 0) then
          unit = j
          return
        end if

      end do

      unit = -1

    end subroutine new_unit

end module system_io
