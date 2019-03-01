!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main

  !use,intrinsic :: ieee_exceptions
  use drivers, only: code
  implicit none

  !call ieee_set_halting_mode ([ieee_invalid, ieee_overflow, ieee_divide_by_zero], .true.)
  call code ()

end program main
