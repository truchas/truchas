!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module test_rade_tools_common

  use,intrinsic :: iso_fortran_env, only: int32
  implicit none

contains

  !! Randomly permutes an array using th Fisher-Yates algorithm
  !! Source: https://www.rosettacode.org/wiki/Knuth_shuffle#Fortran
  subroutine shuffle(a)

    integer, intent(inout) :: a(:)

    integer :: i, randpos, temp
    real :: r

    do i = size(a), 2, -1
      call random_number(r)
      randpos = int(r * i) + 1
      temp = a(randpos)
      a(randpos) = a(i)
      a(i) = temp
    end do

  end subroutine shuffle


  !! Returns a random exactly-representable real number.
  !! Rationals of the form a/2^k, where k is less than the
  !! number of the exponent bits, are exactly representatble.
  elemental impure subroutine random_exact(rand)

    real, intent(out) :: rand

    !! We only want fractions of the form 1/2^k
    integer(int32), parameter :: mask = int(Z'3F800000',int32)
    integer(int32) :: i
    real :: r

    call random_number(r)
    i = transfer(r, i)
    i = iand(i, mask)
    rand = transfer(i, r)

  end subroutine random_exact

end module test_rade_tools_common
