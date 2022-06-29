!!
!! STRING_UTILITIES
!!
!! Neil N. Carlson <nnc@lanl.gov>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This module provides a collection of useful character string-related
!! utilities:
!!
!!    RAISE_CASE(CS) returns a string of the same length and value as CS but
!!      with each lowercase letter replaced by its uppercase counterpart.
!!
!!    LOWER_CASE(CS) returns a string of the same length and value as CS but
!!      with each uppercase letter replaced by its lowercase counterpart.
!!
!!    I_TO_C(N) returns the printed form of the integer N as a string.  The
!!      length of the string is precisely the number of characters required
!!      to express the integer without leading or trailing blanks.
!!

module string_utilities

  implicit none
  private

  public :: raise_case, lower_case, i_to_c

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! RAISE_CASE and LOWER_CASE
 !!
 !!  Functions for mapping strings to upper and lower case.
 !!

  elemental function raise_case (cs)
    character(len=*), intent(in) :: cs
    character(len=len(cs)) :: raise_case
    integer :: i
    do i = 1, len(cs)
      if (iachar(cs(i:i)) >= iachar('a') .and. iachar(cs(i:i)) <= iachar('z')) then
        raise_case(i:i) = achar(iachar(cs(i:i)) - iachar('a') + iachar('A'))
      else
        raise_case(i:i) = cs(i:i)
      end if
    end do
  end function raise_case

  elemental function lower_case (cs)
    character(len=*), intent(in) :: cs
    character(len=len(cs)) :: lower_case
    integer :: i
    do i = 1, len(cs)
      if (iachar(cs(i:i)) >= iachar('A') .and. iachar(cs(i:i)) <= iachar('Z')) then
        lower_case(i:i) = achar(iachar(cs(i:i)) - iachar('A') + iachar('a'))
      else
        lower_case(i:i) = cs(i:i)
      end if
    end do
  end function lower_case

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! I_TO_C -- INTEGER-TO-CHARACTER CONVERSION.
 !!
 !! Function takes an integer and returns its printed form as a string.
 !! The length of the character result is precisely the number of characters
 !! required to express the integer; there are no leading or trailing blanks.
 !! The routine OUTPUT_LEN which follows is a private function that returns
 !! the required length.
 !!
 !! These routines were part of the original EM solver source code imported
 !! into Truchas, and have been moved here for general use.
 !!
 !! Neil N. Carlson, Carlson Science Computing <nnc@newmexico.com>
 !!

  pure function i_to_c (n)result(s)

    integer, intent(in) :: n
    character(len=output_length(n)) :: s

    integer :: absn
    integer :: pos, digit

    absn = abs(n)
    if (absn < 0) then ! Can't treat -huge(0)-1 properly
      s = repeat('*', len(s))
      return
    end if

    if (n < 0) then
      s(1:1) = '-'
    else if (n == 0) then
      s = '0'
      return
    end if

    pos = len(s)
    do while (absn > 0)
      digit = mod(absn, 10) + 1
      s(pos:pos) = '0123456789'(digit:digit)
      pos = pos - 1
      absn = absn / 10
    end do

  end function i_to_c

  pure function output_length(n) result(l)

    integer, intent(in) :: n
    integer :: l

    select case (n)
    case (:-1)  ! negative
      l = int(log10(-dble(n)+0.1d0) + 2)
    case (0)
      l = 1
    case (1:)   ! positive
      l = int(log10(n+0.1d0) + 1)
    end select

  end function output_length

end module string_utilities

