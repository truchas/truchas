!!
!! STRING_UTILITIES
!!
!! Neil N. Carlson <nnc@lanl.gov>
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
!!    CTRIM(STRING) returns STRING with all trailing blanks and null
!!      terminator removed.  Use this to strip off any trailing null
!!      character that a simple TRIM would leave; C-style strings are NT.
!!
!!    CONCAT(ARRAY) returns a atring that is the concatenation of the strings
!!      in ARRAY.  Leading and trailing blanks in each of the string elements
!!      is removed before concatenation, and a single blank character is
!!      interposed between successive strings.  String elements that are
!!      blank are ignored.
!!

module string_utilities

  implicit none
  private
  
  public :: raise_case, lower_case, i_to_c, ctrim, concat
  
contains
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! RAISE_CASE and LOWER_CASE
 !!
 !!  Functions for mapping strings to upper and lower case.
 !!  
 !!  Moved here from EM_input.  NNC
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
  
  pure function i_to_c (n) result (s)
  
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
  
  pure function output_length (n) result (l)
  
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
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! CTRIM -- trim C-style strings
 !!
 !! This function behaves just like the intrinsic TRIM, except that it goes
 !! one step further and removes any trailing sequence of null characters that
 !! a simple TRIM would leave behind.  This is useful when a null-terminated
 !! C string has been copied (including null characters) into a Fortran string.
 !!
 
  function ctrim (string) result (s)
    character(len=*), intent(in) :: string
    character(len=len_ctrim(string)) :: s
    s = string(:len_ctrim(string))
  end function ctrim
  
  pure integer function len_ctrim (string)
    character(len=*), intent(in) :: string
    len_ctrim = len_trim(string)
    if (len_ctrim == 0) return
    do while (string(len_ctrim:len_ctrim) == char(0))
      len_ctrim = len_ctrim - 1
      if (len_ctrim < 1) exit
    end do
  end function len_ctrim

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! CONCAT -- concatenate an array of strings into a single string
 !!
 !! This function returns a string that is the concatenation of the strings
 !! in ARRAY.  Leading and trailing blanks in each of the string elements is
 !! removed before concatenation, and a single blank character is interposed
 !! between successive strings.  String elements that are blank are ignored.
 !!

  pure function concat (array) result (s)
    character(len=*), intent(in) :: array(:)
#ifdef PATHSCALE_COMPILER_WORKAROUND
    character(len=size(array)*len(array)) :: s
#else
    character(len=len_concat(array)) :: s
#endif
    character(len(array)) :: tmp
    integer :: j, m, n
    m = 0
    s = ''
    do j = 1, size(array)
      tmp = adjustl(array(j))
      n = len_trim(tmp)
      if (n == 0) cycle
      if (m > 0) m = m + 1
      s(m+1:m+n) = trim(tmp)
      m = m + n
    end do
  end function concat

#ifndef PATHSCALE_COMPILER_WORKAROUND
  pure integer function len_concat (array) result (m)
    character(len=*), intent(in) :: array(:)
    integer :: j, n
    m = 0
    do j = 1, size(array)
      n = len_trim(adjustl(array(j)))
      if (n == 0) cycle
      if (m > 0) m = m + 1
      m = m + n
    end do
  end function len_concat
#endif

end module string_utilities

