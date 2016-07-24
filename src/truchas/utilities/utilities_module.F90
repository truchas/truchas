!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE UTILITIES_MODULE
  !======================================================================
  ! Purpose(s):
  !
  !   Define general utility routines.
  !
  !   Public Interface(s):
  !
  !     * call STRING_COMPARE (string1, string2, string_match)
  !
  !       Performs a string comparison.
  !
  !     * logical = FSTRING_COMPARE (string1, string2)
  !
  !       Performs a string comparison.
  !
  !     * call TIMESTAMP (date_time)
  !
  !       Returns the date and time in string date_time.
  !
  !     * call MAKE_DIRECTORY (string, status)
  !
  !       an F90 interface to mkdir
  !
  ! Contains: STRING_COMPARE
  !           TIMESTAMP
  !           MAKE_DIRECTORY
  !           MAKE_DIRECTORY_HIERARCHY
  !
  ! Author(s): Douglas B. Kothe, LANL T-3 (dbk@lanl.gov)
  !            Anand V. Reddy, Cateripillar (reddy_anand_v@cat.com)
  !            Sriram Swaminarayan (sriram@lanl.gov)
  !            Neil N. Carlson <nnc@newmexico.com>
  !=======================================================================
  implicit none
  private

  ! public procedures
  public :: STRING_COMPARE,  &
            FSTRING_COMPARE, &
            TIMESTAMP,       &
            MAKE_DIRECTORY,  &
            MAKE_DIRECTORY_HIERARCHY,  &
            BUBBLE_PERMUTE

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

     interface
        subroutine MAKE_DIRECTORY_C(path, status) bind(c, name='make_directory_c')
           use,intrinsic :: iso_c_binding, only: c_char, c_int
           character(kind=c_char), intent(in)  :: path(*)
           integer(c_int), intent(out) :: status
        end subroutine
     end interface

     interface
        subroutine MAKE_DIRECTORY_HIER_C(path, status) bind(c,name='make_directory_hier_c')
           use,intrinsic :: iso_c_binding, only: c_char, c_int
           character(kind=c_char), intent(in) :: path(*)
           integer(c_int), intent(out) :: status
        end subroutine
     end interface

CONTAINS


   function fstring_compare (string1, string2)
      !-------------------------------------------------------------------------
      ! purpose:
      !
      !    case insensitive string comparision - a function interface to a
      !    subroutine
      !-------------------------------------------------------------------------

      ! argument list
      character(*), intent(IN) :: string1
      character(*), intent(IN) :: string2

      ! local variables
      logical :: status

      ! function return
      logical :: fstring_compare

      !-------------------------------------------------------------------------

      call STRING_COMPARE (string1, string2, status)

      fstring_compare = status

   end function fstring_compare

  SUBROUTINE STRING_COMPARE (string1, string2, strings_match)
    !=======================================================================
    ! Purpose(s):
    !
    !   Performs a case insensitive comparison of two character strings,
    !   and returns true or false for the comparison result.
    !
    !=======================================================================
    use string_utilities, only: lower_case

    character(*), intent(in) :: string1, string2
    logical, intent(out) :: strings_match

    strings_match = (lower_case(trim(string1)) == lower_case(trim(string2)))

  END SUBROUTINE STRING_COMPARE

  SUBROUTINE TIMESTAMP (date_time)
    !=======================================================================
    ! Purpose(s):
    !
    !   Subroutine to build a 26-character string containing current
    !   date and time using DATE_AND_TIME intrinsic.
    !
    !                                    12345678901234567890123456
    !   String returned is of the form:  Fri 20 Aug 93 09:33:35.500
    !
    !   Initial routine obtained from ftp.ora.com:/pub/nutshell/fortran90
    !   as part of the examples for the O'Reilly book "Migrating to
    !   Fortran 90" by James F. Kerrigan.
    !
    !=======================================================================

    ! Argument List
    character(*), intent(OUT) :: date_time

    ! Local Variables
    character(3), dimension(0:6) :: days
    character(3), dimension(12)  :: months

    integer, dimension(8) :: elements
    integer :: m, y, w

    ! Define Days and Months
    data days   / 'Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat' /  
    data months / 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', &
                  'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'  /

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! System date and time
    call DATE_AND_TIME (VALUES = elements)

    ! Format date and time.
    INVALID: if (elements(1) /= -HUGE(0)) then
       y = elements(1)
       m = elements(2)

       if (m < 3) then
          m = m + 12
          y = y -  1
       end if

       w = MOD((elements(3) + (13*m-27)/5 + y + y/4 - y/100 + y/400 ), 7)  
       CENTURY: if (elements(1) < 2000 ) then
          elements(1) = elements(1) - 1900
       else
          elements(1) = elements(1) - 2000
       end if CENTURY

       date_time = ' '
       write (date_time,10) days(w), elements(3), months (elements(2)), &
                            elements(1), elements(5), elements(6), &
                            elements(7), elements(8)
10     format (a3, 1x, i2.2, 1x, a3, 1x, i2.2, 1x,   &
            i2.2, ':', i2.2, ':', i2.2, '.', i3.3 ) 
    else INVALID
       date_time = ' '
    end if INVALID

  END SUBROUTINE TIMESTAMP

  !-----------------------------------------------------------------------------

  function MAKE_DIRECTORY (path)

     use parallel_info_module, only: p_info

     character(*) :: path
     integer :: status
     character(1024) :: path2
     integer :: MAKE_DIRECTORY

     if (LEN_TRIM(path) > 1023) then
        MAKE_DIRECTORY = -1
        return
     end if

     ! copy string to temporary string
     path2 = TRIM(ADJUSTL(path))
     ! null terminate the temporary string
     path2(LEN_TRIM(path2)+1:) = ACHAR(0)

     if (p_info%iop) then
        call MAKE_DIRECTORY_C (path2, status)
     else
        status = 0
     end if

     ! return status
     MAKE_DIRECTORY = status

  end function MAKE_DIRECTORY

  function MAKE_DIRECTORY_HIERARCHY (path)

     use parallel_info_module, only: p_info

     character(*) :: path
     integer :: status
     character(1024) :: path2
     integer :: MAKE_DIRECTORY_HIERARCHY

     if (LEN_TRIM(path) > 1023) then
        MAKE_DIRECTORY_HIERARCHY = -1
        return
     end if

     ! copy string to temporary string
     path2 = TRIM(ADJUSTL(path))
     ! null terminate the temporary string
     path2(LEN_TRIM(path2)+1:) = ACHAR(0)

     if (p_info%iop) then
        call MAKE_DIRECTORY_HIER_C (path2, status)
     else
        status = 0
     end if

     ! return status
     MAKE_DIRECTORY_HIERARCHY = status

  end function MAKE_DIRECTORY_HIERARCHY

  !-----------------------------------------------------------------------------

  function BUBBLE_PERMUTE (vector)
     !--------------------------------------------------------------------------
     ! Purpose:
     !    return a permutation vector that will permute an unsorted
     !    integer vector into a vector sorted into ascending order
     !
     !    vector(permute(i)) is the i'th element of the sorted list
     !
     !    warning: this is a simple n^2 sort, and should not be used
     !    if n is greater than about 10
     !
     !--------------------------------------------------------------------------

     ! arguments
     integer, dimension(:) :: Vector

     ! function return
     integer, dimension(SIZE(Vector)) :: BUBBLE_PERMUTE

     ! local variables
     integer :: i, j, n, tmp

     !--------------------------------------------------------------------------

     n = SIZE(Vector)

     ! initialize permutation vector to do nothing
     do i = 1, n
        BUBBLE_PERMUTE(i) = i
     end do

     ! bubble sort on the vector, reordering the permutation indices
     do j = n-1, 1, -1
        do i = 1, j
           if (Vector(BUBBLE_PERMUTE(i)) > Vector(BUBBLE_PERMUTE(i+1))) then
              ! exchange the permutation indices
              tmp               = BUBBLE_PERMUTE(i+1)
              BUBBLE_PERMUTE(i+1) = BUBBLE_PERMUTE(i)
              BUBBLE_PERMUTE(i)   = tmp
           end if
        end do
     end do

  end function BUBBLE_PERMUTE

  !-----------------------------------------------------------------------------

END MODULE UTILITIES_MODULE
