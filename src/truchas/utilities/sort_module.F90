!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE SORT_MODULE
  !=======================================================================
  ! Purpose(s):
  !   Define various procedures for sorting arrays. In general, these
  !   routines were adapted from the text Numerical Recipes in Fortran,
  !   2nd Edition, Volume 2, Fortran 90.
  !
  ! Interface(s):
  !
  !   * SORT(Array_Test, Array_Index, sort_method)
  !
  !     Sort one (real) or two (real and integer) arrays into ascending
  !     numerical order by either the Heapsort method, or by using straight
  !     insertion, or by Shell's method. The sorting is based on the
  !     values of the real array, Array_Test. The elements of the integer
  !     array, Array_Index, are sorted in the same manner as the real array.
  !
  !     Array_Test  - Real rank one array to be sorted. This array is
  !                   modified upon exit of the routine.
  !
  !     Array_Index - Optional interger rank one sorting index array. This
  !                   array can be used to make a sorted pointing to an
  !                   origonal, non-sorted, copy of Array_Test, i.e.
  !                   Array_Test(Array_Index(i)).
  !
  !     sort_method - Optional symbolic constant integer used to specify
  !                   the type of sorting method. Use either heap_sort,
  !                   insert_sort, or shell_sort. If sort_method is not
  !                   supplied then a default method, which is set in the
  !                   SORT routine, is selected.
  !
  !     Note: SORT is designed for single processor applications. Therefore,
  !     it is not appropriate for multi-processor applications. SORT uses the
  !     F90 array syntax. Thus for very large arrays it may require repeated
  !     cache flushes and is therefore intended for relatively small arrays,
  !     i.e. Array_Test should not be on the order of ncells.
  !
  !
  ! Contains: SORT
  !
  !           SORT1_HEAP
  !           SORT2_HEAP
  !
  !           SORT1_INSERT
  !           SORT2_INSERT
  !
  !           SORT1_SHELL
  !           SORT2_SHELL
  !
  !           ----------
  !
  !           SIFT1_DOWN
  !           SIFT2_DOWN
  !
  !           REVERSE_INTEGER
  !           REVERSE_REAL
  !
  !           SWAP_INTEGER
  !           SWAP_INTEGER_R1
  !           SWAP_REAL
  !           SWAP_REAL_R1
  !
  ! Author(s): Jerry S. Brock, LANL T-3 (jsbrock@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  implicit none
  private

  public :: SORT, REVERSE

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Generic Procedure Interfaces
  INTERFACE REVERSE
     MODULE PROCEDURE REVERSE_INTEGER
     MODULE PROCEDURE REVERSE_REAL
  END INTERFACE

  INTERFACE SWAP
     MODULE PROCEDURE SWAP_INTEGER
     MODULE PROCEDURE SWAP_INTEGER_R1
     MODULE PROCEDURE SWAP_REAL
     MODULE PROCEDURE SWAP_REAL_R1
  END INTERFACE

  ! Sort Technique Parameters
  integer, parameter, public :: heap_sort   = 1
  integer, parameter, public :: insert_sort = 2
  integer, parameter, public :: shell_sort  = 3

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE SORT (Array_Test, Array_Index, sort_method)
    !=======================================================================
    ! Purpose(s):
    !   Sort one (real) or two (real and integer) arrays into ascending
    !   numerical order by either the Heapsort method, or by using straight
    !   insertion, or by Shell's method. The sorting is based on the
    !   values of the real array, Array_Test. The elements of the integer
    !   array, Array_Index, are sorted in the same manner as the real array.
    !=======================================================================

    ! Argument List
    integer, intent(IN), optional :: sort_method

    integer, dimension(:), intent(INOUT), optional :: Array_Index
    real(r8), dimension(:), intent(INOUT) :: Array_Test

    ! Local Variables
    logical :: one_array
    integer :: sort_method_flag
    integer :: sort_method_default = insert_sort

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Sort Method
    if (PRESENT(sort_method)) then
       sort_method_flag = sort_method         ! Input Method
    else
       sort_method_flag = sort_method_default ! Default Method
    end if

    ! Number of Arrays
    one_array = .not. PRESENT(Array_Index)

    ! Sort Array(s)
    select case (sort_method_flag)
       ! Heapsort Method
       case (heap_sort)
          if (one_array) then
             call SORT1_HEAP (Array_Test)
          else
             call SORT2_HEAP (Array_Test, Array_Index)
          end if

       ! Straight Insertion
       case (insert_sort)
          if (one_array) then
             call SORT1_INSERT (Array_Test)
          else
             call SORT2_INSERT (Array_Test, Array_Index)
          end if

       ! Shell's Method
       case (shell_sort)
          if (one_array) then
             call SORT1_SHELL (Array_Test)
          else
             call SORT2_SHELL (Array_Test, Array_Index)
          end if
    end select

  END SUBROUTINE SORT

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SORT1_HEAP (Array_Test)
    !=======================================================================
    ! Purpose(s):
    !   Sort a real array into ascending numerical order by the Heapsort
    !   method. This is an N*log(N) routine; its worst case run time is
    !   only 20 percent longer than its average runtime.
    !
    !   This routine was adapted from the text Numerical Recipes in Fortran,
    !   2nd Edition, Volume 2, Fortran 90, page 1171.
    !=======================================================================

    ! Argument List
    real(r8), dimension(:), intent(INOUT) :: Array_Test

    ! Local Variables
    integer :: i
    integer :: num

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Array Size
    num = SIZE(Array_Test)

    ! Heap Creation
    do i = num/2, 1, -1
       call SIFT1_DOWN (Array_Test, i, num)
    end do

    ! Heap Selection
    do i = num, 2, -1
       call SWAP (Array_Test(1), Array_Test(i))
       call SIFT1_DOWN (Array_Test, 1, i-1)
    end do

  END SUBROUTINE SORT1_HEAP

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SORT2_HEAP (Array_Test, Array_Index)
    !=======================================================================
    ! Purpose(s):
    !   Sort two arrays into ascending numerical order by the Heapsort method.
    !   The sorting is based on the values of the real array, Array_Test. The
    !   elements of the integer array, Array_Index, are sorted in the same
    !   manner as the real array. This is an N*log(N) routine; its worst case
    !   run time is only 20 percent longer than its average runtime.
    !
    !   This routine was adapted from the text Numerical Recipes in Fortran,
    !   2nd Edition, Volume 2, Fortran 90, page 1171.
    !=======================================================================

    ! Argument List
    integer,  dimension(:), intent(INOUT) :: Array_Index
    real(r8), dimension(:), intent(INOUT) :: Array_Test

    ! Local Variables
    integer :: i
    integer :: num

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Array Size
    num = SIZE(Array_Test)

    ! Heap Creation
    do i = num/2, 1, -1
       call SIFT2_DOWN (Array_Test, Array_Index, i, num)
    end do

    ! Heap Selection
    do i = num, 2, -1
       call SWAP (Array_Test (1), Array_Test (i))
       call SWAP (Array_Index(1), Array_Index(i))
       call SIFT2_DOWN (Array_Test, Array_Index, 1, i-1)
    end do

  END SUBROUTINE SORT2_HEAP

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SORT1_INSERT (Array_Test)
    !=======================================================================
    ! Purpose(s):
    !   Sort a real array into ascending numerical order by straight
    !   insertion. This is an N**2 routine; it should only be used for
    !   relatively small sorting tasks, i.e. num < 20.
    !
    !   This routine was adapted from the text Numerical Recipes in Fortran,
    !   2nd Edition, Volume 2, Fortran 90, page 1167.
    !=======================================================================

    ! Argument List
    real(r8), dimension(:), intent(INOUT) :: Array_Test

    ! Local Variables
    integer :: i, j
    integer :: num
    real(r8) :: tmp_test

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Array Size
    num = SIZE(Array_Test)

    ! Sort Array
    do j = 2, num
       ! Temporary Value
       tmp_test = Array_Test(j)

       ! Check Array Values
       do i = j-1, 1, -1
          ! Check Value
          if (Array_Test(i) <= tmp_test) exit

          ! Switch Values
          Array_Test(i+1) = Array_Test(i)
       end do

       ! Insert Temporary
       Array_Test(i+1) = tmp_test
    end do

  END SUBROUTINE SORT1_INSERT

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SORT2_INSERT (Array_Test, Array_Index)
    !=======================================================================
    ! Purpose(s):
    !   Sort two arrays into ascending numerical order by straight insertion.
    !   The sorting is based on the values of the real array, Array_Test. The
    !   elements of the integer array, Array_Index, are sorted in the same
    !   manner as the real array. This is an N**2 routine; it should only be
    !   used for relatively small sorting tasks, i.e. num < 20.
    !
    !   This routine was adapted from the text Numerical Recipes in Fortran,
    !   2nd Edition, Volume 2, Fortran 90, page 1167.
    !=======================================================================

    ! Argument List
    integer,  dimension(:), intent(INOUT) :: Array_Index
    real(r8), dimension(:), intent(INOUT) :: Array_Test

    ! Local Variables
    integer :: i, j
    integer :: num
    integer :: tmp_index
    real(r8) :: tmp_test

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Array Size
    num = SIZE(Array_Test)

    ! Sort Arrays
    do j = 2, num
       ! Temporary Values
       tmp_index = Array_Index(j)
       tmp_test  = Array_Test (j)

       ! Check Array Values
       do i = j-1, 1, -1
          ! Check Value
          if (Array_Test(i) <= tmp_test) exit

          ! Switch Values
          Array_Index(i+1) = Array_Index(i)
          Array_Test (i+1) = Array_Test (i)
       end do

       ! Insert Temporarys
       Array_Index(i+1) = tmp_index
       Array_Test (i+1) = tmp_test
    end do

  END SUBROUTINE SORT2_INSERT

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SORT1_SHELL (Array_Test)
    !=======================================================================
    ! Purpose(s):
    !   Sort a real array into ascending numerical order using Shell's
    !   method. This is approximately an N**1.5 routine; it is recommended
    !   for relatively moderate sorting tasks.
    !
    !   This routine was adapted from the text Numerical Recipes in Fortran,
    !   2nd Edition, Volume 2, Fortran 90, page 1168.
    !=======================================================================

    ! Argument List
    real(r8), dimension(:), intent(INOUT) :: Array_Test

    ! Local Variables
    integer :: i, inc, j
    integer :: num
    real(r8) :: tmp_test

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Array Size
    num = SIZE(Array_Test)

    ! Starting Increment
    inc = 1
    do
       inc = 3 * inc + 1
       if (inc > num) exit
    end do

    ! Straight Insertion
    do
       ! Partial Sorts
       inc = inc / 3

       ! Insertion Outer Loop
       do i = inc + 1, num
          ! Temporay Values
          j = i
          tmp_test = Array_Test(i)

          ! Insertion Inner Loop
          do
             ! Check Values
             if (Array_Test(j - inc) <= tmp_test) exit

             ! Switch Values
             Array_Test(j) = Array_Test(j - inc)

             ! Inner Loop Index
             j = j - inc
             if (j <= inc) exit
          end do

          ! Insert Temporarys
          Array_Test(j) = tmp_test
       end do

       ! End Partial Sort
       if (inc <= 1) exit
    end do

  END SUBROUTINE SORT1_SHELL

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SORT2_SHELL (Array_Test, Array_Index)
    !=======================================================================
    ! Purpose(s):
    !   Sort two arrays into ascending numerical order using Shell's method.
    !   The sorting is based on the values of the real array, Array_Test. The
    !   elements of the integer array, Array_Index, are sorted in the same
    !   manner as the real array. This is approximately an N**1.5 routine; it
    !   is recommended for relatively moderate sorting tasks.
    !
    !   This routine was adapted from the text Numerical Recipes in Fortran,
    !   2nd Edition, Volume 2, Fortran 90, page 1168.
    !=======================================================================

    ! Argument List
    integer,  dimension(:), intent(INOUT) :: Array_Index
    real(r8), dimension(:), intent(INOUT) :: Array_Test

    ! Local Variables
    integer :: i, inc, j
    integer :: num
    integer :: tmp_index
    real(r8) :: tmp_test

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Array Size
    num = SIZE(Array_Test)

    ! Starting Increment
    inc = 1
    do
       inc = 3 * inc + 1
       if (inc > num) exit
    end do

    ! Straight Insertion
    do
       ! Partial Sorts
       inc = inc / 3

       ! Insertion Outer Loop
       do i = inc + 1, num
          ! Temporay Values
          j = i
          tmp_index = Array_Index(i)
          tmp_test  = Array_Test (i)

          ! Insertion Inner Loop
          do
             ! Check Values
             if (Array_Test(j - inc) <= tmp_test) exit

             ! Switch Values
             Array_Index(j) = Array_Index(j - inc)
             Array_Test (j) = Array_Test (j - inc)

             ! Inner Loop Index
             j = j - inc
             if (j <= inc) exit
          end do

          ! Insert Temporarys
          Array_Index(j) = tmp_index
          Array_Test (j) = tmp_test
       end do

       ! End Partial Sort
       if (inc <= 1) exit
    end do

  END SUBROUTINE SORT2_SHELL

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SIFT1_DOWN (Array_Test, left, right)
    !=======================================================================
    ! Purpose(s):
    !   Sift one (real) array for the Heapsort method.
    !=======================================================================

    ! Argument List
    integer, intent(IN) :: left, right
    real(r8), dimension(:), intent(INOUT) :: Array_Test

    ! Local Variables
    integer :: j, jold
    real(r8) :: tmp_test

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize Variables
    tmp_test = Array_Test(left)

    ! Initialize Indcies
    jold = left
    j    = left + left

    ! Sift Down Value
    do
       if (j > right) exit ! Do While j <= right
       if (j < right) then ! Compare to Underling
          if (Array_Test(j) < Array_Test(j+1)) j = j + 1
       end if

       ! Test Temporay
       if (tmp_test >= Array_Test(j)) exit

       ! Set Value
       Array_Test(jold) = Array_Test(j)

       ! Array Indcies
       jold = j
       j    = j + j
    end do

    ! Replace Temporary
    Array_Test(jold) = tmp_test

  END SUBROUTINE SIFT1_DOWN

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SIFT2_DOWN (Array_Test, Array_Index, left, right)
    !=======================================================================
    ! Purpose(s):
    !   Sift two (real and integer) arrays for the Heapsort method. The
    !   sifting is based on the values of the real array, Array_Test. The
    !   elements of the integer array, Array_Index, are sifted in the same
    !   manner as the real array.
    !=======================================================================

    ! Argument List
    integer, intent(IN) :: left, right
    integer,  dimension(:), intent(INOUT) :: Array_Index
    real(r8), dimension(:), intent(INOUT) :: Array_Test

    ! Local Variables
    integer :: j, jold
    integer :: tmp_index
    real(r8)   :: tmp_test

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize Variables
    tmp_test  = Array_Test (left)
    tmp_index = Array_Index(left)

    ! Initialize Indcies
    jold = left
    j    = left + left

    ! Sift Down Value
    do
       if (j > right) exit ! Do While j <= right
       if (j < right) then ! Compare to Underling
          if (Array_Test(j) < Array_Test(j+1)) j = j + 1
       end if

       ! Test Temporay
       if (tmp_test >= Array_Test(j)) exit

       ! Set Value
       Array_Test (jold) = Array_Test (j)
       Array_Index(jold) = Array_Index(j)

       ! Array Indcies
       jold = j
       j    = j + j
    end do

    ! Replace Temporary
    Array_Test (jold) = tmp_test
    Array_Index(jold) = tmp_index

  END SUBROUTINE SIFT2_DOWN

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE REVERSE_INTEGER (A)
    !=======================================================================
    ! Purpose(s):
    !   Reverse the order of an integer array.
    !=======================================================================

    ! Argument List
    integer, dimension(:), intent(INOUT) :: A

    ! Local Variables
    integer :: n
    integer :: num
    integer, dimension(SIZE(A)) :: Tmp

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Array Size
    num = SIZE(A)

    ! Reverse Array
    do n = 1, num
       Tmp(n) = A(num - n + 1)
    end do
    A = Tmp

  END SUBROUTINE REVERSE_INTEGER

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE REVERSE_REAL (A)
    !=======================================================================
    ! Purpose(s):
    !   Reverse the order of a real array.
    !=======================================================================

    ! Argument List
    real(r8), dimension(:), intent(INOUT) :: A

    ! Local Variables
    integer :: n
    integer :: num
    real(r8), dimension(SIZE(A)) :: Tmp

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Array Size
    num = SIZE(A)

    ! Reverse Array
    do n = 1, num
       Tmp(n) = A(num - n + 1)
    end do
    A = Tmp

  END SUBROUTINE REVERSE_REAL

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SWAP_INTEGER (a, b)
    !=======================================================================
    ! Purpose(s):
    !   Swap two integer values.
    !=======================================================================

    ! Argument List
    integer, intent(INOUT) :: a
    integer, intent(INOUT) :: b

    ! Local Variables
    integer :: tmp

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Swap Values
    tmp = a
    a = b
    b = tmp

  END SUBROUTINE SWAP_INTEGER

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SWAP_INTEGER_R1 (A, B)
    !=======================================================================
    ! Purpose(s):
    !   Swap two rank one integer arrays.
    !=======================================================================

    ! Argument List
    integer, dimension(:), intent(INOUT) :: A
    integer, dimension(:), intent(INOUT) :: B

    ! Local Variables
    integer, dimension(SIZE(A)) :: Tmp

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Swap Arrays
    Tmp = A
    A = B
    B = Tmp

  END SUBROUTINE SWAP_INTEGER_R1

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SWAP_REAL (a, b)
    !=======================================================================
    ! Purpose(s):
    !   Swap two real values.
    !=======================================================================

    ! Argument List
    real(r8), intent(INOUT) :: a
    real(r8), intent(INOUT) :: b

    ! Local Variables
    real(r8) :: tmp

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Swap Values
    tmp = a
    a = b
    b = tmp

  END SUBROUTINE SWAP_REAL

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE SWAP_REAL_R1 (A, B)
    !=======================================================================
    ! Purpose(s):
    !   Swap two rank one real arrays.
    !=======================================================================

    ! Argument List
    real(r8), dimension(:), intent(INOUT) :: A
    real(r8), dimension(:), intent(INOUT) :: B

    ! Local Variables
    real(r8), dimension(SIZE(A)) :: Tmp

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Swap Arrays
    Tmp = A
    A = B
    B = Tmp

  END SUBROUTINE SWAP_REAL_R1

END MODULE SORT_MODULE
