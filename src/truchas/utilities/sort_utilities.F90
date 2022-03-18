#include "f90_assert.fpp"

module sort_utilities

  use,intrinsic :: iso_fortran_env, only: int8, int16, int32, int64, real32, real64
  use kinds, only: r8
  implicit none
  private

  public :: heap_sort, insertion_sort, quick_sort

  interface heap_sort
    procedure heap_sort_int8, heap_sort_int16, heap_sort_int32, heap_sort_int64
    procedure heap_sort_real32, heap_sort_real64
  end interface heap_sort

  interface insertion_sort
    module procedure insertion_sort_3r8r8, insertion_sort_ir8, insertion_sort_r8
    module procedure insertion_sort_i4
  end interface insertion_sort

  interface quick_sort
    module procedure quick_sort_r8
  end interface quick_sort

contains

  subroutine heap_sort_int8 (array, perm)
    integer(int8), intent(in) :: array(:)
    integer, intent(out) :: perm(:)
#include "heap_sort_body.fpp"
  end subroutine heap_sort_int8

  subroutine heap_sort_int16 (array, perm)
    integer(int16), intent(in) :: array(:)
    integer, intent(out) :: perm(:)
#include "heap_sort_body.fpp"
  end subroutine heap_sort_int16

  subroutine heap_sort_int32 (array, perm)
    integer(int32), intent(in) :: array(:)
    integer, intent(out) :: perm(:)
#include "heap_sort_body.fpp"
  end subroutine heap_sort_int32

  subroutine heap_sort_int64 (array, perm)
    integer(int64), intent(in) :: array(:)
    integer, intent(out) :: perm(:)
#include "heap_sort_body.fpp"
  end subroutine heap_sort_int64

  subroutine heap_sort_real32 (array, perm)
    real(real32), intent(in) :: array(:)
    integer, intent(out) :: perm(:)
#include "heap_sort_body.fpp"
  end subroutine heap_sort_real32

  subroutine heap_sort_real64 (array, perm)
    real(real64), intent(in) :: array(:)
    integer, intent(out) :: perm(:)
#include "heap_sort_body.fpp"
  end subroutine heap_sort_real64

  subroutine insertion_sort_3r8r8 (x,key)
    real(r8), intent(inout) :: x(:,:),key(:)

    real(r8) :: tmp,tmpX(size(x,dim=1))
    integer  :: i,j

    do i = 2,size(key)
      tmp  = key(i)
      tmpX = x(:,i)
      j = i
      do while (j>1)
        ! fortran doesn't guarantee short-circuiting, so we segfault if the two checks are together
        if (key(j-1) <= tmp) exit
        key(j) = key(j-1)
        x(:,j) = x(:,j-1)
        j = j-1
      end do
      key(j) = tmp
      x(:,j) = tmpX
    end do

  end subroutine insertion_sort_3r8r8

  subroutine insertion_sort_ir8 (x,key)
    integer,  intent(inout) :: x(:)
    real(r8), intent(inout) :: key(:)

    integer  :: tmpX
    real(r8) :: tmp
    integer  :: i,j

    do i = 2,size(key)
      tmp  = key(i)
      tmpX = x(i)
      j = i
      do while (j>1)
        ! fortran doesn't guarantee short-circuiting, so we segfault if the two checks are together
        if (key(j-1)>tmp) then
          key(j) = key(j-1)
          x(j)   = x(j-1)
          j = j-1
        else
          exit
        end if
      end do
      key(j) = tmp
      x(j)   = tmpX
    end do

  end subroutine insertion_sort_ir8

  subroutine insertion_sort_r8 (key)
    real(r8), intent(inout) :: key(:)

    real(r8) :: tmp
    integer  :: i,j

    do i = 2,size(key)
      tmp = key(i)
      j = i
      do while (j>1)
        ! fortran doesn't guarantee short-circuiting, so we segfault if the two checks are together
        if (key(j-1) <= tmp) exit
        key(j) = key(j-1)
        j = j-1
      end do
      key(j) = tmp
    end do

  end subroutine insertion_sort_r8

  subroutine insertion_sort_i4 (key)
    integer(int32), intent(inout) :: key(:)

    integer(int32) :: tmp
    integer  :: i,j

    do i = 2,size(key)
      tmp = key(i)
      j = i
      do while (j>1)
        ! fortran doesn't guarantee short-circuiting, so we segfault if the two checks are together
        if (key(j-1) <= tmp) exit
        key(j) = key(j-1)
        j = j-1
      end do
      key(j) = tmp
    end do

  end subroutine


  subroutine quick_sort_r8(arr, perm, n)
    real(r8), intent(inout) :: arr(:)
    integer, intent(inout) :: perm(:)
    integer, intent(in) :: n

    call quicksort_(arr, perm, 1, n)

  contains

    ! quicksort - used for sorting wisp volumes
    recursive subroutine quicksort_(arr, perm, lo, hi)
      real(r8), intent(inout) :: arr(:)
      integer, intent(inout) :: perm(:)
      integer, intent(in) :: lo, hi

      integer :: left, right

      if (lo < hi) then
        call partition_(arr, perm, lo, hi, left, right)
        call quicksort_(arr, perm, lo, left-1)
        call quicksort_(arr, perm, right+1, hi)
      end if
    end subroutine quicksort_

    subroutine swap_(arr, perm, i, j)
      real(r8), intent(inout) :: arr(:)
      integer, intent(inout) :: perm(:)
      integer, intent(in) :: i, j

      real(r8) :: rtmp
      integer :: itmp

      itmp = perm(j)
      perm(j) = perm(i)
      perm(i) = itmp

      rtmp = arr(j)
      arr(j) = arr(i)
      arr(i) = rtmp
    end subroutine swap_

    ! partitioning algorithm designed to be robust for repeated elements
    ! (aka the dutch national flag problem)
    recursive subroutine partition_(arr, perm, lo, hi, left, right)
      real(r8), intent(inout) :: arr(:)
      integer, intent(inout) :: perm(:)
      integer, intent(in) :: lo, hi
      integer, intent(out) :: left, right

      real(r8) :: pivot
      integer :: i, j, k, mid

      ! median of three for pivot choice
      mid = (lo + hi) / 2
      if (arr(mid) < arr(lo)) call swap_(arr, perm, lo, mid)
      if (arr(hi) < arr(lo)) call swap_(arr, perm, lo, hi)
      if (arr(mid) < arr(hi)) call swap_(arr, perm, mid, hi)

      pivot = arr(hi)

      ! at the end of this loop
      ! 1) arr(lo:i-1) < pivot
      ! 2) arr(i:j-1) == pivot
      ! 3) arr(j:hi) > pivot
      i = lo
      j = lo
      k = hi
      do while (j <= k)
        if (arr(j) < pivot) then
          call swap_(arr, perm, i, j)
          i = i + 1
          j = j + 1
        else if (arr(j) > pivot) then
          call swap_(arr, perm, j, k)
          k = k - 1
        else
          j = j + 1
        end if
      end do
      left = i
      right = j-1

    end subroutine partition_

  end subroutine quick_sort_r8

end module sort_utilities
