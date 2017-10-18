#include "f90_assert.fpp"

module sort_utilities

  use,intrinsic :: iso_fortran_env, only: int8, int16, int32, int64, real32, real64
  use kinds, only: r8
  implicit none
  private

  public :: heap_sort, insertion_sort

  interface heap_sort
    procedure heap_sort_int8, heap_sort_int16, heap_sort_int32, heap_sort_int64
    procedure heap_sort_real32, heap_sort_real64
  end interface heap_sort

  interface insertion_sort
    module procedure insertion_sort_3r8r8, insertion_sort_ir8, insertion_sort_r8
  end interface insertion_sort


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

end module sort_utilities
