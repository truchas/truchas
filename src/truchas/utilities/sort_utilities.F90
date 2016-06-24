#include "f90_assert.fpp"

module sort_utilities

  use,intrinsic :: iso_fortran_env, only: int8, int16, int32, int64, real32, real64
  implicit none
  private
  
  public :: heap_sort
  
  interface heap_sort
    procedure heap_sort_int8, heap_sort_int16, heap_sort_int32, heap_sort_int64
    procedure heap_sort_real32, heap_sort_real64
  end interface

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

end module sort_utilities
