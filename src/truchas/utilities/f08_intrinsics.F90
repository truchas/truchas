!!
!! F08_INTRINSICS
!!
!! This module provides Fortran 2008 intrinsic functions the Intel
!! compiler does not yet support. Should be deleted as these functions
!! find support.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! March 2017
!!

module f08_intrinsics

  implicit none
  private

  public :: findloc

  interface findloc
    module procedure findloc_logical, findloc_integer
  end interface findloc

contains

  ! Implements a subset of the features present in the new intrinsic findloc
  pure function findloc_logical (array, value, back) result(i)

    logical, intent(in) :: array(:), value
    logical, intent(in), optional :: back
    integer :: i

    logical :: backh

    backh = .false.
    if (present(back)) backh = back

    if (backh) then
      do i = size(array),1,-1
        if (array(i) .eqv. value) return
      end do
    else
      do i = 1,size(array)
        if (array(i) .eqv. value) return
      end do
    end if
    i = 0

  end function findloc_logical

  pure function findloc_integer (array, value, back) result(i)

    integer, intent(in) :: array(:), value
    logical, intent(in), optional :: back
    integer :: i

    logical :: backh

    backh = .false.
    if (present(back)) backh = back

    if (backh) then
      do i = size(array),1,-1
        if (array(i) == value) return
      end do
    else
      do i = 1,size(array)
        if (array(i) == value) return
      end do
    end if
    i = 0

  end function findloc_integer

end module f08_intrinsics
