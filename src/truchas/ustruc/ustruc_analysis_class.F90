!!
!! USTRUC_ANALYSIS_CLASS
!!
!! This module provides the abstract base class USTRUC_ANALYSIS that defines
!! the interface to the low-level microstructure analysis component used by
!! the microstructure modeling kernel as implemented by the USTRUC_MODEL type.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! August 2014
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module ustruc_analysis_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use,intrinsic :: iso_fortran_env, only: int8
  implicit none
  private

  type, abstract, public :: ustruc_analysis
    integer :: n  ! number of points
  contains
    procedure(update), deferred :: set_state
    procedure(update), deferred :: update_state
    procedure(has),    deferred :: has
    procedure(comp_list), deferred :: get_comp_list
    generic :: get => getl1, geti1, getr1, getr2
    procedure(getl1), deferred :: getl1
    procedure(geti1), deferred :: geti1
    procedure(getr1), deferred :: getr1
    procedure(getr2), deferred :: getr2
    procedure(serialize), deferred :: serialize
    procedure(deserialize), deferred :: deserialize
    procedure, nopass :: vector_magnitude
  end type

  abstract interface
    subroutine update(this, t, temp, temp_grad, frac, invalid)
      import ustruc_analysis, r8
      class(ustruc_analysis), intent(inout) :: this
      real(r8), intent(in) :: t, temp(:), temp_grad(:,:), frac(:)
      logical,  intent(in) :: invalid(:)
    end subroutine
    logical function has(this, name)
      import ustruc_analysis
      class(ustruc_analysis), intent(in) :: this
      character(*), intent(in) :: name
    end function
    subroutine comp_list(this, list)
      import ustruc_analysis
      class(ustruc_analysis), intent(in) :: this
      integer, allocatable, intent(out) :: list(:)
    end subroutine
    subroutine getl1(this, name, array)
      import ustruc_analysis
      class(ustruc_analysis), intent(in) :: this
      character(*), intent(in) :: name
      logical, intent(out) :: array(:)
    end subroutine
    subroutine geti1(this, name, array, invalid)
      import ustruc_analysis, r8
      class(ustruc_analysis), intent(in) :: this
      character(*), intent(in) :: name
      integer, intent(out) :: array(:)
      logical, intent(out), optional :: invalid(:)
    end subroutine
    subroutine getr1(this, name, array, invalid)
      import ustruc_analysis, r8
      class(ustruc_analysis), intent(in) :: this
      character(*), intent(in) :: name
      real(r8), intent(out) :: array(:)
      logical, intent(out), optional :: invalid(:)
    end subroutine
    subroutine getr2(this, name, array, invalid)
      import ustruc_analysis, r8
      class(ustruc_analysis), intent(in) :: this
      character(*), intent(in) :: name
      real(r8), intent(out) :: array(:,:)
      logical, intent(out), optional :: invalid(:)
    end subroutine
    subroutine serialize(this, cid, array)
      import ustruc_analysis, int8
      class(ustruc_analysis), intent(in) :: this
      integer, intent(in) :: cid
      integer(int8), allocatable, intent(out) :: array(:,:)
    end subroutine
    subroutine deserialize(this, cid, array)
      import ustruc_analysis, int8
      class(ustruc_analysis), intent(inout) :: this
      integer, intent(in) :: cid
      integer(int8), intent(in) :: array(:,:)
    end subroutine
  end interface

contains

  !! Numerically robust procedure to compute the magnitude of a 2 or 3-vector.

  pure function vector_magnitude(v) result(vmag)

    real(r8), intent(in) :: v(:)  ! length 2 or 3
    real(r8) :: vmag

    real(r8) :: a, b, c, t

    select case (size(v))
    case (2)
      a = abs(v(1))
      b = abs(v(2))
      !! Swap largest value to A.
      if (b > a) then
        t = a
        a = b
        b = t
      end if
      vmag = 0.0_r8
      if (a > 0.0_r8) vmag = a * sqrt(1.0_r8 + (b/a)**2)
    case (3)
      a = abs(v(1))
      b = abs(v(2))
      c = abs(v(3))
      !! Swap largest value to A.
      if (b > a) then
        if (c > b) then
          t = a
          a = c
          c = t
        else
          t = a
          a = b
          b = t
        end if
      else if (c > a) then
        t = a
        a = c
        c = t
      end if
      vmag = 0.0_r8
      if (a > 0.0_r8) vmag = a * sqrt(1.0_r8 + ((b/a)**2 + (c/a)**2))
    end select

  end function vector_magnitude

end module ustruc_analysis_class
