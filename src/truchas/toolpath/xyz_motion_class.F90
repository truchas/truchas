!!
!! XYZ_MOTION_CLASS
!!
!! This module defines the abstract base class XYZ_MOTION whose instances
!! describe a time dependent path in Cartesion space.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! June 2016
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  Objects of XYZ_MOTION class define a path R(T) on a time interval [T0, T1].
!!  They have the following type bound functions:
!!
!!    COORD(T) returns the position R(T) for T in [T0,T1]; the result is
!!        undefined for other T.
!!    START_TIME() returns T0
!!    FINAL_TIME() returns T1
!!    START_COORD() returns R(T0)
!!    FINAL_COORD() returns R(T1)
!!

module xyz_motion_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: xyz_motion
  contains
    procedure(coord), deferred :: coord
    procedure(start_coord), deferred :: start_coord
    procedure(start_coord), deferred :: final_coord
    procedure(start_time),  deferred :: start_time
    procedure(start_time),  deferred :: final_time
    procedure(partition),   deferred :: partition
  end type xyz_motion

  abstract interface
    pure function coord(this, t) result(r)
      import r8, xyz_motion
      class(xyz_motion), intent(in) :: this
      real(r8), intent(in) :: t
      real(r8) :: r(3)
    end function
    pure function start_coord(this) result(r)
      import r8, xyz_motion
      class(xyz_motion), intent(in) :: this
      real(r8) :: r(3)
    end function
    pure function start_time(this) result(t)
      import r8, xyz_motion
      class(xyz_motion), intent(in) :: this
      real(r8) :: t
    end function
    pure function partition(this, ds) result(t)
      import r8, xyz_motion
      class(xyz_motion), intent(in) :: this
      real(r8), intent(in) :: ds
      real(r8), allocatable :: t(:)
    end function
  end interface

end module xyz_motion_class
