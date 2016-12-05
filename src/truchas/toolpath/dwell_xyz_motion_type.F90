!!
!! DWELL_XYZ_MOTION_TYPE
!!
!! A concrete implementation of the abstract base class XYZ_MOTION.  This
!! implementation defines a pause, or dwell, in the motion for a time interval.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! June 2016
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  Objects of type DWELL_XYZ_MOTION are created using its constructor which
!!  has several possible invocations:
!!
!!    DWELL_XYZ_MOTION(R,T0,DT) - At coordinate R for time interval (T0,T0+DT)
!!    DWELL_XYZ_MOTION(R,T1) - At coordinate R for time interval (-\infty,T1)
!!    DWELL_XYZ_MOTION(R,T0) - At coordinate R for time interval (T0,\infty)
!!
!!  Keywords must be used to specify the T0, T1, DT arguments.  The last two
!!  invocations create special moves that are used to start and end a complete
!!  path.
!!

#include "f90_assert.fpp"

module dwell_xyz_motion_type

  use xyz_motion_class
  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, extends(xyz_motion), public :: dwell_xyz_motion
    private
    real(r8) :: r(3)
    real(r8) :: t0 = -huge(1.0_r8)
    real(r8) :: t1 = huge(1.0_r8)
  contains
    procedure :: coord
    procedure :: start_coord
    procedure :: final_coord => start_coord
    procedure :: start_time
    procedure :: final_time
    procedure :: partition
  end type dwell_xyz_motion

  !! Defined constructor
  interface dwell_xyz_motion
    procedure init
  end interface

  !! Barrier used to force use of argument keywords
  type :: kwarg_barrier
  end type

contains

  function init(r, kwarg, t0, t1, dt) result(this)
    real(r8), intent(in) :: r(:)
    type(kwarg_barrier), optional :: kwarg
    real(r8), intent(in), optional :: t0, t1, dt
    type(dwell_xyz_motion) :: this
    this%r = r
    if (present(t0)) then
      this%t0 = t0
      if (present(dt)) this%t1 = t0 + dt
    else if (present(t1)) then
      this%t1 = t1
    else
      INSIST(.false.)
    end if
  end function init

  pure function coord(this, t) result(r)
    class(dwell_xyz_motion), intent(in) :: this
    real(r8), intent(in) :: t
    real(r8) :: r(3)
    !ASSERT(t >= this%t0 .and. t <= this%t1)
    r = this%r
  end function coord

  pure function start_coord(this) result(r)
    class(dwell_xyz_motion), intent(in) :: this
    real(r8) :: r(3)
    r = this%r
  end function start_coord

  pure function start_time(this) result(t)
    class(dwell_xyz_motion), intent(in) :: this
    real(r8) :: t
    t = this%t0
  end function start_time

  pure function final_time(this) result(t)
    class(dwell_xyz_motion), intent(in) :: this
    real(r8) :: t
    t = this%t1
  end function final_time

  pure function partition(this, ds) result(t)
    class(dwell_xyz_motion), intent(in) :: this
    real(r8), intent(in) :: ds
    real(r8), allocatable :: t(:)
    if (this%t0 == -huge(1.0_r8)) then
      t = [this%t1]
    else if (this%t1 == huge(1.0_r8)) then
      t = [this%t0]
    else
      t = [this%t0, this%t1]
    end if
  end function partition

end module dwell_xyz_motion_type
