!!
!! LINEAR_XYZ_MOTION_TYPE
!!
!! A concrete implementation of the abstract base class XYZ_MOTION.  This
!! implementation defines linear displacement from a point with speed and
!! possible initial acceleration and final deceleration.
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
!!  Objects of type LINEAR_XYZ_MOTION are created using its constructor which
!!  has several possible invocations:
!!
!!  LINEAR_XYZ_MOTION(T0, R0, DR, S) -- Linear path from position R0 at
!!    time T0 to position R+DR with linear speed S.
!!
!!  LINEAR_XYZ_MOTION(T0, R0, DR, S, A) -- The same linear path as above
!!    except that there is an initial constant linear acceleration A to
!!    speed S and a final constant linear deceleration A (same value) to rest.
!!
!!  LINEAR_XYZ_MOTION(T0, R0, DR, S, A, D) -- The same linear path as above
!!    but with a possibly different linear deceleration D.
!!

#include "f90_assert.fpp"

module linear_xyz_motion_type

  use xyz_motion_class
  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, extends(xyz_motion), public :: linear_xyz_motion
    private
    real(r8) :: t0, r0(3), dr(3)  ! independent data
    real(r8) :: ts, td, t1, ca, cd, fs, fd  ! derived data
  contains
    procedure :: coord
    procedure :: start_coord
    procedure :: final_coord
    procedure :: start_time
    procedure :: final_time
    procedure :: partition
  end type linear_xyz_motion

  !! Defined constructor
  interface linear_xyz_motion
    module procedure init
  end interface

contains

  function init(t0, r0, dr, s, a, d) result(this)

    real(r8), intent(in) :: t0, r0(:), dr(:), s
    real(r8), intent(in), optional :: a, d
    type(linear_xyz_motion) :: this

    real(r8) :: speed, accel, decel, dist, tmp, dt_accel, dt_decel

    ASSERT(size(r0) == 3)
    ASSERT(size(dr) == 3)
    ASSERT(s > 0)

    this%t0 = t0
    this%r0 = r0
    this%dr = dr

    dist = vector_length(dr)
    ASSERT(dist > 0)

    speed = s

    if (present(a)) then
      ASSERT(a > 0)
      accel = a
      decel = a
      if (present(d)) then
        ASSERT(d > 0)
        decel = d
      end if
    else
      ASSERT(.not.present(d))
      this%t1 = this%t0 + dist/speed
      this%ts = this%t0
      this%td = this%t1
      this%ca = 0
      this%cd = 0
      this%fs = 0
      this%fd = 1
      return
    end if

    tmp = (2*dist*accel*decel) / (accel + decel)

    if (tmp > speed**2) then ! usual case: move at speed for an interval
        dt_accel = speed/accel
        dt_decel = speed/decel
      this%t1 = this%t0 + (dist/speed + 0.5_r8*(dt_accel + dt_decel))
      this%ts = this%t0 + dt_accel
      this%td = this%t1 - dt_decel
      this%ca = accel/(2*dist)
      this%cd = decel/(2*dist)
      this%fs = this%ca * (this%ts-this%t0)**2
      this%fd = 1 - this%cd * (this%t1-this%td)**2
    else  ! exceptional case: do not attain speed
        tmp = sqrt(tmp)
        dt_accel = tmp/accel
        dt_decel = tmp/decel
      this%t1 = this%t0 + (dt_accel + dt_decel)
      this%td = this%t1 - dt_decel
      this%ts = this%td
      this%ca = accel/(2*dist)
      this%cd = decel/(2*dist)
      this%fs = this%ca * (this%ts-this%t0)**2
      this%fd = this%fs
    end if

  end function init

  pure function coord(this, t) result(r)
    class(linear_xyz_motion), intent(in) :: this
    real(r8), intent(in) :: t
    real(r8) :: r(3), f
    !ASSERT(t >= this%t0 .and. t <= this%t1)
    if (t <= this%ts) then
      f = this%ca*(t-this%t0)**2
    else if (t >= this%td) then
      f = 1 - this%cd*(this%t1-t)**2
    else
      f = this%fs*((this%td-t)/(this%td-this%ts)) + this%fd*((t-this%ts)/(this%td-this%ts))
    end if
    r = this%r0 + f*this%dr
  end function coord

  pure function start_coord(this) result(r)
    class(linear_xyz_motion), intent(in) :: this
    real(r8) :: r(3)
    r = this%r0
  end function start_coord

  pure function final_coord(this) result(r)
    class(linear_xyz_motion), intent(in) :: this
    real(r8) :: r(3)
    r = this%r0 + this%dr
  end function final_coord

  pure function start_time(this) result(t)
    class(linear_xyz_motion), intent(in) :: this
    real(r8) :: t
    t = this%t0
  end function start_time

  pure function final_time(this) result(t)
    class(linear_xyz_motion), intent(in) :: this
    real(r8) :: t
    t = this%t1
  end function final_time

  pure function partition(this, ds) result(t)

    class(linear_xyz_motion), intent(in) :: this
    real(r8), intent(in) :: ds
    real(r8), allocatable :: t(:)

    integer :: n, j, ns, nd
    real(r8) :: f

    n = ceiling(vector_length(this%dr)/ds)
    allocate(t(0:n))

    ns = floor(n*this%fs)
    t(0) = this%t0
    do j = 1, ns
      t(j) = this%t0 + sqrt(j/(n*this%ca))
    end do

    nd = floor(n*this%fd)
    do j = ns+1, nd
      f = j/real(n,r8)
      t(j) = this%ts*((this%fd-f)/(this%fd-this%fs)) + this%td*((f-this%fs)/(this%fd-this%fs))
    end do

    t(n) = this%t1
    do j = 1, n-nd-1
      t(n-j) = this%t1 - sqrt(j/(n*this%cd))
    end do

  end function partition

  pure function vector_length(x) result(len_x)
    real(r8), intent(in) :: x(:)
    real(r8) :: len_x, a, b, c, t
    a = abs(x(1)); b = abs(x(2)); c = abs(x(3))
    if (b > a) then
      if (c > b) then
        t = a; a = c; c = t ! swap A and C
      else
        t = a; a = b; b = t ! swap A and B
      end if
    else if (c > a) then
      t = a; a = c; c = t;  ! swap A and C
    end if
    if (a == 0.0_r8) then
      len_x = 0.0_r8
    else
      len_x = a * sqrt(1.0_r8 + ((b/a)**2 + (c/a)**2))
    end if
  end function vector_length

end module linear_xyz_motion_type
