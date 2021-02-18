!!
!! STATE_HISTORY_TYPE
!!
!! This module provides a structure for maintaining the recent history of an
!! a solution procedure that is characterized by a (time) sequence of state
!! vectors, and methods for performing polynomial interpolation based on that
!! history.
!!
!! Neil N. Carlson <neil.n.carlson@gmail.com> 7 Jul 2003
!! Revision for Fortran 2003, 2 Apr 2011, 27 Dec 2011
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 2003, 2011, 2013  Neil N. Carlson
!!
!! Permission is hereby granted, free of charge, to any person obtaining a
!! copy of this software and associated documentation files (the "Software"),
!! to deal in the Software without restriction, including without limitation
!! the rights to use, copy, modify, merge, publish, distribute, sublicense,
!! and/or sell copies of the Software, and to permit persons to whom the
!! Software is furnished to do so, subject to the following conditions:
!!
!! The above copyright notice and this permission notice shall be included
!! in all copies or substantial portions of the Software.
!!
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
!! THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
!! DEALINGS IN THE SOFTWARE.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! This module defines the derived data type STATE_HISTORY (with private
!! components) with the following methods as type-bound procedures.
!!
!!  INIT (MVEC, {T, X [, XDOT] | VLEN}) initializes the history structure
!!    to maintain up to MVEC state vectors.  In the first variant, the
!!    vector X with time index T is recorded as the initial vector of a new
!!    history.  If the optional vector XDOT is also specified, it is recorded
!!    as the state vector time derivative at the same time index T.  In the
!!    second variant, VLEN specifies the length of the vectors to be maintained
!!    but no state vector is recorded.  It is permissible to re-initialize a
!!    history structure that had been previously initialized.
!!
!!  FLUSH (T, X [, XDOT]) flushes the accumulated state vectors and records the
!!    state vector X with time index T as the initial state vector of a new
!!    history.  If XDOT is specified, it is also recorded as the state vector
!!    time derivative at time index T.  This differs from re-initializing in
!!    that the maximum number of vectors and their length are not changed.
!!
!!  REAL_KIND() returns the kind parameter value expected of all real arguments.
!!
!!  DEPTH() returns the number of state vectors currently stored. This number
!!    will vary between 0 and the value of MVEC used to initialize the object.
!!
!!  MAX_DEPTH() returns the maximum number state vectors that can be stored.
!!    This number is the value of MVEC used to initialize the object.
!!
!!  STATE_SIZE() returns the length of the state vectors stored by the object.
!!    This number will equal the size of the state vector or the value of VLEN
!!    that was used to initialize the object.
!!
!!  LAST_TIME() returns the time index of the last recorded state vector.
!!
!!  GET_LAST_STATE_COPY (COPY) copies the last recorded state vector into
!!    the array COPY, whose length should equal STATE_SIZE().
!!
!!  GET_LAST_STATE_VIEW (VIEW) associates the pointer VIEW with the last
!!    recorded state vector.  THIS SHOULD BE USED WITH GREAT CAUTION.  The
!!    target of this pointer should never be modified or deallocated.  The
!!    pointer will cease to reference the last recorded state vector when
!!    the history is subsequently modified through calls to RECORD_STATE,
!!    FLUSH or INIT.
!!
!!  TIME_DELTAS() returns an array containing the time index differences: the
!!    first element is the difference between the last and penultimate times;
!!    the second element is the difference between the last and antepenultimate
!!    times, and so forth.  The length of the result equals DEPTH()-1.  It is
!!    an error to call this method if DEPTH() is less than 2.
!!
!!  RECORD_STATE (T, X [, XDOT]) records the vector X with time index T as
!!    the last state vector in the history.  If the vector XDOT is present,
!!    it is recorded as the state vector time derivative with the same time
!!    index.  The oldest state vector (or the two oldest in the case XDOT
!!    is present) is discarded once the history is fully populated with MVEC
!!    vectors.  Note that when only one of a X/XDOT pair of vectors is
!!    discarded, it is effectively the derivative vector that gets discarded.
!!
!!  INTERP_STATE (T, X [, FIRST] [, ORDER]) computes the interpolated state
!!    vector at time T from the set of vectors maintained by the history, and
!!    returns the result in the user-supplied array X. Polynomial interpolation
!!    is used, and ORDER, if present, specifies the order using the ORDER+1
!!    most recent vectors; 1 for linear interpolation, 2 for quadratic, etc.
!!    It is an error to request an order for which there is insufficient data.
!!    If not specified, the maximal interpolation order is used given the
!!    available data; once the history is fully populated, the interpolation
!!    order is MVEC-1.  Typically the array X would have the same size as the
!!    stored state vectors, but more generally X may return any contiguous
!!    segment of the interpolated state starting at index FIRST (default 1)
!!    and length the size of X.
!!
!!  DEFINED() returns true if the object is well-defined; otherwise it returns
!!    false.  Defined means that the data components of the object are properly
!!    and consistently defined.  DEFINED() should return true at any time after
!!    the INIT method has been called.  This method is primarily intended to be
!!    used in debugging situations and is used internally in assertion checks.
!!
!! REMARKS
!!
!! Oldest/newest refers to the sequence in which vectors are recorded and not
!! to the time values supplied.  Typically, the corresponding time sequence
!! will be strictly increasing; anything else would be perverse.
!!
!! When recording the time derivative (XDOT) of a solution vector, this vector
!! counts as an additional vector when considered by DEPTH() and MAX_DEPTH().
!!
!! IMPLEMENTATION NOTES
!!
!! The sequence of solution vectors is actually stored as a table of
!! divided differences, which is easily updated and which makes polynomial
!! interpolation particularly simple to express.
!!
!! History shifting (induced by recording a new vector) takes place through
!! shifting pointers; vectors themselves are not copied.  The natural
!! implementation would use an array of pointers to vectors, but this is not
!! possible to do directly.  The indirect means is to box up the pointer in
!! an auxillary type (PTR_BOX) and declare an array of this type.
!!

#include "f90_assert.fpp"

module new_state_history_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use vector_class
  implicit none
  private

  type :: vec_box
    class(vector), allocatable :: vec
  end type vec_box

  type, public :: state_history
    private
    integer :: nvec = 0   ! number of vectors
    integer :: mvec = 0   ! maximum number of vectors
    real(r8), allocatable :: t(:) ! vector times
    type(vec_box), allocatable :: d(:)  ! divided differences
  contains
    procedure, private :: init_vector, init_state, init_copy
    generic   :: init => init_vector, init_state, init_copy
    procedure :: defined
    procedure :: flush
    procedure :: real_kind
    procedure :: depth
    procedure :: max_depth
    procedure :: last_time
    procedure :: get_last_state_view
    procedure :: get_last_state_copy
    procedure :: time_deltas
    procedure :: record_state
    procedure :: interp_state
    generic :: assignment(=) => copy
    procedure, private :: copy
  end type state_history

contains

  subroutine init_vector(this, mvec, vec)
    class(state_history), intent(out) :: this
    integer, intent(in) :: mvec
    class(vector), intent(in) :: vec
    integer :: j
    ASSERT(mvec > 0)
    this%mvec = mvec
    this%nvec = 0
    allocate(this%t(mvec), this%d(mvec))
    do j = 1, mvec
      call vec%clone(this%d(j)%vec)
    end do
  end subroutine init_vector

  subroutine init_state(this, mvec, t, x, xdot)
    class(state_history), intent(out) :: this
    integer, intent(in) :: mvec
    real(r8), intent(in) :: t
    class(vector), intent(in) :: x
    class(vector), intent(in), optional :: xdot
    call init_vector(this, mvec, x)
    call record_state(this, t, x, xdot)
  end subroutine init_state

  subroutine init_copy(this, src)
    class(state_history), intent(out) :: this
    class(state_history), intent(in)  :: src
    integer :: j
    this%mvec = src%mvec
    this%nvec = src%nvec
    allocate(this%t(this%mvec), this%d(this%mvec))
    do j = 1, this%mvec
      call src%d(1)%vec%clone(this%d(j)%vec)
    end do
    do j = 1, this%nvec
      this%t(j) = src%t(j)
      call this%d(j)%vec%copy(src%d(j)%vec)
    end do
  end subroutine init_copy

  subroutine flush (this, t, x, xdot)
    class(state_history), intent(inout) :: this
    real(r8), intent(in) :: t
    class(vector), intent(in) :: x
    class(vector), intent(in), optional :: xdot
    this%nvec = 0
    call record_state(this, t, x, xdot)
  end subroutine flush

  pure integer function real_kind (this)
    class(state_history), intent(in) :: this
    real_kind = kind(this%t)
  end function real_kind

  integer function depth (this)
    class(state_history), intent(in) :: this
    depth = this%nvec
  end function depth

  integer function max_depth (this)
    class(state_history), intent(in) :: this
    max_depth = 0
    if (allocated(this%t)) max_depth = size(this%t)
  end function max_depth

  subroutine get_last_state_view (this, view)
    class(state_history), intent(in), target :: this
    class(vector), pointer :: view
    ASSERT(this%nvec > 0)
    view => this%d(1)%vec
  end subroutine get_last_state_view

  subroutine get_last_state_copy (this, copy)
    class(state_history), intent(in) :: this
    class(vector), intent(inout) :: copy
    ASSERT(this%nvec > 0)
    !ASSERT(size(copy) == this%vlen)
    !copy = this%d(1)%vec
    call copy%copy(this%d(1)%vec)
  end subroutine get_last_state_copy

  function last_time (this) result (t)
    class(state_history), intent(in) :: this
    real(r8) :: t
    ASSERT(this%nvec > 0)
    t = this%t(1)
  end function last_time

  function time_deltas (this) result (deltas)
    class(state_history), intent(in) :: this
    real(r8) :: deltas(this%nvec-1)
    ASSERT(this%nvec > 1)
    deltas = this%t(1) - this%t(2:this%nvec)
  end function time_deltas

  subroutine interp_state (this, t, x, order)

    class(state_history), intent(in)  :: this
    real(r8), intent(in)  :: t
    class(vector), intent(inout) :: x
    integer, intent(in), optional :: order

    integer :: k, interp_order

    ASSERT(defined(this))
    ASSERT(this%nvec > 0)

    !! Set the interpolation order.
    if (present(order)) then
      ASSERT(order >= 0 .and. order < this%nvec)
      interp_order = order
    else ! do maximal order of interpolation with available data.
      interp_order = this%nvec - 1
    end if

!    do j = 1, size(x)
      !value = this%d(interp_order+1)%vec(j)
      call x%copy(this%d(interp_order+1)%vec)
      do k = interp_order, 1, -1
        !value = this%d(k)%vec(j) + (t-this%t(k)) * value
        call x%update(1.0_r8, this%d(k)%vec, (t-this%t(k)))
      end do
!      x(j) = value
!    end do

  end subroutine interp_state

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! RECORD_STATE
 !!
 !! This subroutine records the vector X with time index T as the most recent
 !! state vector in the history structure THIS.  If the vector XDOT is
 !! present, it is recorded as the solution vector time derivative at the same
 !! time index.  The oldest solution vector (or two oldest in the case XDOT is
 !! present) is discarded once the history is fully populated.  Note that it
 !! is actually the table of vector divided differences that is stored, rather
 !! than the vectors themselves.
 !!
 !! The treatment of the derivative data XDOT can be viewed as recording a
 !! second state vector at a slightly greater time index than T such that
 !! the first divided difference is XDOT, and then taking the limit as this
 !! time index approaches T.  This is implemented here by first recording X
 !! as a new state vector (and updating all the divided differences) and
 !! then introducing X again as a new state vector at the same time index
 !! but with the first divided difference set to XDOT, and then computing the
 !! the higher order divided differences as usual.
 !!
 !! Interpolation using the resulting history works as expected: evaluating the
 !! interpolant at the vector's time index yields the vector, and evaluating
 !! the derivative of the interpolant at that time index would yield XDOT.
 !!
 !! Note that when only one of a X/XDOT pair of vectors is discarded, it is the
 !! derivative vector (which exists in the higher order divided difference) that
 !! gets discarded.
 !!

  subroutine record_state (this, t, x, xdot)

    class(state_history), intent(inout) :: this
    real(r8), intent(in) :: t
    class(vector), intent(in) :: x
    class(vector), intent(in), optional :: xdot

    integer :: j
    type(vec_box) :: tmp

    ASSERT(defined(this))

    this%nvec = min(1+this%nvec, this%mvec)  ! update the number of vectors

    !! Shift the divided differences
    !tmp = this%d(this%nvec)  ! storage for oldest gets recycled for newest
    call move_alloc(this%d(this%nvec)%vec, tmp%vec)
    do j = this%nvec, 2, -1
      this%t(j) = this%t(j-1)
      !this%d(j) = this%d(j-1)
      call move_alloc(this%d(j-1)%vec, this%d(j)%vec)
    end do

    !! Insert the new vector
    this%t(1) = t
    !this%d(1) = tmp
    !this%d(1)%vec = x
    call move_alloc(tmp%vec, this%d(1)%vec)
    call this%d(1)%vec%copy(x)

    !! Update the divided differences
    do j = 2, this%nvec
      !this%d(j)%vec = (this%d(j-1)%vec - this%d(j)%vec) / (this%t(1) - this%t(j))
      !a = 1.0_r8 / (this%t(1) - this%t(j))
      !call this%d(j)%vec%update(a, this%d(j-1)%vec, -a)
      call this%d(j)%vec%update(-1.0_r8, this%d(j-1)%vec)
      call this%d(j)%vec%scale(1.0_r8/(this%t(j) - this%t(1)))
    end do

    if (.not.present(xdot)) return
    if (this%mvec == 1) return  ! no room to store derivative info

    this%nvec = min(1+this%nvec, this%mvec)  ! update the number of vectors

    !! Shift the divided differences, except the first; the new vector and
    !! time index are the same as the most recent.
    !tmp = this%d(this%nvec)  ! storage for oldest gets recycled for newest
    call move_alloc(this%d(this%nvec)%vec, tmp%vec) ! storage for oldest gets recycled for newest
    do j = this%nvec, 3, -1
      this%t(j) = this%t(j-1)
      !this%d(j) = this%d(j-1)
      call move_alloc(this%d(j-1)%vec, this%d(j)%vec)
    end do

    !! The first divided difference (same time index) is the specified derivative.
    this%t(2) = this%t(1)
    !this%d(2) = tmp
    !this%d(2)%vec = xdot
    call move_alloc(tmp%vec, this%d(2)%vec)
    call this%d(2)%vec%copy(xdot)

    !! Update the rest of the divided differences.
    do j = 3, this%nvec
      !this%d(j)%vec = (this%d(j-1)%vec - this%d(j)%vec) / (this%t(1) - this%t(j))
      call this%d(j)%vec%update(-1.0_r8, this%d(j-1)%vec)
      call this%d(j)%vec%scale(1.0_r8/(this%t(j) - this%t(1)))
    end do

  end subroutine record_state

  logical function defined (this)
    class(state_history), intent(in) :: this
    integer :: j
    defined = .false.
    if (.not.allocated(this%t)) return
    if (.not.allocated(this%d)) return
    if (size(this%t) /= this%mvec) return
    if (size(this%d) /= this%mvec) return
    if (this%nvec < 0 .or. this%nvec > this%mvec) return
    do j = 1, this%mvec
      if (.not.allocated(this%d(j)%vec)) return
    end do
    defined = .true.
  end function defined

  !! Defined assignment subroutine
  subroutine copy(lhs, rhs)
    class(state_history), intent(inout) :: lhs
    class(state_history), intent(in) :: rhs
    integer :: j
    if (lhs%mvec /= rhs%mvec) then
      call init_copy(lhs, rhs)
    else
      lhs%nvec = rhs%nvec
      do j = 1, lhs%nvec
        lhs%t(j) = rhs%t(j)
        !lhs%d(j)%vec = rhs%d(j)%vec
        call lhs%d(j)%vec%copy(rhs%d(j)%vec)
      end do
    end if
  end subroutine copy

end module new_state_history_type
