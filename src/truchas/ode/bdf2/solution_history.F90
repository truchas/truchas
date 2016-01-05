!!
!!  SOLUTION_HISTORY
!!
!!  Neil N. Carlson <nnc@newmexico.com> 7 Jul 2003
!!  Last revised 9 Aug 2006
!!
!!  This module provides a structure for maintaining the recent history of
!!  a solution procedure that is characterized by a (time) sequence of solution
!!  vectors, and procedures for performing polynomial interpolation based
!!  on that history.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  The module provides the derived data type HISTORY (with private components)
!!  and the following procedures:
!!
!!    CALL CREATE_HISTORY (THIS, MVEC, {T, X [, XDOT] | N})
!!
!!      TYPE(HISTORY), INTENT(OUT) :: THIS
!!      REAL(KIND=R8), INTENT(IN) :: T, X(:), XDOT(:)
!!      INTEGER, INTENT(IN) :: N, MVEC
!!
!!      Creates a new history structure THIS capable of maintaining MVEC
!!      solution vectors.  In the first variant, the vector X with time
!!      index T is recorded as the initial solution vector.  If the optional
!!      vector XDOT is also specified, it is recorded as the solution vector
!!      time derivative at the same time index T.  In the second variant, N
!!      specifies the length of the vectors to be maintained, but no solution
!!      vector is recorded.
!!
!!    CALL FLUSH_HISTORY (THIS, T, X [, XDOT])
!!
!!      TYPE(HISTORY), INTENT(INOUT) :: THIS
!!      REAL(KIND=R8), INTENT(IN)    :: T, X(:), XDOT(:)
!!
!!      Flushes the accumulated solution vectors from an existing history
!!      structure THIS, and records the solution vector X with time index T as
!!      the initial solution vector of a new history.  If XDOT is specified
!!      it is also recorded as the solution vector time derivative at time
!!      index T.  It is an error to pass an undefined history structure.
!!
!!    DEFINED(THIS) returns the value true if the history variable THIS is
!!      defined; otherwise the value false is returned.  Defined means that
!!      the components of the variable are consistently defined and allocated.
!!      This function is primarly intended to be used in assertion checks.
!!
!!    CALL DESTROY (THIS) deallocates all storage associated with the history
!!      structure THIS, and returns the variable to its default initialization
!!      state.
!!
!!    CALL RECORD_SOLUTION (THIS, T, X [, XDOT])
!!
!!      TYPE(HISTORY), INTENT(INOUT) :: THIS
!!      REAL(KIND=R8), INTENT(IN) :: T, X(:), XDOT(:)
!!
!!      Records the vector X with time index T as the most recent solution
!!      vector in the history structure THIS.  If the vector XDOT is present,
!!      it is recorded as the solution vector time derivative at the same time
!!      index.  The oldest solution vector (or the two oldest in the case XDOT
!!      is present) is discarded once the history is fully populated with MVEC
!!      vectors.  Note that when only one of a X/XDOT pair of vectors is
!!      discarded, it is the derivative vector that gets discarded.
!!
!!    CALL REVISE_HISTORY (THIS, INDEX[, X][, XDOT])
!!
!!      TYPE(HISTORY), INTENT(INOUT) :: THIS
!!      INTEGER, INTENT(IN) :: INDEX
!!      REAL(KIND=R8), INTENT(IN) :: X, XDOT
!!
!!      Revises the history of a selected vector component of the solution
!!      procedure.  INDEX is the component, X is the new most recent value
!!      of that component and XDOT is the new first divided difference.
!!      Either or both can be specified.  All unspecified divided differences
!!      for the component are set to zero.
!!
!!    CALL INTERPOLATE_SOLUTION (THIS, T, X, ORDER)
!!
!!      TYPE(HISTORY), INTENT(IN) :: THIS
!!      REAL(KIND=R8), INTENT(IN)  :: T
!!      REAL(KIND=R8), INTENT(OUT) :: X(:)
!!      INTEGER, INTENT(IN), OPTIONAL :: ORDER
!!
!!      Computes the interpolated (or extrapolated) vector X at time T using
!!      polynomial interpolation from the set of solution vectors maintained
!!      by the history THIS.  ORDER, if present, specifies the interpolation
!!      order using the ORDER + 1 most recent solution vectors; 1 for linear
!!      interpolation, 2 for quadratic, etc.  It is an error to request an
!!      order for which there is insufficient data.  If not specified, the
!!      maximal interpolation order is used given the available data; once
!!      the history is fully populated, the order of interpolation is MVEC-1.
!!
!!    MOST_RECENT_SOLUTION (THIS) RESULT (X)
!!
!!      TYPE(HISTORY), INTENT(IN) :: THIS
!!      REAL(KIND=R8), POINTER :: X(:)
!!
!!      Function returns a pointer to the most recent solution vector
!!      maintained by the history THIS.  NB: The target of this pointer
!!      should never be modified
!!
!!    MOST_RECENT_TIME (THIS) RESULT (T)
!!
!!      TYPE(HISTORY), INTENT(IN) :: THIS
!!      REAL(KIND=R8) :: T
!!
!!      Function returns the the time index T associated with the most
!!      recent solution vector maintained by the history THIS.
!!
!!    LAST_TIME_DELTA (THIS) RESULT (H)
!!
!!      TYPE(HISTORY), INTENT(IN) :: THIS
!!      REAL(KIND=R8) :: H
!!
!!      Function returns the difference between the most recent time index and
!!      the penultimate time index.
!!
!!    TIME_DELTAS (THIS) RESULT (H)
!!
!!      TYPE(HISTORY), INTENT(IN) :: THIS
!!      REAL(KIND=R8) :: H(NVEC)
!!
!!      Function returns an array H of time index differences.  The first
!!      element of H is the difference between the most recent time and the
!!      penultimate time.  The second element is the difference between the
!!      most recent time and the antepenultimate time, and so forth.  The
!!      length of the result equals one less than the number of solution
!!      vectors being maintained by the history THIS
!!
!!    HISTORY_SIZE (THIS) returns the number of solution vectors currently
!!      maintained in the history structure THIS.  The number will be
!!      between 0 and the value of MVEC used to create the structure.
!!
!!  IMPLEMENTATION NOTES:
!!  (1) The sequence of solution vectors is actually stored as a table of
!!      divided differences, which is easily updated, and which makes
!!      polynomial interpolation particularly simple to express.
!!  (2) Oldest/newest refers to the sequence in which vectors are recorded
!!      and not to the time values supplied.  Typically, the corresponding
!!      time sequence will be strictly increasing; anything else would be
!!      perverse.
!!  (3) History shifting (induced by recording a new vector) takes place
!!      through shifting pointers; vectors themselves are not copied.
!!      The natural implementation would then use an array of pointers (length
!!      MVEC), but this is not possible to do directly.  The indirect means
!!      is to box up the pointer in an auxillary type (PTR_BOX) and declare
!!      an array of this type.
!!

#include "f90_assert.fpp"

module solution_history

  implicit none
  private

  public :: history, defined, destroy
  public :: create_history, flush_history, revise_history
  public :: record_solution, interpolate_solution
  public :: most_recent_time, most_recent_solution, time_deltas, history_size
  public :: last_time_delta

  integer, parameter :: r8 = kind(1.0d0)

  type :: history
    private
    integer :: n = 0      ! vector length
    integer :: nvec = 0   ! number of vectors
    real(kind=r8), pointer :: t(:) => null()  ! vector times
    type(ptr_box), pointer :: d(:) => null()  ! divided differences
  end type history

  type :: ptr_box
    real(kind=r8), pointer :: ptr(:) => null()
  end type ptr_box

  interface create_history
    module procedure create_history_vec, create_history_size
  end interface

  interface destroy
    module procedure destroy_history
  end interface

  interface defined
    module procedure defined_history
  end interface

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! Specific procedures for CREATE_HISTORY
 !!
 !! These routines configure the history structure THIS to handle a maximum
 !! of MVEC solution vectors.  The first form records the vector X with time
 !! index T as the initial vector.  In the second form, only the length N of
 !! the vectors to be maintained is specified; no solution vector is recorded.
 !!

  subroutine create_history_vec (this, mvec, t, x, xdot)

    type(history), intent(out) :: this
    integer, intent(in) :: mvec
    real(kind=r8), intent(in) :: t, x(:)
    real(kind=r8), intent(in), optional :: xdot(:)

    call create_history_size (this, mvec, size(x))
    call record_solution (this, t, x, xdot)

  end subroutine create_history_vec

  subroutine create_history_size (this, mvec, n)

    type(history), intent(out) :: this
    integer, intent(in) :: mvec, n

    integer :: j

    ASSERT( mvec > 0 )
    ASSERT( n >= 0 )

    this%n = n
    this%nvec = 0
    allocate(this%t(mvec), this%d(mvec))
    do j = 1, mvec
      allocate(this%d(j)%ptr(this%n))
    end do

  end subroutine create_history_size

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! Specific procedure for generic DEFINED
 !!
 !! Return the value true if the history variable THIS is defined;  otherwise
 !! return the value false.  Defined means that the components of the variable
 !! are consistently defined and allocated.
 !!
 !! This function is primarly intended to be used in assertion checks.
 !!

  logical function defined_history (this)

    type(history), intent(in) :: this

    integer :: j

    defined_history = .false.
    if (.not. associated(this%t)) return
    if (.not. associated(this%d)) return
    if ( size(this%t) /= size(this%d) ) return
    if (this%nvec < 0 .or. this%nvec > size(this%t)) return
    if (this%n < 0) return
    do j = 1, size(this%d)
      if (.not.associated(this%d(j)%ptr)) return
      if (size(this%d(j)%ptr) /= this%n)  return
    end do
    defined_history = .true.

  end function defined_history

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! Specific procedure for generic DESTROY
 !!
 !! This routine deallocates any allocated components of the history structure
 !! THIS, and assigns default values to the remaining components, returning the
 !! variable to its default initialization state.
 !!

  subroutine destroy_history (this)

    type(history), intent(inout) :: this

    integer :: j
    type(history) :: default

    if (associated(this%t)) deallocate(this%t)
    if (associated(this%d)) then
      do j = 1, size(this%d)
        if (associated(this%d(j)%ptr)) deallocate(this%d(j)%ptr)
      end do
      deallocate(this%d)
    end if

    this = default  ! assign default initialization values

  end subroutine destroy_history

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! FLUSH_HISTORY
 !!
 !! This routine flushes the accumulated solution vectors from the history
 !! structure THIS, and records the solution vector X with time index T as
 !! the initial solution of a new history.  If the optional vector XDOT is
 !! specified, it is recorded as the solution vector time derivative at the
 !! same time index.
 !!

  subroutine flush_history (this, t, x, xdot)

    type(history), intent(inout) :: this
    real(kind=r8), intent(in) :: t, x(:)
    real(kind=r8), intent(in), optional :: xdot(:)

    ASSERT( defined(this) )

    this%nvec = 0
    call record_solution (this, t, x, xdot)

  end subroutine flush_history

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! RECORD_SOLUTION
 !!
 !! This subroutine records the vector X with time index T as the most recent
 !! solution vector in the history structure THIS.  If the vector XDOT is
 !! present, it is recorded as the solution vector time derivative at the same
 !! time index.  The oldest solution vector (or two oldest in the case XDOT is
 !! present) is discarded once the history is fully populated.  Note that it
 !! is actually the table of vector divided differences that is stored, rather
 !! than the vectors themselves.
 !!
 !! The treatment of the derivative data XDOT can be viewed as recording a
 !! second solution vector at a slightly greater time index than T such that
 !! the first divided difference is XDOT, and then taking the limit as this
 !! time index approaches T.  This is implemented here by first recording X
 !! as a new solution vector (and updating all the divided differences) and
 !! then introducing X again as a new solution vector at the same time index
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

  subroutine record_solution (this, t, x, xdot)

    type(history), intent(inout) :: this
    real(kind=r8), intent(in) :: t, x(:)
    real(kind=r8), intent(in), optional :: xdot(:)

    integer :: j
    type(ptr_box) :: tmp

    ASSERT( defined(this) )
    ASSERT( size(x) == this%n )

    this%nvec = min(1+this%nvec, size(this%d))  ! update the number of vectors

    !! Shift the divided differences
    tmp = this%d(this%nvec)  ! storage for oldest gets recycled for newest
    do j = this%nvec, 2, -1
      this%t(j) = this%t(j-1)
      this%d(j) = this%d(j-1)
    end do

    !! Insert the new vector
    this%t(1) = t
    this%d(1) = tmp
    this%d(1)%ptr = x

    !! Update the divided differences
    do j = 2, this%nvec
      this%d(j)%ptr = (this%d(j-1)%ptr - this%d(j)%ptr) / (this%t(1) - this%t(j))
    end do

    if (.not.present(xdot)) return

    ASSERT( size(xdot) == this%n )

    this%nvec = min(1+this%nvec, size(this%d))  ! update the number of vectors

    !! Shift the divided differences, except the first; the new vector and
    !! time index are the same as the most recent.
    tmp = this%d(this%nvec)  ! storage for oldest gets recycled for newest
    do j = this%nvec, 3, -1
      this%t(j) = this%t(j-1)
      this%d(j) = this%d(j-1)
    end do

    !! The first divided difference (same time index) is the specified derivative.
    this%t(2) = this%t(1)
    this%d(2) = tmp
    this%d(2)%ptr = xdot

    !! Update the rest of the divided differences.
    do j = 3, this%nvec
      this%d(j)%ptr = (this%d(j-1)%ptr - this%d(j)%ptr) / (this%t(1) - this%t(j))
    end do

  end subroutine record_solution

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! REVISE_HISTORY
 !!
 !! This subroutine revises the history of a selected vector component of the
 !! solution procedure.  INDEX is the component, X is the new most recent value
 !! of that component and XDOT, if present, is the new first divided difference.
 !! All higher order divided differences for the component are set to zero.
 !!

  subroutine revise_history (this, index, x, xdot)
  
    type(history), intent(inout) :: this
    integer, intent(in) :: index
    real(kind=r8), intent(in), optional :: x, xdot
    
    integer :: j
    
    ASSERT( defined(this) )
    ASSERT( index >= 1 .and. index <= this%n )
    ASSERT( this%nvec >= 1 )
    
    if (present(x)) then
      this%d(1)%ptr(index) = x
      do j = 2, this%nvec
        this%d(j)%ptr(index) = 0.0_r8
      end do
    end if
  
    if (present(xdot)) then
      ASSERT( this%nvec >= 2 )
      this%d(2)%ptr(index) = xdot
      do j = 3, this%nvec
        this%d(j)%ptr(index) = 0.0_r8
      end do
    end if
    
  end subroutine revise_history

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! INTERPOLATE_SOLUTION
 !!
 !! This subroutine computes the interpolated (or extrapolated) vector X at
 !! time T using polynomial interpolation from the set of solution vectors
 !! maintained by the history THIS.  ORDER, if present, specifies the
 !! interpolation order using the ORDER+1 most recent solution vectors; the
 !! default is the maximum order corresponding to the number of available
 !! vectors.
 !!

  subroutine interpolate_solution (this, t, x, order)

    type(history), intent(in)  :: this
    real(kind=r8),  intent(in)  :: t
    real(kind=r8),  intent(out) :: x(:)
    integer, intent(in), optional :: order

    integer :: j, k, interp_order
    real(kind=r8) :: value

    ASSERT( defined(this) )
    ASSERT( this%nvec > 0 )
    ASSERT( size(x) == this%n )

    !! Set the interpolation order.
    if (present(order)) then
      ASSERT( order >= 0 .and. order < this%nvec )
      interp_order = order
    else ! do maximal order of interpolation with available data.
      interp_order = this%nvec - 1
    end if

    do j = 1, size(x)
      value = this%d(interp_order+1)%ptr(j)
      do k = interp_order, 1, -1
        value = this%d(k)%ptr(j) + (t-this%t(k)) * value
      end do
      x(j) = value
    end do

  end subroutine interpolate_solution

  function most_recent_solution (this) result (x)
    type(history), intent(in) :: this
    real(kind=r8), pointer :: x(:)
    ASSERT( this%nvec > 0 )
    x => this%d(1)%ptr
  end function most_recent_solution

  function most_recent_time (this) result (t)
    type(history), intent(in) :: this
    real(kind=r8) :: t
    ASSERT( this%nvec > 0 )
    t = this%t(1)
  end function most_recent_time
  
  function last_time_delta (this) result (h)
    type(history), intent(in) :: this
    real(kind=r8) :: h
    ASSERT( this%nvec > 1 )
    h = this%t(1) - this%t(2)
  end function last_time_delta

  function time_deltas (this) result (h)
    type(history), intent(in) :: this
    real(kind=r8) :: h(this%nvec-1)
    ASSERT( this%nvec > 1 )
    h = this%t(1) - this%t(2:this%nvec)
  end function time_deltas

  integer function history_size (this)
    type(history), intent(in) :: this
    history_size = this%nvec
  end function history_size

end module solution_history
