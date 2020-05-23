!!
!! INVERSE_FUNC_CLASS
!!
!! This module defines an abstract base class for inverting a strictly
!! monotonic function.  Concrete extensions implement the function in the
!! deferred type bound function.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! April 2017
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  The abstract INVERSE_FUNC type provides a framework for inverting a function
!!  g(x); that is, given y, solving y = g(x) for x.  It does this by finding the
!!  root of the function f(x) = g(x) - y, and for that it extends the RIDDERS
!!  abstract type.  To apply the framework, define an extension of INVERSE_FUNC
!!  that defines the deferred type bound function G to implement g(x).
!!
!!  The derived type TOFH has the following type bound procedures
!!
!!  INIT(EPS[,MAXITR][,MAXADJ]) initializes the object. EPS is error tolerance
!!    used by COMPUTE in solving for x. The optional MAXITR argument specifies
!!    the maximum number of Ridders iterations in solving for x. The optional
!!    MAXADJ argument specifies the maximum number of attempts to discover a
!!    bracketing interval when the interval provided to COMPUTE is invalid.
!!    The default for MAXITR is 100 and the default for MAXADJ is 0.  See
!!    COMPUTE for details.
!!
!!  COMPUTE(Y, XMIN, XMAX, X [,STAT [,ERRMSG]]) solves g(x) = y for x. The
!!    given interval [XMIN, XMAX] should bracket the solution X. If MAXADJ > 0
!!    and the interval fails to bracket X, then the interval will be expanded
!!    as many as MAXADJ times seeking an interval that does.  The interval is
!!    shifted and its length doubled on each attempt.  STAT returns 0 if X
!!    was computed to the desired accuracy.  It returns 1 if the solution
!!    iteration failed to converge within a set maximum number of iterations.
!!    It returns -1 if it was unable to discover a bracketing interval within
!!    the allowed number of attempts.
!!
!!  The object maintains several performance counters as public data components:
!!
!!    %NUM_CALL -- number of COMPUTE calls
!!    %MAX_ITR  -- maximum number of solution iterations taken (single call)
!!    %NUM_ITR  -- total number of solution iterations (all calls)
!!    %NUM_REC  -- number of COMPUTE calls undertaking bracketing recovery
!!    %MAX_ADJ  -- maximum number of interval adjustments (single call)
!!    %NUM_ADJ  -- total number of interval adjustments (all calls)
!!
!!  Associated with the last COMPUTE call are the public data components:
!!
!!    %ERROR    -- computed bound on the solution error
!!    %NUMITR   -- number of iterations
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module inverse_func_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use ridders_class
  implicit none
  private

  type, abstract, extends(ridders), public :: inverse_func
    real(r8), private :: y = 0.0_r8
    !! Parameters for the algorithm that wraps Ridders root finding
    integer  :: maxadj = 0
    !! Performance counters
    integer :: num_call = 0 ! number of successful calls
    integer :: max_itr = 0  ! max number of single-call Ridders iterations
    integer :: num_itr = 0  ! total number of Ridders iterations
    integer :: num_rec = 0  ! number of calls requiring bracketing recovery
    integer :: max_adj = 0  ! max number of single-call interval adjustments
    integer :: num_adj = 0  ! total number of interval adjustments
  contains
    procedure :: init
    procedure :: compute
#ifdef GNU_PR88043
    procedure :: f ! deferred procedure from ridders, f(x) = g(x) - y
#else
    procedure, non_overridable :: f ! deferred procedure from ridders, f(x) = g(x) - y
#endif
    procedure(func),  deferred :: g
  end type inverse_func

  abstract interface
    function func (this, x) result (gx)
      import inverse_func, r8
      class(inverse_func), intent(in) :: this
      real(r8), intent(in) :: x
      real(r8) :: gx
    end function
  end interface

contains

  function f(this, x) result(fx)
    class(inverse_func), intent(in) :: this
    real(r8), intent(in) :: x
    real(r8) :: fx
    fx = this%g(x) - this%y
  end function

  subroutine init(this, eps, maxitr, maxadj)
    class(inverse_func), intent(out) :: this
    real(r8), intent(in) :: eps
    integer, intent(in), optional :: maxitr, maxadj
    this%eps = eps
    this%maxitr = 100
    if (present(maxitr)) this%maxitr = maxitr
    if (present(maxadj)) this%maxadj = maxadj
  end subroutine init

  subroutine compute(this, y, xmin, xmax, x, stat, errmsg)
#ifdef NAGFOR
    use,intrinsic :: f90_unix, only: exit
#endif
    use,intrinsic :: iso_fortran_env, only: error_unit
    class(inverse_func), intent(inout) :: this
    real(r8), intent(in)  :: y, xmin, xmax
    real(r8), intent(out) :: x
    integer, intent(out), optional :: stat
    character(:), intent(out), allocatable, optional :: errmsg
    integer :: ierr
    character(100) :: string
    real(r8) :: a, b
    this%y = y
    a = xmin; b = xmax
    call this%find_root(a, b, x, ierr)
    if (ierr < 0 .and. this%maxadj > 0) then ! root not bracketed; attempt to recover
      call this%bracket_root(a, b, this%maxadj, ierr)
      if (ierr == 0) then ! root bracketed -- try again
        call this%find_root (a, b, x, ierr)
        if (ierr == 0) then
          this%num_rec = this%num_rec + 1
          this%num_adj = this%num_adj + this%numadj
          this%max_adj = max(this%max_adj, this%numadj)
        end if
      end if
    end if
    if (present(stat)) stat = ierr
    if (ierr == 0) then
      this%num_call = this%num_call + 1
      this%num_itr = this%num_itr + this%numitr
      this%max_itr = max(this%max_itr, this%numitr)
    else
      if (ierr < 0) then
        write(string,'(2(a,es21.14),a)') 'root not bracketed: [', a, ',', b, ']'
      else
        write(string,'(a,es10.4,2(a,es21.14))') &
          'convergence failure: error=', this%error, ', x=', x, ', g(x)-y=', this%f(x)
      end if
      if (present(stat)) then
        if (present(errmsg)) errmsg = trim(string)
      else
        write(error_unit,'(2a)') 'inverse_func%compute: ', trim(string)
        call exit(1)
      end if
    end if

  end subroutine compute

end module inverse_func_class
