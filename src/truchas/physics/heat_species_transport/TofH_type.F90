!!
!! TOFH_TYPE
!!
!! This module defines a type for representing the relation for temperature as
!! a function of enthalpy density, T(H), that derives from its given inverse
!! property mesh function H(T).
!!
!! Neil N. Carlson <nnc@lanl.gov.
!! Adapted for Fortran 2008, May 2015
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  The derived type TOFH has the following type bound procedures
!!
!!  INIT (HOFT, EPS [,MAX_TRY [,DELTA]]) initializes the object to represent
!!    the inverse of the given property mesh function HOFT, which is assumed
!!    to describe an increasing enthalpy density as a function of temperature
!!    relation.  EPS is the accuracy to which the temperature will be computed
!!    by COMPUTE.  The optional arguments control the recovery algorithm used
!!    when the root bracketing interval provided to COMPUTE is invalid; see
!!    COMPUTE for details.
!!
!!  COMPUTE (CELL, H, TMIN, TMAX, T) computes the temperature T as a function
!!    of enthalpy density H for the given cell index CELL.
!!    (Note that the mixture of materials can vary from cell to cell.)
!!    The computation involves finding the root of the function H - H(T).
!!    The given interval [TMIN, TMAX] must bracket the root T.  If MAX_TRY
!!    and DELTA were given to TOFH_INIT and the interval fails to bracket
!!    the temperature, then the interval will be expanded as many as MAX_TRY
!!    times seeking an interval that does.  The errant endpoint is shifted
!!    successively starting with the increment DELTA, and increasing by a
!!    factor of 10 each additional try.
!!
!! N.B. The given enthalpy relation H(T) is assumed to be a function of
!!      temperature alone; this will need to be extended to handle systems
!!      involving species concentrations.
!!

#include "f90_assert.fpp"

module TofH_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use prop_mesh_func_type
  use ridders_class
  implicit none
  private

  type, extends(ridders), public :: TofH
    private
    type(prop_mesh_func), pointer :: HofT => null()
    real(r8) :: H
    integer  :: cell
    !! Parameters for the algorithm that wraps Ridders root finding
    integer  :: max_try
    real(r8) :: delta
    !! Performance counters
    integer :: num_call = 0 ! number of successful calls
    integer :: max_itr = 0  ! max number of single-call Ridders iterations
    integer :: num_itr = 0  ! total number of Ridders iterations
    integer :: num_rec = 0  ! number of calls requiring bracketing recovery
    integer :: max_adj = 0  ! max number of single-call interval adjustments
    integer :: num_adj = 0  ! total number of interval adjustments
  contains
    procedure :: f          ! deferred procedure from ridders
    procedure :: init
    procedure :: compute
    procedure :: get_metrics
  end type TofH

contains

  function f (this, x) result (fx)
    class(TofH), intent(in) :: this
    real(r8), intent(in) :: x
    real(r8) :: fx
    call this%HofT%compute_value(this%cell, [x], fx)
    fx = this%H - fx
  end function

  subroutine init (this, HofT, eps, max_try, delta)
    class(TofH), intent(out) :: this
    type(prop_mesh_func), target :: HofT
    real(r8), intent(in) :: eps
    integer, intent(in), optional :: max_try
    real(r8), intent(in), optional :: delta
    this%HofT => HofT
    !call set_context_HofT (HofT)
    this%eps = eps
    this%maxitr = 100
    if (present(max_try)) then
      this%max_try = max_try
    else
      this%max_try = 1
    end if
    if (present(delta)) then
      this%delta = delta
    else
      this%delta = eps
    end if
  end subroutine init

  subroutine get_metrics (this, avg_itr, max_itr, rec_rate, avg_adj, max_adj)
    use parallel_communication, only: global_sum, global_maxval
    class(TofH), intent(in) :: this
    integer, intent(out), optional :: max_itr, max_adj
    real, intent(out), optional :: avg_itr, rec_rate, avg_adj
    if (present(avg_itr))  avg_itr  = real(global_sum(this%num_itr)) / max(1,global_sum(this%num_call))
    if (present(max_itr))  max_itr  = global_maxval(this%max_itr)
    if (present(rec_rate)) rec_rate = real(global_sum(this%num_rec)) / max(1,global_sum(this%num_call))
    if (present(avg_adj))  avg_adj  = real(global_sum(this%num_adj)) / max(1,global_sum(this%num_rec))
    if (present(max_adj))  max_adj  = global_maxval(this%max_adj)
  end subroutine get_metrics

  subroutine compute (this, cell, H, Tmin, Tmax, T)
    use truchas_logging_services, only: TLS_fatal
    class(TofH), intent(inout) :: this
    integer,  intent(in)  :: cell
    real(r8), intent(in)  :: H, Tmin, Tmax
    real(r8), intent(out) :: T
    integer :: n, stat
    character(100) :: errmsg
    real(r8) :: a, b, d
    this%cell = cell
    this%H = H
    !call set_context_cell_H (cell, H) ! parameters used by f
    call this%find_root (Tmin, Tmax, T, stat)
    if (stat < 0) then ! root not bracketed -- attempt to recover
      a = Tmin; b = Tmax; d = this%delta
      if (this%f(a) < 0.0_r8) then ! shift a to the left
        do n = 1, this%max_try
          a = a - d
          if (this%f(a) > 0.0_r8) exit
          d = 10*d
        end do
      else  ! then f(b) > 0; shift b to the right
        do n = 1, this%max_try
          b = b + d
          if (this%f(b) < 0.0_r8) exit
          d = 10*d
        end do
      end if
      if (n <= this%max_try) then ! root bracketed -- try again
        call this%find_root (a, b, T, stat)
        if (stat == 0) then
          this%num_rec = this%num_rec + 1
          this%num_adj = this%num_adj + n
          this%max_adj = max(this%max_adj, n)
        end if
      end if
    end if
    if (stat == 0) then
      this%num_call = this%num_call + 1
      this%num_itr = this%num_itr + this%numitr
      this%max_itr = max(this%max_itr, this%numitr)
    else if (stat < 0) then
      write(errmsg,'(2(a,es21.14),a)') 'root not bracketed: [', a, ',', b, ']'
      call TLS_fatal ('TofH_compute: ' // trim(errmsg))
    else
      write(errmsg,'(a,es10.4,2(a,es21.14))') &
        'convergence failure: error=', this%error, ', T=', T, ', H-H(T)=', this%f(T)
      call TLS_fatal ('TofH_compute: ' // trim(errmsg))
    end if
  end subroutine compute

end module TofH_type

