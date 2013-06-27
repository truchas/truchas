!!
!! TOFH_TYPE
!!
!! This module defines a type for representing the relation for temperature as
!! a function of enthalpy density, T(H), that derives from its given inverse
!! property mesh function H(T).
!! 
!!  CALL TOFH_INIT (THIS, HOFT, EPS[, MAX_TRY[, DELTA]]) initializes THIS to
!!    represent the inverse of the given property mesh function HOFT, which is
!!    assumed to describe an increasing enthalpy density as a function of
!!    temperature relation.  EPS is the accuracy to which the temperature will
!!    be computed by TOFH_COMPUTE.  The optional arguments control the recovery
!!    algorithm used when the root bracketing interval provided to TOFH_COMPUTE
!!    is invalid; see TOFH_COMPUTE for details.
!!
!!  CALL TOFH_COMPUTE (THIS, CELL, H, TMIN, TMAX, T) computes the temperature
!!    T as a function of enthalpy density H for the given cell index CELL.
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

  use kinds, only: r8
  use property_mesh_function
  use ridders_type
  use TofH_callback
  implicit none
  private
  
  type, public :: TofH
    private
    type(prop_mf), pointer :: HofT => null()
    type(ridders) :: root_finder
    integer  :: max_try
    real(r8) :: delta
    !! Performance counters
    integer :: num_call = 0 ! number of successful calls
    integer :: max_itr = 0  ! max number of single-call Ridders iterations
    integer :: num_itr = 0  ! total number of Ridders iterations
    integer :: num_rec = 0  ! number of calls requiring bracketing recovery
    integer :: max_adj = 0  ! max number of single-call interval adjustments
    integer :: num_adj = 0  ! total number of interval adjustments
  end type TofH
  
  public :: TofH_init
  public :: TofH_compute
  public :: TofH_eps
  public :: TofH_get_metrics
  
contains

  subroutine TofH_init (this, HofT, eps, max_try, delta)
    type(TofH), intent(out) :: this
    type(prop_mf), target :: HofT
    real(r8), intent(in) :: eps
    integer, intent(in), optional :: max_try
    real(r8), intent(in), optional :: delta
    this%HofT => HofT
    call set_context_HofT (HofT)
    this%root_finder%eps = eps
    this%root_finder%maxitr = 100
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
  end subroutine TofH_init
  
  real(r8) function TofH_eps (this)
    type(TofH), intent(in) :: this
    TofH_eps = this%root_finder%eps
  end function TofH_eps
  
  subroutine TofH_get_metrics (this, avg_itr, max_itr, rec_rate, avg_adj, max_adj)
    use parallel_communication, only: global_sum, global_maxval
    type(TofH), intent(in) :: this
    integer, intent(out), optional :: max_itr, max_adj
    real, intent(out), optional :: avg_itr, rec_rate, avg_adj
    if (present(avg_itr))  avg_itr  = real(global_sum(this%num_itr)) / max(1,global_sum(this%num_call))
    if (present(max_itr))  max_itr  = global_maxval(this%max_itr)
    if (present(rec_rate)) rec_rate = real(global_sum(this%num_rec)) / max(1,global_sum(this%num_call))
    if (present(avg_adj))  avg_adj  = real(global_sum(this%num_adj)) / max(1,global_sum(this%num_rec))
    if (present(max_adj))  max_adj  = global_maxval(this%max_adj)
  end subroutine TofH_get_metrics

!  subroutine TofH_compute (this, cell, H, Tmin, Tmax, T)
!    use truchas_logging_services, only: TLS_fatal
!    type(TofH), intent(inout) :: this
!    integer,  intent(in)  :: cell
!    real(r8), intent(in)  :: H, Tmin, Tmax
!    real(r8), intent(out) :: T
!    integer :: stat
!    character(80) :: errmsg
!    call set_context_cell_H (cell, H)
!    call ridder_find_root (this%root_finder, f, Tmin, Tmax, T, stat)
!    if (stat /= 0) then
!      if (stat < 0) then
!        write(errmsg,'(2(a,es21.14),a)') 'root not bracketed: [', Tmin, ',', Tmax, ']'
!      else if (stat > 0) then
!        write(errmsg,'(a,es10.4,2(a,es21.14))') &
!          'convergence failure: error=', this%root_finder%error, ', T=', T, ', H-H(T)=', f(T)
!      end if
!      call TLS_fatal ('TofH_compute: ' // trim(errmsg))
!    end if
!write(*,'(i0,1x)',advance='no') this%root_finder%numitr
!  end subroutine TofH_compute

  subroutine TofH_compute (this, cell, H, Tmin, Tmax, T)
    use truchas_logging_services, only: TLS_fatal
    type(TofH), intent(inout) :: this
    integer,  intent(in)  :: cell
    real(r8), intent(in)  :: H, Tmin, Tmax
    real(r8), intent(out) :: T
    integer :: n, stat
    character(100) :: errmsg
    real(r8) :: a, b, d
    call set_context_cell_H (cell, H) ! parameters used by f
    call ridders_find_root (this%root_finder, f, Tmin, Tmax, T, stat)
    if (stat < 0) then ! root not bracketed -- attempt to recover
      a = Tmin; b = Tmax; d = this%delta
      if (f(a) < 0.0_r8) then ! shift a to the left
        do n = 1, this%max_try
          a = a - d
          if (f(a) > 0.0_r8) exit
          d = 10*d
        end do
      else  ! then f(b) > 0; shift b to the right
        do n = 1, this%max_try
          b = b + d
          if (f(b) < 0.0_r8) exit
          d = 10*d
        end do
      end if
      if (n <= this%max_try) then ! root bracketed -- try again
        call ridders_find_root (this%root_finder, f, a, b, T, stat)
        if (stat == 0) then
          this%num_rec = this%num_rec + 1
          this%num_adj = this%num_adj + n
          this%max_adj = max(this%max_adj, n)
        end if
      end if
    end if
    if (stat == 0) then
      this%num_call = this%num_call + 1
      this%num_itr = this%num_itr + this%root_finder%numitr
      this%max_itr = max(this%max_itr, this%root_finder%numitr)
    else if (stat < 0) then
      write(errmsg,'(2(a,es21.14),a)') 'root not bracketed: [', a, ',', b, ']'
      call TLS_fatal ('TofH_compute: ' // trim(errmsg))
    else
      write(errmsg,'(a,es10.4,2(a,es21.14))') &
        'convergence failure: error=', this%root_finder%error, ', T=', T, ', H-H(T)=', f(T)
      call TLS_fatal ('TofH_compute: ' // trim(errmsg))
    end if
  end subroutine TofH_compute

end module TofH_type

  
  
