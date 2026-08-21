!!
!! HT_2D_TOFH_TYPE
!!
!! This module defines the HT_2D_TOFH type for computing cell temperatures
!! from enthalpy densities by inverting a cell material property mesh function
!! H(T).
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module ht_2d_tofh_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64, output_unit
  use new_mesh_func_class
  use ridders_class
  implicit none
  private

  !! Per-cell inverse of an increasing enthalpy-temperature relation. The
  !! enthalpy relation may differ by cell because material composition varies.
  type, extends(ridders), public :: ht_2d_tofh
    private
    class(new_mesh_func), pointer :: HofT => null()
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
  end type ht_2d_tofh

contains

  function f (this, x) result (fx)
    class(ht_2d_tofh), intent(in) :: this
    real(r8), intent(in) :: x
    real(r8) :: fx
    call this%HofT%compute_value_cell(this%cell, x, fx)
    fx = this%H - fx
  end function

  !! Initialize the inverse of H(T). EPS is the temperature root-finding
  !! tolerance. MAX_TRY and DELTA control recovery when COMPUTE receives an
  !! interval that does not bracket the root.

  subroutine init (this, HofT, eps, max_try, delta)
    class(ht_2d_tofh), intent(out) :: this
    class(new_mesh_func), target :: HofT
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

  !! Return collective root-finding and bracket-recovery performance metrics.
  subroutine get_metrics (this, avg_itr, max_itr, rec_rate, avg_adj, max_adj)
    use parallel_communication, only: global_sum, global_maxval
    class(ht_2d_tofh), intent(in) :: this
    integer, intent(out), optional :: max_itr, max_adj
    real, intent(out), optional :: avg_itr, rec_rate, avg_adj
    if (present(avg_itr))  avg_itr  = real(global_sum(this%num_itr)) / max(1,global_sum(this%num_call))
    if (present(max_itr))  max_itr  = global_maxval(this%max_itr)
    if (present(rec_rate)) rec_rate = real(global_sum(this%num_rec)) / max(1,global_sum(this%num_call))
    if (present(avg_adj))  avg_adj  = real(global_sum(this%num_adj)) / max(1,global_sum(this%num_rec))
    if (present(max_adj))  max_adj  = global_maxval(this%max_adj)
  end subroutine get_metrics

  !! Compute temperature T for cell CELL and enthalpy density H by finding
  !! the root of H - H(T). TMIN and TMAX should bracket the root. If they do
  !! not, the errant endpoint is shifted up to MAX_TRY times, starting with
  !! DELTA and increasing the shift by a factor of 10 after each attempt.
  !!
  !! H(T) is presently assumed to depend on temperature alone. Coupled
  !! thermal-species systems require an extended interface.

  subroutine compute (this, cell, H, Tmin, Tmax, T)
    class(ht_2d_tofh), intent(inout) :: this
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
      call panic('TofH_compute: ' // trim(errmsg))
    else
      write(errmsg,'(a,es10.4,2(a,es21.14))') &
        'convergence failure: error=', this%error, ', T=', T, ', H-H(T)=', this%f(T)
      call panic('TofH_compute: ' // trim(errmsg))
    end if
  end subroutine compute

  !! Terminate immediately after reporting an unrecoverable local failure.
  !! TofH evaluation is performed independently on each process, so this
  !! must not require a collective error path.

  subroutine panic(message)
    use parallel_communication, only: this_PE, abort_parallel_communication
    character(*), intent(in) :: message
    write(output_unit,'(a,i0,2a)') 'PANIC[', this_PE, ']: ', trim(message)
    flush(output_unit)
    call abort_parallel_communication
  end subroutine panic

end module ht_2d_tofh_type
