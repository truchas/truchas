!!
!! HT_2D_SOLVER_TYPE
!!
!! This module defines the solver for a 2D thermal transport simulation. It
!! owns the preconditioner, correction norm, integrator adapter, integrator,
!! and current solution vector, and coordinates initialization, time stepping,
!! and access to thermal state data for its sim-owned model.
!!
!! David Neill-Asanza <davidhneill@gmail.com>, July 2020
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module ht_2d_solver_type

  use,intrinsic :: iso_fortran_env, only: int64, r8 => real64
  use ht_2d_model_type
  use ht_2d_precon_type
  use ht_2d_norm_type
  use ht_2d_ic_solver_type
  use ht_2d_vector_type
  use ht_2d_idaesol_model_type
  use new_idaesol_type
  use parameter_list_type
  use simulation_environment_type
  use time_step_sync_type
  implicit none
  private

  type, public :: ht_2d_solver
    private
    type(ht_2d_model), pointer :: model => null()   ! reference only -- do not own
    type(ht_2d_precon) :: precon
    type(ht_2d_norm) :: norm
    type(ht_2d_idaesol_model) :: integ_model
    type(idaesol) :: integ
    !! Pending/current state
    real(r8) :: t
    type(ht_2d_vector) :: u
    logical :: step_is_pending = .false.
    type(parameter_list) :: ic_params
    !! Time-step policy
    integer(int64) :: nstep = 0_int64
    real(r8) :: tlast, hlast, hnext
    real(r8) :: dt_init, dt_min, dt_max, dt_grow
    integer :: max_try
    character(8) :: hnext_cause = 'init'
    type(time_step_sync) :: ts_sync
    logical :: time_stepper_initialized = .false.
  contains
    procedure :: init
    procedure :: init_time_stepper
    procedure :: set_initial_state
    procedure :: integrate
    procedure :: step
    procedure :: commit_step
    procedure :: reject_step
    procedure :: last_time
    procedure :: num_steps
    procedure :: initial_time_step
    procedure :: init_temporal_output
    procedure :: set_temporal_output
    procedure :: get_cell_heat_soln
    procedure :: get_cell_temp_soln
    procedure :: write_metrics
    procedure :: set_ext_enthalpy_rate
  end type

contains

  subroutine init(this, env, model, params, stat, errmsg)

    class(ht_2d_solver), intent(out), target :: this
    type(simulation_environment), intent(in) :: env
    type(ht_2d_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: plist
    character(:), allocatable :: context
    real(r8) :: rel_tol
    integer :: max_itr

    this%model => model
    call this%u%init(this%model%mesh)

    !! Create the preconditioner
    context = 'processing ' // params%path() // ': '
    if (params%is_sublist('preconditioner')) then
      plist => params%sublist('preconditioner')
      call this%precon%init(this%model, plist, stat, errmsg)
      if (stat /= 0) return
    else
      stat = 1
      errmsg = context//'missing "preconditioner" sublist parameter'
      return
    end if

    !! Create the error norm
    if (params%is_sublist('error-norm')) then
      plist => params%sublist('error-norm')
      call this%norm%init(this%model, plist, stat, errmsg)
      if (stat /= 0) return
    else
      stat = 1
      errmsg = context//'missing "error-norm" sublist parameter'
      return
    end if

    !! Create the IDAESOL model
    call this%integ_model%init(this%model, this%precon, this%norm)

    !! Create the IDAESOL integrator
    if (params%is_sublist('integrator')) then
      plist => params%sublist('integrator')
      call this%integ%init(this%integ_model, plist, stat, errmsg)
      if (stat /= 0) return
    else
      stat = 1
      errmsg = context//'missing "integrator" sublist parameter'
      return
    end if

    call this%ic_params%set('rel-tol', 1.0e-6_r8)
    call this%ic_params%set('max-iter', 100)
    if (params%is_sublist('initial-condition')) then
      plist => params%sublist('initial-condition')
      call plist%get('rel-tol', rel_tol, default=1.0e-6_r8, stat=stat, errmsg=errmsg)
      if (stat /= 0) then
        errmsg = context // errmsg
        return
      else if (rel_tol <= 0.0_r8) then
        stat = 1
        errmsg = context // '"rel-tol" must be > 0.0'
        return
      end if
      call this%ic_params%set('rel-tol', rel_tol)
      call plist%get('max-iter', max_itr, default=100, stat=stat, errmsg=errmsg)
      if (stat /= 0) then
        errmsg = context // errmsg
        return
      else if (max_itr <= 0) then
        stat = 1
        errmsg = context // '"max-iter" must be > 0'
        return
      end if
      call this%ic_params%set('max-iter', max_itr)
    end if
    stat = 0

  end subroutine init


  subroutine set_initial_state(this, env, t, temp, stat, errmsg, dt)
    class(ht_2d_solver), intent(inout), target :: this
    type(simulation_environment), intent(in) :: env
    real(r8), intent(in) :: t, temp(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    real(r8), intent(in), optional :: dt
    type(ht_2d_ic_solver) :: ic
    type(ht_2d_vector) :: udot
    real(r8) :: dt_ic

    if (present(dt)) then
      dt_ic = dt
    else
      ASSERT(this%time_stepper_initialized)
      dt_ic = this%dt_init
    end if

    this%t = t
    call udot%init(this%u)
    call this%ic_params%set('dt', dt_ic)
    call ic%init(this%model, this%ic_params)
    call ic%compute(env, t, temp, this%u, udot, stat, errmsg)
    if (stat /= 0) return
    call this%integ%set_initial_state(t, this%u, udot)
    this%nstep = 0_int64
    if (this%time_stepper_initialized) then
      this%tlast = t
      this%hnext = this%dt_init
      this%hnext_cause = 'init'
      this%hlast = this%hnext
    end if
  end subroutine set_initial_state


  !! Initialize the solver-owned time-step policy from simulation controls.
  subroutine init_time_stepper(this, params, stat, errmsg)
    class(ht_2d_solver), intent(inout) :: this
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer :: lookahead

    call params%get('initial-time-step', this%dt_init, stat, errmsg)
    if (stat /= 0) return
    call params%get('min-time-step', this%dt_min, stat, errmsg)
    if (stat /= 0) return
    call params%get('max-time-step', this%dt_max, default=huge(1.0_r8), stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call params%get('time-step-growth', this%dt_grow, default=5.0_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call params%get('time-step-lookahead', lookahead, default=3, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call params%get('max-try-at-step', this%max_try, default=10, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (this%dt_init <= 0.0_r8 .or. this%dt_min <= 0.0_r8 .or. this%dt_min > this%dt_init .or. &
        this%dt_init > this%dt_max .or. this%dt_grow < 1.0_r8 .or. lookahead < 1 .or. this%max_try <= 0) then
      stat = 1
      errmsg = 'require 0 < min-time-step <= initial-time-step <= max-time-step, ' // &
          'time-step-growth >= 1, time-step-lookahead >= 1, and max-try-at-step > 0'
      return
    end if
    this%ts_sync = time_step_sync(lookahead)
    this%time_stepper_initialized = .true.
    stat = 0
  end subroutine init_time_stepper

  !! Returns the current integration time.

  real(r8) function last_time(this)
    class(ht_2d_solver), intent(in) :: this
    last_time = this%integ%last_time()
  end function

  !! Returns the current cell enthalpy solution.

  subroutine get_cell_heat_soln(this, enth)
    class(ht_2d_solver), intent(in) :: this
    real(r8), intent(inout) :: enth(:)
    ASSERT(size(enth) == this%model%mesh%ncell_onP)
    enth = this%u%hc(:this%model%mesh%ncell_onP)
  end subroutine

  !! Returns the current cell temperature solution.

  subroutine get_cell_temp_soln(this, temp)
    class(ht_2d_solver), intent(in) :: this
    real(r8), intent(inout) :: temp(:)
    ASSERT(size(temp) == this%model%mesh%ncell_onP)
    temp = this%u%tc(:this%model%mesh%ncell_onP)
  end subroutine

  subroutine write_metrics(this, string)
    class(ht_2d_solver), intent(in) :: this
    character(*), intent(out) :: string(:)
    ASSERT(size(string) == 2)
    call this%integ%write_metrics(string)
  end subroutine

  !! Set the cell-integrated external enthalpy rate used in thermal residual
  !! evaluation.
  subroutine set_ext_enthalpy_rate(this, enthalpy_rate)
    class(ht_2d_solver), intent(inout) :: this
    real(r8), intent(in) :: enthalpy_rate(:)

    call this%model%set_ext_enthalpy_rate(enthalpy_rate)
  end subroutine

  !! Integrate the thermal state from its current time to TOUT. The solver
  !! owns the adaptive time-step policy and retries failed thermal steps.
  subroutine integrate(this, env, tout, stat, errmsg)
    use signal_handler, only: read_signal, SIGURG

    class(ht_2d_solver), intent(inout) :: this
    type(simulation_environment), intent(inout) :: env
    real(r8), intent(in) :: tout
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8) :: t_n, t_np1, hproposed, hthermal
    logical :: sig_rcvd
    integer :: n
    character(8) :: cause

    stat = 0
    ASSERT(this%time_stepper_initialized)
    t_n = this%tlast
    ASSERT(tout >= t_n)
    do while (t_n < tout)
      hproposed = this%hnext
      cause = this%hnext_cause
      t_np1 = this%ts_sync%next_time(tout, t_n, this%hlast, this%hnext)
      if (t_np1 < t_n + hproposed) cause = 'output'
      if (t_np1 - t_n < this%dt_min) then
        stat = -1
        errmsg = 'next time step is too small'
        return
      end if
      do n = 1, this%max_try
        call this%step(env, t_n, t_np1, stat, errmsg, cause, n, hthermal)
        if (stat == 0) exit
        t_np1 = t_n + hthermal
        if (t_np1 - t_n < this%dt_min) then
          stat = -1
          errmsg = 'next time step is too small'
          return
        end if
      end do
      if (stat /= 0) then
        stat = -2
        errmsg = 'unable to take a thermal time step'
        return
      end if
      call this%commit_step()
      t_n = t_np1
      this%hnext = min(hthermal, this%dt_grow*this%hlast, this%dt_max)
      this%hnext_cause = 'thermal'
      if (this%hnext == this%dt_grow*this%hlast) this%hnext_cause = 'growth'
      if (this%hnext == this%dt_max) this%hnext_cause = 'max'

      call read_signal(SIGURG, sig_rcvd)
      if (sig_rcvd) then
        stat = 1
        errmsg = 'received SIGURG signal'
        return
      end if
    end do
  end subroutine

  !! Attempt a step from the current committed state to time T. On success,
  !! the tentative solution is stored in THIS%U and remains pending until
  !! COMMIT_STEP is called; HNEXT is the suggested size of the next step. On
  !! failure, THIS%U is restored to the last committed state and STAT is
  !! nonzero.

  subroutine step(this, env, t_n, t_np1, stat, errmsg, step_cause, attempt, hnext)

    class(ht_2d_solver), intent(inout) :: this
    type(simulation_environment), intent(inout) :: env
    real(r8), intent(in) :: t_n, t_np1
    integer, intent(out) :: stat
    character(:), allocatable, optional, intent(out) :: errmsg
    character(*), optional, intent(in) :: step_cause
    integer, optional, intent(in) :: attempt
    real(r8), optional, intent(out) :: hnext

    character(256) :: line
    character(8) :: cause
    integer :: iattempt

    ASSERT(t_np1 > t_n)
    ASSERT(this%integ%last_time() == t_n)
    cause = 'explicit'
    if (present(step_cause)) cause = step_cause
    iattempt = 1
    if (present(attempt)) iattempt = attempt
    write(line,'(a,i0,a,i0,a,es0.5,a,es0.5,a,a)') 'step=', this%nstep + 1_int64, &
        ' attempt=', iattempt, ' t0=', t_n, ' dt=', t_np1 - t_n, ' cause=', trim(cause)
    call env%simlog%begin_section(trim(line))
    call env%timer%start('thermal/transport')

    call this%integ%step(t_np1, this%u, this%hnext, stat)
    call env%timer%stop('thermal/transport')
    call write_step_metrics()
    if (stat == 0) then
      this%t = t_np1
      this%step_is_pending = .true.
      call env%simlog%end_section('step-end status=accepted')
    else ! failed -- restore last good state before returning
      call this%integ%get_last_state_copy(this%u)
      this%t = t_n
      this%step_is_pending = .false.
      if (present(errmsg)) errmsg = 'thermal integrator step failed'
      call env%simlog%end_section('step-end status=failed')
    end if
    if (present(hnext)) hnext = this%hnext

  contains

    subroutine write_step_metrics()

      integer :: counters(6)
      character(128) :: line
      character(6) :: status

      call this%integ%get_stepping_statistics(counters)
      if (stat == 0) then
        status = 'ok'
      else
        status = 'failed'
      end if
      write(line,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,a)') 'thermal nres=', counters(1), &
          ' npc=', counters(2), ' nnr=', counters(4), ' nnf=', counters(5), &
          ' nsr=', counters(6), ' status=', trim(status)
      call env%simlog%info(trim(line))

    end subroutine write_step_metrics

  end subroutine step

  !! Commit the tentative solution produced by a successful STEP, making it the
  !! current state of the DAE system. This has no effect if no step is pending.

  subroutine commit_step(this)
    class(ht_2d_solver), intent(inout) :: this
    if (this%step_is_pending) then
      call this%integ%commit_state(this%t, this%u)
      if (this%time_stepper_initialized) then
        this%hlast = this%t - this%tlast
        this%tlast = this%t
      end if
      this%nstep = this%nstep + 1_int64
      this%step_is_pending = .false.
    end if
  end subroutine

  !! Reject the tentative solution produced by a successful STEP, restoring
  !! the last committed solution and time.
  subroutine reject_step(this)
    class(ht_2d_solver), intent(inout) :: this

    if (this%step_is_pending) then
      call this%integ%get_last_state_copy(this%u)
      this%t = this%integ%last_time()
      this%step_is_pending = .false.
    end if
  end subroutine


  integer(int64) function num_steps(this)
    class(ht_2d_solver), intent(in) :: this
    num_steps = this%nstep
  end function


  real(r8) function initial_time_step(this)
    class(ht_2d_solver), intent(in) :: this
    ASSERT(this%time_stepper_initialized)
    initial_time_step = this%dt_init
  end function


  subroutine init_temporal_output(this, data)
    class(ht_2d_solver), intent(in) :: this
    type(parameter_list), intent(inout) :: data
    call data%set('NStep', this%nstep)
  end subroutine


  subroutine set_temporal_output(this, data)
    class(ht_2d_solver), intent(in) :: this
    type(parameter_list), intent(inout) :: data
    call data%set('NStep', this%nstep)
  end subroutine

end module ht_2d_solver_type
