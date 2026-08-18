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

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use ht_2d_model_type
  use ht_2d_precon_type
  use ht_2d_norm_type
  use ht_2d_ic_solver_type
  use ht_2d_vector_type
  use ht_2d_idaesol_model_type
  use new_idaesol_type
  use parameter_list_type
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
  contains
    procedure :: init
    procedure :: set_initial_state
    procedure :: integrate
    procedure :: step
    procedure :: commit_step
    procedure :: last_time
    procedure :: get_cell_heat_soln
    procedure :: get_cell_temp_soln
    procedure :: write_metrics
    procedure :: set_integrator_log
  end type

contains

  subroutine init(this, model, params, stat, errmsg)

    class(ht_2d_solver), intent(out), target :: this
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


  subroutine set_initial_state(this, t, dt, temp, stat, errmsg)
    class(ht_2d_solver), intent(inout), target :: this
    real(r8), intent(in) :: t, dt, temp(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    type(ht_2d_ic_solver) :: ic
    type(ht_2d_vector) :: udot
    this%t = t
    call udot%init(this%u)
    call this%ic_params%set('dt', dt)
    call ic%init(this%model, this%ic_params)
    call ic%compute(t, temp, this%u, udot, stat, errmsg)
    if (stat /= 0) return
    call this%integ%set_initial_state(t, this%u, udot)
  end subroutine set_initial_state

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


  subroutine set_integrator_log(this, unit)
    class(ht_2d_solver), intent(inout) :: this
    integer, intent(in) :: unit
    call this%integ%set_verbose_stepping(unit)
  end subroutine

  !! This delegates to the IDAESOL integration driver.  A target time (TOUT)
  !! and/or (maximum) number of steps (NSTEP) is specified and the driver
  !! integrates until the target time or number of steps has been reached.
  !! The driver will adjust the time step as needed, and attempt to recover
  !! from failed steps by decreasing the time step if necessary.  The minimum
  !! and maximum step sizes (HMIN/HMAX) can be specified; if not, there is no
  !! limit.  The maximum number of attempts (MTRY) at a time step can also be
  !! specified; it defaults to a reasonable value.  The integration status is
  !! returned in STATUS; the possible values from IDAESOL are exported (see
  !! above).  The input value of HNEXT is the initial time step the driver
  !! will attempt to use.  Its return value is the time step the driver would
  !! use on the next step if it were continuing to integrate.  For the first
  !! call, HNEXT should be set to the (user-specified) initial time step, but
  !! thereafter the return value should normally be used for the next call.
  !! It permissible to change it, but there is little reason to do so in this
  !! multi-step driver scenario.

  subroutine integrate(this, hnext, status, nstep, tout, hmin, hmax, mtry)
    class(ht_2d_solver), intent(inout) :: this
    real(r8), intent(inout) :: hnext
    integer, intent(out) :: status
    integer,  intent(in), optional :: nstep, mtry
    real(r8), intent(in), optional :: tout, hmin, hmax
    call this%integ%integrate(hnext, status, nstep, tout, hmin, hmax, mtry)
    call this%integ%get_last_state_copy(this%u)
  end subroutine

  !! Attempt a step from the current committed state to time T. On success,
  !! the tentative solution is stored in THIS%U and remains pending until
  !! COMMIT_STEP is called; HNEXT is the suggested size of the next step. On
  !! failure, THIS%U is restored to the last committed state and STAT is
  !! nonzero.

  subroutine step(this, t, hnext, stat)

    class(ht_2d_solver), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(out) :: hnext
    integer, intent(out) :: stat

    call this%integ%step(t, this%u, hnext, stat)
    if (stat == 0) then
      this%t = t
      this%step_is_pending = .true.
    else ! failed -- restore last good state before returning
      call this%integ%get_last_state_copy(this%u)
      this%step_is_pending = .false.
    end if

  end subroutine step

  !! Commit the tentative solution produced by a successful STEP, making it the
  !! current state of the DAE system. This has no effect if no step is pending.

  subroutine commit_step(this)
    class(ht_2d_solver), intent(inout) :: this
    if (this%step_is_pending) then
      call this%integ%commit_state(this%t, this%u)
      this%step_is_pending = .false.
    end if
  end subroutine

end module ht_2d_solver_type
