!!
!! T2D_THERMAL_SOLVER_TYPE
!!
!! This module defines the thermal solver for a 2D thermal transport
!! simulation. It owns the preconditioner, correction norm, IDAESOL adapter,
!! and current solution vector, and performs individual thermal steps for its
!! sim-owned model.
!!
!! David Neill-Asanza <davidhneill@gmail.com>, July 2020
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module t2d_thermal_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use ht_2d_model_type
  use ht_2d_precon_type
  use ht_2d_norm_type
  use ht_2d_ic_solver_type
  use ht_2d_vector_type
  use ht_2d_idaesol_model_type
  use new_idaesol_type
  use parameter_list_type
  use simulation_environment_type
  implicit none
  private

  type, public :: t2d_thermal_solver
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
    procedure :: step
    procedure :: commit_step
    procedure :: reject_step
    procedure :: last_time
    procedure :: get_cell_heat_soln
    procedure :: get_cell_temp_soln
    procedure :: write_metrics
    procedure :: set_ext_enthalpy_rate
  end type

contains

  subroutine init(this, env, model, params, stat, errmsg)

    class(t2d_thermal_solver), intent(out), target :: this
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
    class(t2d_thermal_solver), intent(inout), target :: this
    type(simulation_environment), intent(in) :: env
    real(r8), intent(in) :: t, temp(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    real(r8), intent(in) :: dt
    type(ht_2d_ic_solver) :: ic
    type(ht_2d_vector) :: udot
    real(r8) :: dt_ic

    dt_ic = dt

    this%t = t
    call udot%init(this%u)
    call this%ic_params%set('dt', dt_ic)
    call ic%init(this%model, this%ic_params)
    call ic%compute(env, t, temp, this%u, udot, stat, errmsg)
    if (stat /= 0) return
    call this%integ%set_initial_state(t, this%u, udot)
  end subroutine set_initial_state

  !! Returns the current integration time.

  real(r8) function last_time(this)
    class(t2d_thermal_solver), intent(in) :: this
    last_time = this%integ%last_time()
  end function

  !! Returns the current cell enthalpy solution.

  subroutine get_cell_heat_soln(this, enth)
    class(t2d_thermal_solver), intent(in) :: this
    real(r8), intent(inout) :: enth(:)
    ASSERT(size(enth) == this%model%mesh%ncell_onP)
    enth = this%u%hc(:this%model%mesh%ncell_onP)
  end subroutine

  !! Returns the current cell temperature solution.

  subroutine get_cell_temp_soln(this, temp)
    class(t2d_thermal_solver), intent(in) :: this
    real(r8), intent(inout) :: temp(:)
    ASSERT(size(temp) == this%model%mesh%ncell_onP)
    temp = this%u%tc(:this%model%mesh%ncell_onP)
  end subroutine

  subroutine write_metrics(this, string)
    class(t2d_thermal_solver), intent(in) :: this
    character(*), intent(out) :: string(:)
    ASSERT(size(string) == 2)
    call this%integ%write_metrics(string)
  end subroutine

  !! Set the cell-integrated external enthalpy rate used in thermal residual
  !! evaluation.
  subroutine set_ext_enthalpy_rate(this, enthalpy_rate)
    class(t2d_thermal_solver), intent(inout) :: this
    real(r8), intent(in) :: enthalpy_rate(:)

    call this%model%set_ext_enthalpy_rate(enthalpy_rate)
  end subroutine


  !! Attempt a step from the current committed state to time T. On success,
  !! the tentative solution is stored in THIS%U and remains pending until
  !! COMMIT_STEP is called; HNEXT is the suggested size of the next step. On
  !! failure, THIS%U is restored to the last committed state and STAT is
  !! nonzero.

  subroutine step(this, env, t_n, t_np1, stat, errmsg, hnext)

    class(t2d_thermal_solver), intent(inout) :: this
    type(simulation_environment), intent(inout) :: env
    real(r8), intent(in) :: t_n, t_np1
    integer, intent(out) :: stat
    character(:), allocatable, optional, intent(out) :: errmsg
    real(r8), optional, intent(out) :: hnext
    real(r8) :: hnext_local

    ASSERT(t_np1 > t_n)
    ASSERT(this%integ%last_time() == t_n)
    call env%timer%start('thermal/transport')

    call this%integ%step(t_np1, this%u, hnext_local, stat)
    call env%timer%stop('thermal/transport')
    call write_step_metrics()
    if (stat == 0) then
      this%t = t_np1
      this%step_is_pending = .true.
    else ! failed -- restore last good state before returning
      call this%integ%get_last_state_copy(this%u)
      this%t = t_n
      this%step_is_pending = .false.
      if (present(errmsg)) errmsg = 'thermal integrator step failed'
    end if
    if (present(hnext)) hnext = hnext_local

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
    class(t2d_thermal_solver), intent(inout) :: this
    if (this%step_is_pending) then
      call this%integ%commit_state(this%t, this%u)
      this%step_is_pending = .false.
    end if
  end subroutine

  !! Reject the tentative solution produced by a successful STEP, restoring
  !! the last committed solution and time.
  subroutine reject_step(this)
    class(t2d_thermal_solver), intent(inout) :: this

    if (this%step_is_pending) then
      call this%integ%get_last_state_copy(this%u)
      this%t = this%integ%last_time()
      this%step_is_pending = .false.
    end if
  end subroutine


end module t2d_thermal_solver_type
