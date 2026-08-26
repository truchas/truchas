!!
!! NS_HT_2D_SOLVER_TYPE
!!
!! This module defines NS_HT_2D_SOLVER, which coordinates an attempted
!! two-dimensional incompressible Navier--Stokes/thermal-transport step. It
!! advances material transport, converts its fluxes to an enthalpy rate,
!! attempts thermal transport, and then advances flow momentum and pressure.
!! The mesh, material composition, and models remain sim-owned. The flow
!! solver owns the flow state and provides accessors for data needed by the
!! coupled simulation driver.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module ns_ht_2d_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use simulation_environment_type
  use parameter_list_type
  use material_model_type
  use material_class
  use material_composition_type
  use flow_2d_model_type
  use flow_2d_material_transport_type
  use ns_2d_solver_type
  use ht_2d_model_type
  use ht_2d_solver_type
  use ns_ht_2d_enthalpy_advector_type
  use time_step_sync_type
  implicit none
  private

  type, public :: ns_ht_2d_solver
    private
    type(material_composition), pointer :: matl_comp => null() ! unowned reference
    integer, allocatable :: flow_material_ids(:)
    type(ns_2d_solver) :: flow
    type(flow_2d_material_transport) :: material_transport
    type(ns_ht_2d_enthalpy_advector) :: enthalpy_advector
    type(ht_2d_solver), pointer :: thermal => null()
    real(r8), allocatable :: temp(:), enthalpy_increment(:)
    integer :: ncell_onP
    real(r8) :: dt_init, dt_min, dt_max, dt_grow, courant_number, hnext, hlast
    integer :: max_try
    type(time_step_sync) :: ts_sync
  contains
    procedure :: init
    procedure :: set_initial_state
    procedure :: integrate
    procedure :: last_time
    procedure :: get_cell_flow_soln
    procedure :: get_face_velocity
    procedure :: get_cell_heat_soln
    procedure :: get_cell_temp_soln
    procedure :: write_metrics
    final :: delete
  end type

contains

  subroutine init(this, env, flow_model, ht_model, matl_model, matl_comp, &
      params, stat, errmsg)

    class(ns_ht_2d_solver), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(flow_2d_model), target, intent(in) :: flow_model
    type(ht_2d_model), target, intent(in) :: ht_model
    type(material_model), intent(in) :: matl_model
    type(material_composition), target, intent(in) :: matl_comp
    type(parameter_list), target, intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    type(parameter_list), pointer :: flow_params, momentum_params, projection_params, thermal_params
    class(material), pointer :: matl

    stat = 0
    ASSERT(size(matl_comp%vfrac,1) == matl_model%nmatl)
    if (matl_model%nmatl /= 1 .or. matl_model%nphase_real /= 1) then
      stat = 1
      errmsg = 'current non-isothermal flow requires one single-phase material'
      return
    end if
    call matl_model%get_matl_ref(1, matl)
    if (.not.matl%has_attr('fluid')) then
      stat = 1
      errmsg = 'current non-isothermal flow material must have the fluid attribute'
      return
    end if
    allocate(this%flow_material_ids(1), source=1)
    if (size(this%flow_material_ids) /= size(flow_model%density)) then
      stat = 1
      errmsg = 'flow material properties do not match fluid materials'
      return
    end if

    if (.not.params%is_sublist('flow') .or. .not.params%is_sublist('thermal')) then
      stat = 1
      errmsg = 'solver requires flow and thermal sublists'
      return
    end if
    flow_params => params%sublist('flow')
    thermal_params => params%sublist('thermal')
    call params%get('initial-time-step', this%dt_init, stat, errmsg)
    if (stat /= 0) return
    call params%get('min-time-step', this%dt_min, stat, errmsg)
    if (stat /= 0) return
    call params%get('max-time-step', this%dt_max, default=huge(1.0_r8), stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call params%get('time-step-growth', this%dt_grow, default=1.05_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call params%get('courant-number', this%courant_number, default=0.5_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call params%get('max-try-at-step', this%max_try, default=10, stat=stat, errmsg=errmsg)
    if (stat /= 0) return

    if (this%dt_init <= 0.0_r8 .or. this%dt_min <= 0.0_r8 .or. this%dt_min > this%dt_init .or. &
        this%dt_init > this%dt_max .or. this%dt_grow < 1.0_r8 .or. this%courant_number <= 0.0_r8 .or. &
        this%courant_number > 1.0_r8) then
      stat = 1
      errmsg = 'invalid coupled time-step controls'
      return
    end if
    if (this%max_try <= 0) then
      stat = 1
      errmsg = 'maximum coupled step attempts must be > 0'
      return
    end if
    if (.not.flow_params%is_sublist('projection-solver')) then
      stat = 1
      errmsg = 'solver.flow requires a projection-solver sublist'
      return
    end if
    projection_params => flow_params%sublist('projection-solver')

    this%ncell_onP = flow_model%mesh%ncell_onP
    this%matl_comp => matl_comp
    this%ts_sync = time_step_sync(4)
    allocate(this%temp(flow_model%mesh%ncell_onP), this%enthalpy_increment(flow_model%mesh%ncell_onP))
    if (flow_model%inviscid) then
      call this%flow%init(env, flow_model, projection_params=projection_params)
    else
      if (.not.flow_params%is_sublist('momentum-solver')) then
        stat = 1
        errmsg = 'viscous flow requires a momentum-solver sublist'
        return
      end if
      momentum_params => flow_params%sublist('momentum-solver')
      call this%flow%init(env, flow_model, momentum_params, projection_params)
    end if
    call this%material_transport%init(env, flow_model%mesh, size(this%flow_material_ids))
    if (allocated(ht_model%bc_inflow)) then
      call this%enthalpy_advector%init(flow_model%mesh, matl_model, this%flow_material_ids, stat, errmsg, &
          inflow_temperature=ht_model%bc_inflow)
    else
      call this%enthalpy_advector%init(flow_model%mesh, matl_model, this%flow_material_ids, stat, errmsg)
    end if
    if (stat /= 0) return
    allocate(this%thermal)
    call this%thermal%init(env, ht_model, thermal_params, stat, errmsg)
  end subroutine


  subroutine delete(this)
    type(ns_ht_2d_solver), intent(inout) :: this

    if (associated(this%thermal)) deallocate(this%thermal)
  end subroutine


  !! Set the initial flow and thermal states at TIME. The configured initial
  !! time step is supplied to their respective initial-condition procedures.
  subroutine set_initial_state(this, env, time, velocity, temp, stat, errmsg)

    class(ns_ht_2d_solver), intent(inout) :: this
    type(simulation_environment), intent(in) :: env
    real(r8), intent(in) :: time, velocity(:,:), temp(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    call this%thermal%set_initial_state(env, time, this%dt_init, temp, stat, errmsg)
    if (stat /= 0) return
    call this%thermal%get_cell_temp_soln(this%temp)
    call this%flow%set_buoyancy_temperature(this%temp)
    call this%flow%set_initial_state(time, this%dt_init, velocity, stat)
    if (stat /= 0) then
      errmsg = 'initializing flow state failed'
      return
    end if
    this%hnext = this%dt_init
    this%hlast = this%dt_init
  end subroutine


  !! Integrate the coupled state from its current time to TOUT. This owns
  !! output-time synchronization, thermal retry, time-step growth, and the
  !! convective Courant restriction.
  subroutine integrate(this, tout, stat, errmsg)
    use signal_handler, only: read_signal, SIGURG

    class(ns_ht_2d_solver), intent(inout) :: this
    real(r8), intent(in) :: tout
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8) :: t_n, t_np1, t, thermal_hnext
    logical :: sig_rcvd

    stat = 0
    t_n = this%thermal%last_time()
    ASSERT(tout >= t_n)
    do while (t_n < tout)
      t_np1 = this%ts_sync%next_time(tout, t_n, this%hlast, this%hnext)
      call attempt_step(this, t_n, t_np1, t, thermal_hnext, stat, errmsg)
      if (stat /= 0) return
      this%hlast = t - t_n
      this%hnext = min(thermal_hnext, this%dt_grow*this%hlast, this%dt_max, &
          this%flow%courant_time_step(this%courant_number))
      t_n = t
      call read_signal(SIGURG, sig_rcvd)
      if (sig_rcvd) then
        stat = 1
        errmsg = 'received SIGURG signal'
        return
      end if
    end do
  end subroutine


  !! Attempt a coupled step from T_N toward T_NP1. T is the accepted endpoint:
  !! it normally equals T_NP1, but is smaller after thermal recovery. A flow
  !! failure is non-recoverable.
  subroutine attempt_step(this, t_n, t_np1, t, hnext, stat, errmsg)

    class(ns_ht_2d_solver), intent(inout) :: this
    real(r8), intent(in) :: t_n, t_np1
    real(r8), intent(out) :: t, hnext
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    real(r8) :: dt, t_try
    real(r8), pointer :: face_velocity(:)

    ASSERT(t_np1 > t_n)
    ASSERT(this%thermal%last_time() == t_n)
    t = t_n
    t_try = t_np1
    do n = 1, this%max_try
      dt = t_try - t_n
      call this%flow%get_face_velocity(face_velocity)
      call this%material_transport%advance(t_n, t_try, face_velocity)
      call this%thermal%get_cell_temp_soln(this%temp)
      call this%enthalpy_advector%get_advected_enthalpy(t_n, this%temp, &
          this%material_transport%flux_volumes, this%enthalpy_increment)
      call this%thermal%set_ext_enthalpy_rate(this%enthalpy_increment/dt)
      call this%thermal%step(t_try, hnext, stat)
      if (stat /= 0) then
        t_try = t_n + hnext
        if (t_try - t_n < this%dt_min) then
          stat = -1
          errmsg = 'next coupled time step is too small'
          return
        end if
        cycle
      end if
      call this%thermal%get_cell_temp_soln(this%temp)
      call this%flow%set_buoyancy_temperature(this%temp)
      call this%flow%advance_momentum(t_n, t_try, this%material_transport%flux_volumes, stat, errmsg)
      if (stat /= 0) then
        stat = -3
        if (.not.allocated(errmsg)) errmsg = 'flow momentum update failed'
        return
      end if
      call this%thermal%commit_step
      t = t_try
      if (stat == 0) exit
    end do

    if (stat == 0) return

    stat = -2
    errmsg = 'unable to take a coupled time step'
  end subroutine


  real(r8) function last_time(this)
    class(ns_ht_2d_solver), intent(in) :: this

    last_time = this%thermal%last_time()
  end function


  !! Returns the current local cell pressure and velocity, including ghosts.
  subroutine get_cell_flow_soln(this, pressure, velocity)
    class(ns_ht_2d_solver), target, intent(in) :: this
    real(r8), pointer, intent(out) :: pressure(:), velocity(:,:)

    call this%flow%get_cell_flow_soln(pressure, velocity)
  end subroutine


  !! Returns the current face-normal velocity, including ghost faces.
  subroutine get_face_velocity(this, velocity)
    class(ns_ht_2d_solver), target, intent(in) :: this
    real(r8), pointer, intent(out) :: velocity(:)

    call this%flow%get_face_velocity(velocity)
  end subroutine


  subroutine get_cell_heat_soln(this, enth)
    class(ns_ht_2d_solver), intent(in) :: this
    real(r8), intent(inout) :: enth(:)

    call this%thermal%get_cell_heat_soln(enth)
  end subroutine


  subroutine get_cell_temp_soln(this, temp)
    class(ns_ht_2d_solver), intent(in) :: this
    real(r8), intent(inout) :: temp(:)

    call this%thermal%get_cell_temp_soln(temp)
  end subroutine


  subroutine write_metrics(this, string)
    class(ns_ht_2d_solver), intent(in) :: this
    character(*), intent(out) :: string(:)

    call this%thermal%write_metrics(string)
  end subroutine


end module ns_ht_2d_solver_type
