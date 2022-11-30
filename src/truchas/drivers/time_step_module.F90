!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE TIME_STEP_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Time step and cycle information.
  !
  !   Public Interface(s):
  !
  !     * call TIME_STEP ()
  !
  !         Compute the new time step.
  !
  !
  ! Contains: TIME_STEP
  !           TIME_STEP_COURANT
  !           TIME_STEP_VISCOUS
  !
  !           TIME_STEP_DISTANCE
  !           TIME_STEP_DISTANCE_SQUARED
  !
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov)
  !
  !=======================================================================
  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  ! Public Subroutines
  public :: TIME_STEP

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Namelist Input Variables
  ! Current Simulation Time
  real(r8), save, public :: t = 0.0_r8

  ! Cycle Numbers
  integer, save, public :: cycle_max
  integer, save, public :: cycle_number

  ! Time Step Numbers
  real(r8), save, public :: dt_constant
  real(r8), save, public :: dt_init, dt_grow
  real(r8), save, public :: dt_max, dt_min

  ! Derived Quantities
  character(LEN = 80),      save, public :: dt_constraint
  logical, save, public :: constant_dt ! Constant dt flag

  integer,               save, public :: cycle_number_restart

  real(r8), save, public :: dt            ! Current time step
  real(r8), save, public :: dt_old        ! Previous time step
  real(r8), save, public :: dt_ds = huge(1.0d0)        ! diffusion solver time step limit
  real(r8), save, public :: t1, t2        ! Pre and Post-Cycle times
  real(r8), save, public :: dt_surften    ! surface tension time step limit

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE TIME_STEP ()
    !=======================================================================
    ! Purpose(s):
    !   Compute the time step. The time step is limited by constraints for
    !   fluid-flow (advection and viscosity) and heat-transfer (conduction).
    !=======================================================================
    use physics_module,           only: flow
    use restart_variables,        only: restart
    use zone_module,              only: Zone
    use diffusion_solver_data,    only: ds_enabled
    use flow_driver,              only: flow_timestep
    use truchas_logging_services
    use truchas_timers

    ! Local Variables
    integer :: n, s
    real(r8) :: dt_next, dt_growth, dt_flow
    character(:), allocatable :: dt_flow_constraint
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Start Time-Step Timer
    call start_timer("Time Step")

    ! Initialize all old time data.
    ! ZJIBBEN 5/16/22: This subroutine is called after fields have already been
    !   updated. The temp_old field is needed for output during the current cycle,
    !   so it is instead updated in ds_accept.
    do n = 1, 3
       Zone%Vc_Old(n) = Zone%Vc(n)
    end do
    Zone%Rho_Old      = Zone%Rho
    !Zone%Temp_Old     = Zone%Temp
    Zone%Enthalpy_Old = Zone%Enthalpy

    ! Initialize Variables
    dt_old = dt

    ! Time Step Limits
    ! For the next time step, take the minimum of all constraints.

    ! For the next time step, take the minimum of all constraints.
    dt_growth = dt_grow*dt
    dt_next = MIN(dt_growth, dt_max)

    ! Diffusion solver time step limit.  The DT_DS value was returned by
    ! the solver on last step with a huge default initialized value to start.
    if (ds_enabled) dt_next = min(dt_next, dt_ds)

    ! Fluid-Flow Time Step: This must be done here because the ENTHALPY module
    !                       may have changed the fluid properties after NAVIER_STOKES
    !                       completed
    if (flow) then
      call flow_timestep(dt_flow, dt_flow_constraint)
      dt_next = min(dt_next, dt_flow)
    endif

    ! Set a character string according to the constraint that has been activated
    if (dt_next == dt_growth) then
      dt_constraint = 'growth'
    else if (dt_next == dt_max) then
      dt_constraint = 'maximum'
    else if (dt_next == dt_ds) then
      dt_constraint = 'diffusion solver'
    else if (dt_next == dt_flow) then
      dt_constraint = dt_flow_constraint
    end if

    ! (Non) Constant Time Step
    ! Constant Time Step
    if (constant_dt) then

       if (dt_next < dt_constant) then
          write(message,10) dt_constant, TRIM(dt_constraint), dt_next
10        format ('Constant time step of ',1pe13.5,' > ', &
                  a,' time step constraint of ',1pe13.5)
          call TLS_warn (message)
       end if
       dt = dt_constant             ! Constant Time Step
       dt_constraint = 'constant'   ! Time Step Constraint

    ! Non-Constant Time Step
    else
       if (cycle_number == cycle_number_restart .and. .not.restart) then
          ! First cycle; use initial time step
          dt = dt_init
          dt_constraint = 'initial'
       else
          ! Non-First Cycle; use computed time step
          dt = dt_next
       end if

    end if

    ! Minimum Time Step
    if (dt < dt_min) then
      write(message,'(3a,es13.5,a)') 'Time step too small: dt(', trim(dt_constraint), &
                                     ') = ', dt, ' < dt_min'
      call TLS_fatal(message)
    end if

    call stop_timer("Time Step")

  end subroutine time_step

end module time_step_module
