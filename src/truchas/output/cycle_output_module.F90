!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CYCLE_OUTPUT_MODULE
  !=======================================================================
  ! Purpose:
  !
  !   Define procedures that perform output each cycle.
  !
  !   Public Interface(s):
  !
  !     * call CYCLE_OUTPUT_PRE ()
  !     * call CYCLE_OUTPUT_POST ()
  !
  !         Output key cycle information such as the time step, cycle
  !         number, and iteration count, to stdout and other output files.
  !         use _pre before the cycle computation, and _post after.
  !
  !     * call CYCLE_OUTPUT_DRIVER (quit)
  !
  !         Main driver for all cyclic output functions.
  !
  ! Contains: CYCLE_OUTPUT
  !           CYCLE_OUTPUT_DRIVER
  !
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov)
  !=======================================================================
  use truchas_logging_services
  use truchas_timers
  implicit none
  private

  public :: CYCLE_OUTPUT_PRE, CYCLE_OUTPUT_POST

CONTAINS

  SUBROUTINE CYCLE_OUTPUT_PRE ()
    !=======================================================================
    ! Purpose:
    !
    !   write cycle information that is known before the cycle begins
    !   (time step, cycle number) to stdout and various output files
    !=======================================================================
    use time_step_module, only: cycle_number, dt, dt_constraint, t

    character(128) :: string

    !! Log the time step and cycle number.
    write(string,20) cycle_number, t, trim(dt_constraint), dt
    20 format (1x,i10,': t = ',1pe13.5,', dt(',a,') = ',1pe13.5)
    call TLS_info ('')
    call TLS_info (string)

  END SUBROUTINE CYCLE_OUTPUT_PRE


  SUBROUTINE CYCLE_OUTPUT_POST ()
    !=======================================================================
    ! Purpose:
    !
    !   write cycle information that is known after the cycle ends
    !   (iteration counts) to stdout and various output files
    !=======================================================================
    use kinds
    use process_info_module,    only: get_process_size
    use parallel_communication
    use solid_mechanics_module, only: thermo_elastic_iterations, viscoplastic_iterations
    use physics_module, only: solid_mechanics => legacy_solid_mechanics
    use time_step_module,       only: cycle_number
    use flow_driver, only: flow_enabled, flow_vel_cc_view

    ! Local variables.
    integer :: vmsize, rssize, dsize
    character(128) :: string
    real(r8), pointer :: vel_cc(:,:)
    real(r8) :: x(3)

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
110 format ('SOLVER="',a,'" COUNT="',i10,'"')

    if (solid_mechanics) then
      write (string, 10) thermo_elastic_iterations,'Solid mechanics iterations (linear)'
      call TLS_info (string)
      write (string, 10) viscoplastic_iterations,'Solid mechanics iterations (nonlinear)'
      call TLS_info (string)
    end if
10  format (8x,i5,1x,a)

    ! Output diagnostics by physics
    if (flow_enabled()) then
      vel_cc => flow_vel_cc_view()
      x(1) = global_minval(vel_cc(1,:))
      x(2) = global_minval(vel_cc(2,:))
      x(3) = global_minval(vel_cc(3,:))
      call TLS_info ('')
      write(string, 99) x
99    format(12x,'Min Velocity: (',1p,e11.4,', ',e11.4,', ',e11.4,')')
      call TLS_info (string)
      x(1) = global_maxval(vel_cc(1,:))
      x(2) = global_maxval(vel_cc(2,:))
      x(3) = global_maxval(vel_cc(3,:))
      write(string, 98) x
98    format(12x,'Max Velocity: (',1p,e11.4,', ',e11.4,', ',e11.4,')')
      call TLS_info (string)
    endif

    ! If debug, write out additional memory usage info.
    if (TLS_verbosity >= TLS_VERB_NOISY) then
      call get_process_size (vmsize, rssize, dsize)
      if (vmsize /= -1) Then
        write (string, 20) global_maxval(vmsize), global_sum(vmsize)
20      format (8x,'vmsize, largest, total: ',i12,', ',i12,' kb')
        call TLS_info (string, TLS_VERB_NOISY)
      end if
    end if

  END SUBROUTINE CYCLE_OUTPUT_POST

END MODULE CYCLE_OUTPUT_MODULE
