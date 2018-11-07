!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
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

  public :: CYCLE_OUTPUT_PRE, CYCLE_OUTPUT_POST, CYCLE_OUTPUT_DRIVER

CONTAINS

  SUBROUTINE CYCLE_OUTPUT_PRE ()
    !=======================================================================
    ! Purpose:
    !
    !   write cycle information that is known before the cycle begins
    !   (time step, cycle number) to stdout and various output files
    !=======================================================================
    use time_step_module, only: cycle_number, dt, dt_constraint, t

    integer :: iStatus
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
    use debug_control_data
    use fluid_data_module,      only: fluid_flow, minVel, maxVel
    use process_info_module,    only: get_process_size
    use pgslib_module,          only: PGSLIB_GLOBAL_MAXVAL, PGSLIB_GLOBAL_SUM
    use parallel_communication
    use projection_data_module, only: mac_projection_iterations, &
        prelim_projection_iterations
    use viscous_data_module,    only: inviscid,             &
        viscous_implicitness, &
        viscous_iterations,   &
        prelim_viscous_iterations
    use solid_mechanics_module, only: thermo_elastic_iterations, viscoplastic_iterations
    use solid_mechanics_input,  only: solid_mechanics
    use time_step_module,       only: cycle_number
    use fluid_utilities_module
    use flow_driver, only: flow_enabled, flow_vel_cc_view

    ! Local variables.
    integer :: vmsize, rssize, dsize
    character(128) :: string
    real(r8), pointer :: vel_cc(:,:)
    real(r8) :: x(3)

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
110 format ('SOLVER="',a,'" COUNT="',i10,'"')

    ! Count the linear and nonlinear iterations.
    if (fluid_flow) then
      if(cycle_number == 1) then
        write (string, 10) prelim_projection_iterations,'Preliminary projection iterations (linear)'
        call TLS_info (string)
        if (.not. inviscid .and. viscous_implicitness > 0) then
          write (string, 10) prelim_viscous_iterations,'Preliminary Viscous iterations (linear)'
          call TLS_info (string)
        end if
      endif
      write (string, 10) mac_projection_iterations,'Projection iterations (linear)'
      call TLS_info (string)
      if (.not. inviscid .and. viscous_implicitness > 0) then
        write (string, 10) viscous_iterations,'Viscous iterations (linear)'
        call TLS_info (string)
      end if
    end if
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
    else if(fluid_flow) then
      ! Output the min/max velocities
      call calcVelLimits()
      call TLS_info ('')
      write(string, 30) minVel(1), minVel(2), minVel(3)
30    format(12x,'Min Velocity: (',1p,e11.4,', ',e11.4,', ',e11.4,')')
      call TLS_info (string)
      write(string, 40) maxVel(1), maxVel(2), maxVel(3)
40    format(12x,'Max Velocity: (',1p,e11.4,', ',e11.4,', ',e11.4,')')
      call TLS_info (string)
    endif

    ! If debug, write out additional memory usage info.
    if (debug >= DEBUG_NOISY) then
      call get_process_size (vmsize, rssize, dsize)
      if (vmsize /= -1) Then
        write (string, 20) PGSLIB_GLOBAL_MAXVAL(vmsize), PGSLIB_GLOBAL_SUM(vmsize)
20      format (8x,'vmsize, largest, total: ',i12,', ',i12,' kb')
        call TLS_info (string)
      end if
    end if

  END SUBROUTINE CYCLE_OUTPUT_POST

  SUBROUTINE CYCLE_OUTPUT_DRIVER (quit, cycle)
    !=======================================================================
    ! Purpose:
    !
    !   Main driver for all cyclic output functions.
    !=======================================================================
    use edit_module,             only: EDIT_SHORT, short_edit, Short_Output_Dt_Multiplier
    use interface_output_module, only: Int_Output_Dt_Multiplier, interface_dump, &
                                       time_for_int_dump
    use output_control,          only: next_op, nops, Output_Dt, Output_T
    use time_step_module,        only: cycle_number, cycle_number_restart, t1, t2, &
                                       cycle_max, t
    use output_control,          only: Output_Dt_Multiplier, retain_last_step
    use diagnostics_module,      only: DIAGNOSTICS
    use probe_output_module,     only: PROBES_OUTPUT, Probe_Output_Cycle_Multiplier
    use truchas_danu_output,     only: TDO_write_timestep

    ! Argument List
    integer, intent(IN)    :: cycle
    logical, intent(INOUT) :: quit

    ! Local Variables
    integer :: idiff, last, next
    logical :: do_edit_short, do_xml_output

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Start the outputs timer
    call start_timer("Output")

    ! Calculate diagnostics quantities needed for output
    call DIAGNOSTICS

    ! Initialize dump flags
    time_for_int_dump  = .false.
    do_edit_short      = .true.
    do_xml_output      = .true.

    ! Initial time output
    INITIAL_OUTPUT: if (cycle_number  ==  cycle_number_restart) then

       ! Initial edits
       if (short_edit) then
          call EDIT_SHORT ()
          do_edit_short = .false.
       end if

       ! Initial output
       call TDO_write_timestep
       do_xml_output = .false.

       if (retain_last_step) then
          !TODO: need Danu version of "last step" output
       end if

       ! Interface dump is not yet integrated into main output
       ! Interfaces flag
       if (interface_dump) time_for_int_dump = .true.


    end if INITIAL_OUTPUT

    ! Current time output
10  continue

    ! Determine if we need to quit.
    if (t2 < Output_T(next_op)) then
       if (next_op == 1) then
          quit = .true.
          go to 30
       end if
       next_op = next_op - 1
       go to 10
    end if

15  continue

    if (t2 > Output_T(next_op+1)) then
       if (next_op >= nops) then
          go to 20
       end if
       next_op = next_op + 1
       go to 15
    end if

    ! If there is a new output specification since last time
    ! step and the user has not just started a new run,  check
    ! all output multipliers and perform all outputs as may be
    ! requested.  Perform outputs if either the new or the old
    ! frequency is non-zero.
    TIME_OUTPUT: if (t1 < Output_T(next_op) .and. &
                     cycle_number /= cycle_number_restart) then

       if (Short_Output_Dt_Multiplier(next_op) /= 0 .or. &
            Short_Output_Dt_Multiplier(MAX(next_op-1,1)) /= 0) then
          call EDIT_SHORT ()
          do_edit_short = .false.
       end if

       if (Output_Dt_Multiplier (next_op) /= 0 .or. &
            Output_Dt_Multiplier (MAX(next_op-1,1)) /= 0) then
          call TDO_write_timestep
          do_xml_output = .false.
       end if

       if (Int_Output_Dt_Multiplier(next_op) /= 0 .or. &
            Int_Output_Dt_Multiplier(MAX(next_op-1,1)) /= 0) then
          time_for_int_dump = .true.
       end if

       go to 25 ! Skip any other output checks

    end if TIME_OUTPUT

    ! If output cycle count between last and next (this) time step
    ! spans an output multiplier or more than an output multiplier,
    ! check all output frequencies and perform outputs as requested.
20  continue

    last = (t1 - Output_T(next_op))/Output_Dt(next_op)  ! Number op cycles for last time
    next = (t2 - Output_T(next_op))/Output_Dt(next_op)  ! Number op cycles for this time
    idiff = next - last

    MULTIPLIER_OUTPUT: if (idiff > 0 .and. cycle_number /= cycle_number_restart) then

       if (Short_Output_Dt_Multiplier(next_op) > 0) then
          if (MOD(next, Short_Output_Dt_Multiplier(next_op)) <= &
               MOD(last, Short_Output_Dt_Multiplier(next_op)) .or. &
               idiff >= Short_Output_Dt_Multiplier(next_op)) then
             call EDIT_SHORT ()
             do_edit_short = .false.
          end if
       end if

       if (Output_Dt_Multiplier(next_op) > 0) then
          if (MOD(next, Output_Dt_Multiplier(next_op)) <= &
               MOD(last, Output_Dt_Multiplier(next_op)) .or. &
               idiff >= Output_Dt_Multiplier(next_op)) then
             call TDO_write_timestep
             do_xml_output = .false.
          end if
       end if

       if (Int_Output_Dt_Multiplier(next_op) > 0) then
          if (MOD(next, Int_Output_Dt_Multiplier(next_op)) <= &
               MOD(last, Int_Output_Dt_Multiplier(next_op)) .or. &
               idiff >= Int_Output_Dt_Multiplier(next_op)) then
             time_for_int_dump = .true.
          end if
       end if


    end if MULTIPLIER_OUTPUT

    ! Check cycle outputs
    EVERY_CYCLE_OUTPUT: if (cycle_number /= cycle_number_restart) then

       if (Short_Output_Dt_Multiplier(next_op) < 0) then
          call EDIT_SHORT ()
          do_edit_short = .false.
       end if

       if (Output_Dt_Multiplier(next_op) < 0) then
          call TDO_write_timestep
          do_xml_output=.false.
       end if

       if (Probe_Output_Cycle_Multiplier < 0) then
          call PROBES_OUTPUT ()
       end if

       if (Int_Output_Dt_Multiplier(next_op) < 0) then
          time_for_int_dump = .true.
       end if

       if (retain_last_step) then
       end if


       if (Probe_Output_Cycle_Multiplier > 0) then
          if (MOD(cycle_number, Probe_Output_Cycle_Multiplier) <= &
               MOD(cycle_number-1, Probe_Output_Cycle_Multiplier)) then
             call PROBES_OUTPUT()
          end if
       end if

    end if EVERY_CYCLE_OUTPUT

25  continue

    ! Determine whether or not it's time to quit. If so, then perform
    ! those specified outputs if they haven't already been done.
    TERMINATION_OUTPUT: if (t2 > Output_T(next_op+1) .or. &
                            t2 == Output_T(nops+1)   .or. &
                            cycle > cycle_max) then

       ! Set the quit flag.
       quit = .true.

       ! Last editsout
       if (short_edit .and. do_edit_short) then
          call EDIT_SHORT ()
       end if

       if (do_xml_output) then
          call TDO_write_timestep
       end if

    end if TERMINATION_OUTPUT

30  continue

    ! Stop the outputs timer
    call stop_timer("Output")

  END SUBROUTINE CYCLE_OUTPUT_DRIVER

END MODULE CYCLE_OUTPUT_MODULE
