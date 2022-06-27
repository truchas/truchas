!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

MODULE DRIVERS
  !-----------------------------------------------------------------------------
  ! Purpose:
  !
  !    drivers to control the high level flow of control
  !
  ! Public Interface:
  !
  !    call CODE ()
  !
  ! Contains: CODE
  !           CYCLE_DRIVER
  !           CYCLE_INIT
  !           PROGRAM_SPECIFICATIONS
  !           PROCESS_COMMAND_LINE
  !
  ! Author(s): Bryan Lally (lally@lanl.gov)
  !            Douglas B. Kothe (dbk@lanl.gov)
  !-----------------------------------------------------------------------------
  use process_info_module
  implicit none
  private

  ! Public Procedures
  public :: CODE

  logical :: mem_on = .false.

CONTAINS

  SUBROUTINE CODE ()
    !---------------------------------------------------------------------------
    ! Truchas: A three-dimensional, high-resolution Eulerian algorithm
    ! designed to model incompressible flows of multiple immiscible
    ! fluids, coupled with heat transfer, phase change, and other physics. The
    ! fluids are delineated with tracked interfaces that can have surface
    ! tension.
    !---------------------------------------------------------------------------
    use input_driver,           only: READ_INPUT
    use truchas_env,            only: input_file, title
    use parallel_communication, only: init_parallel_communication
    use setup_module,           only: SETUP
    use signal_handler
    use output_utilities,       only: announce
    use truchas_logging_services
    use truchas_timers

    !---------------------------------------------------------------------------

    call init_parallel_communication

!   you can use this to debug in parallel under Linux with LAM and Totalview
!   see the comments in src/utility/wait_for_debugger.tv
!   call wait_for_debugger ()

    ! parse the command line
    call PROCESS_COMMAND_LINE

    ! arrange to catch signals
    call init_signal_handler(SIGUSR2)

    ! start overall timing
    call start_timer ("Total")

    ! initialize logging services
    call TLS_initialize

    ! announce
    call ANNOUNCE ('PROGRAM INFORMATION')
    call PROGRAM_SPECIFICATIONS ()

    ! read the data file
    call ANNOUNCE ('INPUT')
    call READ_INPUT (input_file, title)

    ! set up the problem
    call ANNOUNCE ('INITIALIZATION')
    call SETUP ()

! now go off and run my own code
!#define HIJACK
#ifdef HIJACK
call hijack_truchas ()
#else

    ! cycle through the problem
    call ANNOUNCE ('EXECUTION')
    call CYCLE_DRIVER ()
#endif

    ! prepare to terminate
    call ANNOUNCE ('TERMINATION')

    !Stop main timer
    call stop_timer("Total")


    ! Clean up
    call CLEANUP ()

    ! normal termination
    call TLS_exit

  END SUBROUTINE CODE

  SUBROUTINE CYCLE_DRIVER ()
    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !   Cycle through each time step
    !---------------------------------------------------------------------------
    use cycle_output_module,      only: CYCLE_OUTPUT_PRE, CYCLE_OUTPUT_POST
    use edit_module,              only: edit_short
    use EM,                       only: INDUCTION_HEATING
    use parallel_communication,   only: global_any
    use signal_handler
    use time_step_module,         only: cycle_number, cycle_max, dt, dt_old, t, t1, t2, dt_ds, &
        TIME_STEP, constant_dt, dt_constraint, dt_min
    use diffusion_solver,         only: ds_step, ds_accept, ds_restart, ds_get_face_temp_view, update_moving_vf
    use diffusion_solver_data,    only: ds_enabled
    use ustruc_driver,            only: ustruc_update
    use flow_driver, only: flow_enabled, flow_step, flow_accept, flow_vel_fn_view, &
        flow_set_pre_solidification_density
    use vtrack_driver, only: vtrack_update, vtrack_enabled, vtrack_vof_view, vtrack_flux_vol_view, &
        get_vof_from_matl
    use solid_mechanics_driver, only: solid_mechanics_enabled, solid_mechanics_step
    use string_utilities, only: i_to_c
    use truchas_danu_output, only: TDO_write_timestep
    use sim_event_queue_type
    use simulation_event_queue
    use time_step_sync_type
    use truchas_logging_services
    use truchas_timers
    use probes_driver, only: probes_write
    use kinds

    ! Local Variables
    Logical :: sig_rcvd, restart_ds
    integer :: c, errc, lookahead, num_try
    integer, parameter :: MAX_TRY = 10  !TODO: make expert parameter
    type(time_step_sync) :: ts_sync
    type(action_list), allocatable :: actions
    class(event_action), allocatable :: action
    real(r8), pointer :: vel_fn(:), vof(:,:), flux_vol(:,:), temperature_fc(:) => null()
    real(r8) :: tout, t_write, dt_new
    !---------------------------------------------------------------------------

    if (cycle_max == 0) then
      call TLS_info('')
      call TLS_info('Maximum number of cycles completed; writing time step data and terminating')
      return
    end if

    if (mem_on) call mem_diag_open

    call init_sim_event_queue(dt_min)
    call params%get('event-lookahead', lookahead, default=5)
    ts_sync = time_step_sync(lookahead)

    call start_timer('Main Cycle')

    call mem_diag_write('Before main loop:')

    call TDO_write_timestep
    t_write = t
    call probes_write(t)  ! Write initial probe info.

    t1 = t
    restart_ds = .false.
    call time_step  ! Does stuff not related to dt

    call event_queue%fast_forward(t)

    c = 0
    MAIN_CYCLE: do
      if (event_queue%is_empty()) exit

      tout = event_queue%next_time()

      do ! time step until reaching TOUT

        c = c + 1
        cycle_number = cycle_number + 1

        if (constant_dt) then
          t2 = t1 + dt
        else
          t2 = ts_sync%next_time(tout, t1, dt_old, dt) ! soft landing on TOUT
          if (t2 < t1 + dt) dt_constraint = 'time'
          dt = t2 - t1
        end if

        call mem_diag_write('Cycle ' // i_to_c(cycle_number) // ': Before output cycle:')

        call cycle_output_pre

        ! Evaluate the Joule heat source for the enthalpy calculation.
        call mem_diag_write('Cycle ' // i_to_c(cycle_number) // ': before induction heating:')
        call stop_timer('Main Cycle')
        call start_timer('electromagnetics')
        call induction_heating(t1, t2)
        call stop_timer('electromagnetics')
        call start_timer('Main Cycle')

        do num_try = 1, MAX_TRY

          ! move materials and associated quantities
          call mem_diag_write('Cycle ' // i_to_c(cycle_number) // ': before advection:')

          if (vtrack_enabled() .and. flow_enabled()) then
            vel_fn => flow_vel_fn_view()
            call vtrack_update(t, dt, vel_fn)
            vof => vtrack_vof_view()
            flux_vol => vtrack_flux_vol_view()
            call flow_set_pre_solidification_density(vof)
          end if

          ! solve heat transfer and phase change
          call mem_diag_write ('Cycle ' // i_to_c(cycle_number) // ': before heat transfer/species diffusion:')

          ! Diffusion solver: species concentration and/or heat.
          if (ds_enabled) then
            if (restart_ds) call ds_restart(t2 - t1)
            call ds_step(dt, dt_ds, errc)
            if (errc == 0) then
              call ds_accept
              exit
            end if
            dt = dt_ds
            t2 = t1 + dt
            call TLS_info('Diffusion solver step failed; retrying with reduced step size')
          else
            exit
          end if

        end do

        if (num_try > max_try) then
          call TLS_info('Too many repeated failures to take a step; giving up')
          exit MAIN_CYCLE
        end if

        ! calculate new velocity field
        call mem_diag_write('Cycle ' // i_to_c(cycle_number) // ': before fluid flow:')
        if (flow_enabled()) then
          if (ds_enabled) call ds_get_face_temp_view(temperature_fc)

          ! This updates the volume tracker's internal variable for the vof, so
          ! we can give the flow the current post-heat-transfer volume fractions.
          if (vtrack_enabled() .and. ds_enabled) call get_vof_from_matl(vof)

          call flow_step(t,dt,vof,flux_vol,temperature_fc)
          ! since this driver doesn't know any better, always accept
          call flow_accept()
        end if

        call mem_diag_write('Cycle ' // i_to_c(cycle_number) // ': before thermomechanics:')
        if (solid_mechanics_enabled()) call solid_mechanics_step(t, dt)

        ! output iteration information
        call cycle_output_post

        ! post-processing modules (no side effects)
        call mem_diag_write('Cycle ' // i_to_c(cycle_number) // ': before microstructure:')
        call ustruc_update(t2) ! microstructure modeling

        t = t2 ! set current time
        restart_ds = .false.

        call probes_write(t)

        ! set beginning cycle time (= previous cycle's end time)
        t1 = t2

        call time_step  ! next dt, plus other stuff it should not be doing

        ! See if a signal was caught.
        call read_signal(SIGUSR2, sig_rcvd)
        if (global_any(sig_rcvd)) then
          call TLS_info('')
          call TLS_info('Received signal USR2; writing time step data and terminating')
          call TDO_write_timestep
          exit MAIN_CYCLE
        end if

        if (t >= tout) exit
        if (c >= cycle_max) exit

      end do

      if (t >= tout) then
        ! Handle the event actions
        call event_queue%pop_actions(actions)
        dt_new = dt
        do
          call actions%get_next_action(action)
          if (.not.allocated(action)) exit
          select type (action)
          type is (output_event)
            call TDO_write_timestep
            t_write = t
          type is (short_edit_event)
            call edit_short
          type is (phase_event)
            dt_new = min(dt_new, action%init_dt(dt_old, dt))
            restart_ds = .true.
          type is (toolpath_event)
            call action%next_toolpath_segment
          type is (toolhead_event)
            call action%update_toolhead(t1)
            dt_new = min(dt_new, action%init_dt(dt_old, dt))
          type is (vf_event)
            call update_moving_vf
            dt_new = min(dt_new, action%init_dt(dt_old, dt))
          type is (stop_event)
            exit MAIN_CYCLE
          class default
            INSIST(.false.)
          end select
        end do
        dt = dt_new
      end if

      if (c >= cycle_max) then
        call TLS_info('')
        if (t == t_write) then
          call TLS_info('Maximum number of cycles completed; terminating')
        else
          call TLS_info('Maximum number of cycles completed; writing time step data and terminating')
          call TDO_write_timestep
        end if
        exit
      end if

    end do MAIN_CYCLE

    call stop_timer('Main Cycle')

  END SUBROUTINE CYCLE_DRIVER


  SUBROUTINE CLEANUP ()
    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !   clean up prior to termination, print reports
    !---------------------------------------------------------------------------
    use zone_module, only: zone_free
    use matl_module, only: matl_free
!NNC    use flow_driver, only: flow_destroy
    use time_step_module,       only: t, cycle_number
    use diffusion_solver,       only: ds_delete
    use truchas_logging_services
    use truchas_timers

    character(128) :: message

    !---------------------------------------------------------------------------

    !deallocate the fluidvof array, and others
!NNC    call flow_destroy()

    ! deallocate the base types
    call zone_free
    call matl_free

    ! free the diffusion solver resources
    call ds_delete ()

    ! end of run; print out diagnostics
    Write (message, 1) t, cycle_number
1   Format(17x,'Final Time: ',1pe11.4, ' after ', i5,' steps')
    call TLS_info (message)
    call TLS_info ('')

    ! report the timing info
    call write_timer_data

    call mem_diag_write ('before termination')
    call mem_diag_close

    ! report the timing info
    block
      use base_mesh_class
      use mesh_manager, only: named_mesh_ptr
      class(base_mesh), pointer :: mesh
      mesh => named_mesh_ptr('MAIN')
      call REPORT_MEMORY(mesh%cell_imap%global_size)
    end block

  END SUBROUTINE CLEANUP

  SUBROUTINE PROGRAM_SPECIFICATIONS ()
    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !   print build time and run time information
    !---------------------------------------------------------------------------
    use parallel_communication, only: nPE, IO_PE
    use utilities_module,     only: TIMESTAMP
    use truchas_logging_services
    use string_utilities, only: i_to_c
    use version_info

    character(LEN=48)  :: run_date, run_host, run_architecture
    character(:), allocatable :: version_str

    interface
      subroutine getrunhostinfo(n, arch, host) bind(c, name='getrunhostinfo')
        use,intrinsic :: iso_c_binding, only: c_char
        integer, value :: n
        character(kind=c_char), intent(out) :: arch(*), host(*)
      end subroutine
    end interface

    run_architecture = ""
    run_host         = ""
    call getrunhostinfo (len(run_architecture), run_architecture, run_host)
    call TIMESTAMP (run_date)

    call version(version_str)

    call TLS_info ('')
    call TLS_info ('   code:                ' // 'Truchas ' // version_str)
    call TLS_info ('   build architecture:  ' // ARCHITECTURE)
    call TLS_info ('   build date/time:     ' // BUILD_DATE)
    call TLS_info ('   build flags:         ' // COMPILER_FLAGS)
    call TLS_info ('   build host:          ' // HOST_NAME)
    call TLS_info ('   run architecture:    ' // run_architecture)
    call TLS_info ('   run host:            ' // run_host)
    call TLS_info ('   run date/time:       ' // run_date(5:22))
    if (nPE > 1) then
      call TLS_info ('   processors:          ' // i_to_c(nPE) // &
                     ' (processor ' // i_to_c(IO_PE) // ' is performing I/O)')
    else
      call TLS_info ('   processors:          ' // i_to_c(nPE))
    end if

  END SUBROUTINE PROGRAM_SPECIFICATIONS

  SUBROUTINE PROCESS_COMMAND_LINE
    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !    Scan the command line for arguments, and set flags accordingly.
    !
    !    A truchas command line looks like this:
    !
    !       truchas [-option -option] inputfile
    !
    !    Options are single letters.  They can be combined into one
    !    argument.  Accepted command line options are:
    !
    !       -q   quiet, eliminate most trace output, especially to the terminal
    !       -n   normal, write standard amount of trace output (default)
    !       -v   verbose, write lots of trace output
    !       -r   restart (not yet implemented here)
    !       -m   turn mem diagnostics off/on
    !
    !    We have two global integer variables, verbose and debug, that
    !    have values of 0, 1, 2.  0 is quiet, 1 is normal, 2 is
    !    verbose.  debug is 0 unless debug is explicitly turned on,
    !    then it takes the value 0, 1, or 2, according to the other
    !    flags.
    !
    !    if you change the command line options, please change the
    !    usage string declared below
    !
    !---------------------------------------------------------------------------

    use parallel_communication, only: is_IOP, halt_parallel_communication, broadcast
    use truchas_env,          only: input_dir, output_dir, prefix, &
                                    input_file, overwrite_output
    !use truchas_logging_services
    use utilities_module,     only: MAKE_DIRECTORY_HIERARCHY
    use parameter_module,     only: string_len
    use restart_variables,    only: restart_file, restart
    use file_utility,         only: count_tokens, get_token
    use string_utilities,     only: i_to_c
    use truchas_danu_output_data, only: io_group_size
    use truchas_logging_services, only: TLS_set_verbosity, TLS_VERB_SILENT, &
                                        TLS_VERB_NORMAL, TLS_VERB_NOISY

    ! local variables
    integer :: status, ios
    logical :: file_exist
    character (LEN=1024) :: string
    character (LEN=1024) :: token
    character (LEN=8)   :: file_read
    integer :: n_arg, n_tokens
    integer :: i
    integer :: indx
    integer :: j
    logical :: v
    logical :: o
    logical :: r
    logical :: h
    logical :: f
    logical :: g
    character (LEN=string_len), dimension(10) :: usage = (/                       &
    'usage: truchas [options] infile                                       ', &
    '                                                                      ', &
    'options:                                                              ', &
    '  -v:n          verbose level (0, 1, 2)                               ', &
    '  -o:filename   output filename root                                  ', &
    '  -r:filename   restart path/filename                                 ', &
    '  -m            turn on memory diagnostics                            ', &
    '  -g:n          output group size                                     ', &
    '  -f            force overwrite of output directory contents          ', &
    '  -h            help                                                  ' /)
    
    !---------------------------------------------------------------------------

    ! mark each flag as "unset"
    v = .false.                         ! verbose
    o = .false.                         ! output path
    r = .false.                         ! restart file
    h = .false.                         ! help
    f = .false.                         ! input file
    g = .false.                         ! h5 output group size
    overwrite_output = .false.
    ! there must be at least one argument
    n_arg = command_argument_count()
    if (n_arg < 1) then
       call TLS_error ('insufficient arguments')
       call error_check (.true., usage, 'PROCESS_COMMAND_LINE')
    end if

    ! process the arguments
    do i = 1, n_arg

       ! get the next argument
       call get_command_argument(i, string)

       ! arguments must start with '-', or be a filename (last argument)
       if (string(1:1) /= '-') then
          if (i /= n_arg) then
             call TLS_error ('filename must be last argument or argument must start with "-": ' // trim(string))
             call error_check (.true., usage, 'PROCESS_COMMAND_LINE')
          else
             input_file = string
             f = .true.
          end if

       else if (string == '--version') then

          block
             use version_info
             use,intrinsic :: iso_fortran_env, only: output_unit
             character(:), allocatable :: string
             if (is_IOP) then
                call version(string)
                write(output_unit,'(a)') string
             end if
             call halt_parallel_communication
             stop
          end block

       else

          ! arguments must have at least two characters
          if (LEN_TRIM(string) < 2) then
             call TLS_error ('invalid argument: ' // trim(string))
             call error_check (.true., usage, 'PROCESS_COMMAND_LINE')
          end if

          ! arguments of length > 2 must have a ':' at 3
          if (LEN_TRIM(string) > 2) then
             if (string(3:3) /= ':') then
                call TLS_error ('invalid argument: ' // trim(string))
                call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
             end if
          end if

          ! arguments can not have 3 characters (-c:)
          if (LEN_TRIM(string) == 3) then
             call TLS_error ('invalid argument (missing value): ' // TRIM(string))
             call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
          end if

          ! two tokens max (seperated by ':')
          n_tokens = COUNT_TOKENS(string,':')
          if (n_tokens > 2) then
             call TLS_error ('invalid argument: ' // trim(string))
             call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
          end if

          ! all tokens must have non-zero lengths
          do j = 1, n_tokens
             call GET_TOKEN(token,j,string,':')
             if (LEN_TRIM(token) == 0) then
                call TLS_error ('invalid argument: ' // trim(string))
                call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
             end if
          end do

          ! get the option token
          call GET_TOKEN(token,1,string,':')

          ! act on it
          select case (token(2:2))
          case ('m')
             ! Turn memory diagnostics on. - Added August 3rd, 2007.
             mem_on = .true.
          case ('g')
             if (g) then
                call TLS_error('repeated argument: ' // trim(string))
                call ERROR_CHECK(.true., usage, 'PROCESS_COMMAND_LINE')
             end if
             g = .true.
             if (n_tokens == 1) then
                call TLS_error('missing value for -g: ' // trim(string))
                call ERROR_CHECK(.true., usage, 'PROCESS_COMMAND_LINE')
             else
                call GET_TOKEN(token,2,string,':')
                read(token,*,iostat=ios) io_group_size
                if (ios /= 0) then
                   call TLS_error('invalid value for -g: ' // trim(token))
                   call ERROR_CHECK(.true., usage, 'PROCESS_COMMAND_LINE')
                else if (io_group_size < 0) then
                   call TLS_error('invalid value for -g: ' // i_to_c(io_group_size))
                   call ERROR_CHECK(.true., usage, 'PROCESS_COMMAND_LINE')
                end if
             end if
          case ('h')
             if (h) then
                call TLS_error ('repeated argument: ' // trim(string))
                call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
             end if
             h = .true.
          case ('v')
             if (v) then
                call TLS_error ('repeated argument: ' // trim(string))
                call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
             end if
             v = .true.
             if (n_tokens == 1) then
                call TLS_error ('invalid argument (missing verbosity level): ' // trim(string))
                call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
             else
                call GET_TOKEN(token,2,string,':')
                if (LEN_TRIM(token) /= 1) then
                   call TLS_error ('invalid argument: ' // trim(string))
                   call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
                end if
                select case (token(1:1))
                case ('0')
                   call TLS_set_verbosity(TLS_VERB_SILENT)
                case ('1')
                   call TLS_set_verbosity(TLS_VERB_NORMAL)
                case ('2')
                   call TLS_set_verbosity(TLS_VERB_NOISY)
                case default
                   call TLS_error ('invalid argument: ' // trim(string))
                   call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
                end select
             end if
          case ('o')
             if (o) then
                call TLS_error ('repeated argument: ' // trim(string))
                call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
             end if
             o = .true.
             if (n_tokens == 1) then
                call TLS_error ('invalid argument (pathname required): ' // trim(string))
                call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
             else
                call GET_TOKEN(prefix,2,string,':')
             end if
          case ('r')
             if (r) then
                call TLS_error ('repeated argument: ' // trim(string))
                call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
             end if
             r = .true.
             if (n_tokens == 1) then
                call TLS_error ('invalid argument (pathname required): ' // trim(string))
                call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
             else
                call GET_TOKEN(restart_file,2,string,':')
             end if
          case ('f')
             overwrite_output = .true.
          case default
             call TLS_error ('invalid argument ' // trim(string))
             call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
          end select
       end if
    end do

    ! act on what we found

    if (h) then
       call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
    end if

    if (.not. g) then
       io_group_size = 1
    end if

    if (.not. v) then
       call TLS_set_verbosity(TLS_VERB_NORMAL)
    end if

    if (.not. f) then
       call TLS_error ('no input file specified')
       call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
    end if

    ! fix up file name so that input_file ends in .inp
    indx = INDEX(input_file, '.inp', .true.)
    if (indx == 0 .or. indx+3 /= LEN_TRIM(input_file)) then
       input_file = Trim(input_file) // '.inp'
    end if

    ! check input_file for existence
    if (is_IOP) then
       inquire (FILE=Trim(input_file) , EXIST=file_exist)
    end if
    call broadcast (file_exist)
    if (.not. file_exist) then
       call TLS_error ('input file ' // trim(input_file) // ' does not exist')
       call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
    end if

    ! check input_file for readability
    if (is_IOP) then
       Inquire (FILE=TRIM(input_file) , READ=file_read)
    end if
    call broadcast (file_read)
    if (Trim(file_read) == 'NO') then
       call TLS_error ('input file ' // trim(input_file) // ' cannot be read')
       call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
    end if

    ! set restart flag
    if (r) then

       restart = .true.

       ! check restart_file for existence
       if (is_IOP) then
          inquire (FILE=Trim(restart_file) , EXIST=file_exist)
       end if
       call broadcast (file_exist)
       if (.not. file_exist) then
          call TLS_error ('restart file ' // trim(restart_file) // ' does not exist')
          call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
       end if

       ! check restart_file for readability
       if (is_IOP) then
          Inquire (FILE=TRIM(restart_file) , READ=file_read)
       end if
       call broadcast (file_read)
       if (Trim(file_read) == 'NO') then
          call TLS_error ('restart file ' // trim(restart_file) // ' cannot be read')
          call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
       end if
    else
       restart = .false.
    end if

    ! set prefix
    if (.not. o) then
       prefix = input_file(1:LEN_TRIM(input_file)-4)
       output_dir = TRIM(ADJUSTL(prefix)) // '_output/'
    else
       output_dir = trim(ADJUSTL(prefix))
       if (prefix(LEN(TRIM(prefix)):LEN_TRIM(prefix)) /= '/') then
          ! add a trailing slash if it was left out
          output_dir = trim(output_dir)//'/'
       end if
    end if

    ! Make directory hierarchy
    status=MAKE_DIRECTORY_HIERARCHY(trim(output_dir))

    indx = INDEX(input_file, '/', .true.)
    if(indx > 0) then
       input_dir = input_file(1:indx)
    else
       input_dir = './'
    end if

    if ( indx > 0 ) then
       prefix = TRIM(ADJUSTL(output_dir)) //  input_file(indx+1:LEN_TRIM(input_file)-4)
    else
       prefix = TRIM(ADJUSTL(output_dir)) //  input_file(1:LEN_TRIM(input_file)-4)
    end if

    return

  CONTAINS

    !! NNC, April 2012.  These are temporary stand-ins for ERROR_MODULE:ERROR_CHECK
    !! and TRUCHAS_LOGGING_SYSTEM:TLS_ERROR.  We can't use the new logging facilities
    !! because the output directory hasn't been created yet (that's one of the things
    !! this routine sets) and the logging can't be initialized until that's done
    !! (assuming we continue to insist upon logging to a file there, and not just stdout).
    !! This routine is processing the command line and any errors ought to be sent
    !! directly to stderr.  THIS ROUTINE NEEDS TO BE REFACTORED.

    subroutine ERROR_CHECK (flag, message, name)
      use,intrinsic :: iso_fortran_env, only: error_unit
      use parallel_communication, only: is_IOP, halt_parallel_communication
      logical, intent(in) :: flag ! ignored
      integer :: n
      character(*), intent(in) :: message(:)
      character(*), intent(in) :: name  ! ignored
      if (is_IOP) then
        do n = 1, size(message)
          write(error_unit,'(a)') message(n)(:len_trim(message(n)))
        end do
      end if
      call halt_parallel_communication
      stop
    end subroutine ERROR_CHECK

    subroutine TLS_error (message)
      use,intrinsic :: iso_fortran_env, only: error_unit
      use parallel_communication, only: is_IOP
      character(*), intent(in) :: message
      if (is_IOP) write(error_unit,'(a)') 'ERROR: ' // message(:len_trim(message))
    end subroutine TLS_error

    subroutine TLS_warn (message)
      use,intrinsic :: iso_fortran_env, only: error_unit
      use parallel_communication, only: is_IOP
      character(*), intent(in) :: message
      if (is_IOP) write(error_unit,'(a)') 'Warning: ' // message(:len_trim(message))
    end subroutine TLS_warn

  END SUBROUTINE PROCESS_COMMAND_LINE

END MODULE DRIVERS
