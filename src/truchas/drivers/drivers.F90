!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  !           PREPROCESSOR_DEFS
  !
  ! Author(s): Bryan Lally (lally@lanl.gov)
  !            Douglas B. Kothe (dbk@lanl.gov)
  !-----------------------------------------------------------------------------
  use parameter_module,  only: string_len
  use process_info_module
  implicit none
  private

  ! Public Procedures
  public :: CODE

  ! Disclaimer Notice
  character (LEN=string_len), dimension(2), target :: disclaimer = (/       &
     '   This Truchas release is registered with the Los Alamos National ', &
     '   Laboratory (LANL) as Los Alamos Computer Code LA-CC-15-097.     ' /)

  ! Copyright Notice
  character (LEN=string_len), dimension(39),target :: copyright = [                       &
     '   Copyright (c) 2007-2015. Los Alamos National Security, LLC.                   ', &
     '   All rights reserved.                                                          ', &
     '                                                                                 ', &
     '   This software was produced under U.S. Government contract DE-AC52-06NA25396   ', &
     '   for Los Alamos National Laboratory (LANL), which is operated by Los Alamos    ', &
     '   National Security, LLC for the U.S. Department of Energy. The U.S. Government ', &
     '   has rights to use, reproduce, and distribute this software.  NEITHER THE      ', &
     '   GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS  ', &
     '   OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE. If software', &
     '   is modified to produce derivative works, such modified software should be     ', &
     '   clearly marked, so as not to confuse it with the version available from LANL. ', &
     '                                                                                 ', &
     '   Additionally, redistribution and use in source and binary forms, with or      ', &
     '   without modification, are permitted provided that the following conditions    ', &
     '   are met:                                                                      ', &
     '                                                                                 ', &
     '   1. Redistributions of source code must retain the above copyright notice,     ', &
     '      this list of conditions and the following disclaimer.                      ', &
     '                                                                                 ', &
     '   2. Redistributions in binary form must reproduce the above copyright notice,  ', &
     '      this list of conditions and the following disclaimer in the documentation  ', &
     '      and/or other materials provided with the distribution.                     ', &
     '                                                                                 ', &
     '   3. Neither the name of Los Alamos National Security, LLC, Los Alamos National ', &
     '      Laboratory, LANL, the U.S. Government, nor the names of its contributors   ', &
     '      may be used to endorse or promote products derived from this software      ', &
     '      without specific prior written permission.                                 ', &
     '                                                                                 ', &
     '   THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND            ', &
     '   CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,        ', &
     '   BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS     ', &
     '   FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS         ', &
     '   NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,    ', &
     '   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT      ', &
     '   NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,     ', &
     '   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY         ', &
     '   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT           ', &
     '   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF      ', &
     '   THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.             ']

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
    use parallel_util_module,   only: PARALLEL_INIT
    use parallel_communication, only: init_parallel_communication
    use setup_module,           only: SETUP
    use timing_tree
    use truchas_timing
    use random_module,          only: INITIALIZE_RANDOM
    use signal_module,          only: SignalSet
    use pgslib_module,          only: PGSLib_CL_MAX_TOKEN_LENGTH
    use output_utilities,       only: announce
    use truchas_logging_services
    use truchas_danu_output, only: TDO_open, TDO_close
    implicit none

    character(len=PGSLib_CL_MAX_TOKEN_LENGTH), dimension(:), pointer :: argv
    integer :: i

    !---------------------------------------------------------------------------

    i = 0
    ! assign preprocessor definitions to global variables
    call PREPROCESSOR_DEFS ()

    ! initialize parallelism
    call PARALLEL_INIT (argv)
    call init_parallel_communication ()

!   you can use this to debug in parallel under Linux with LAM and Totalview
!   see the comments in src/utility/wait_for_debugger.tv
!   call wait_for_debugger ()

    ! parse the command line
    call PROCESS_COMMAND_LINE (argv)

    ! arrange to catch signals
    call SignalSet ()

    ! start overall timing
    call start_timer ("Total")

    ! initialize logging services
    call TLS_initialize

    ! announce
    call ANNOUNCE ('PROGRAM INFORMATION')
    call PROGRAM_SPECIFICATIONS ()

    call ANNOUNCE ('COPYRIGHT')
    call TLS_info (copyright)

    call ANNOUNCE ('DISCLAIMER')
    call TLS_info (disclaimer)

    ! open the danu output file
    call TDO_open

    ! initialize the random number generator
    call INITIALIZE_RANDOM()
    
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
    
    ! close the danu output file
    call TDO_close

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
    use advection_module,         only: ADVECT_MASS
    use cycle_output_module,      only: CYCLE_OUTPUT_PRE, CYCLE_OUTPUT_POST, &
                                        CYCLE_OUTPUT_DRIVER
    use fluid_flow_module,        only: FLUID_FLOW_DRIVER
    use EM,                       only: INDUCTION_HEATING
    use solid_mechanics_module,   only: THERMO_MECHANICS
    use pgslib_module,            only: PGSLib_GLOBAL_ANY
    use restart_variables,        only: restart
    use signal_module,            only: SignalInquire
    use time_step_module,         only: cycle_number, cycle_max, dt, t, t1, t2, dt_ds, &
                                        TIME_STEP
    use timing_tree
    use diffusion_solver,         only: ds_step
    use diffusion_solver_data,    only: ds_enabled
    use truchas_logging_services
    use string_utilities, only: i_to_c
    use truchas_danu_output, only: TDO_write_timestep
    use probe_output_module, only: probe_init_danu

    ! Local Variables
    Logical :: quit = .False.
    integer :: errc
    Integer :: c
    Integer :: HUP                      ! signal flag
    Integer :: USR2                     ! signal flag
    Integer :: URG                      ! signal flag

    !---------------------------------------------------------------------------
    
    if (mem_on) call mem_diag_open

    call PROBE_INIT_DANU  ! for Tbrook output this was done in TBU_writebasicdata

    ! start the cycle timer
    call start_timer ("Main Cycle")

    ! Prepass to initialize a solenoidal velocity field for Advection
    if(.not.restart) call Fluid_Flow_Driver (t)
  
    call mem_diag_write ('Before main loop:')

    ! Main computation loop.
    MAIN_CYCLE: do c = 1, cycle_max+1
       
       ! See if a signal was caught.
       call SignalInquire (HUP, USR2, URG)
       ! signal actions
       if (PGSLib_Global_Any(HUP  /= 0)) call TLS_info ('received signal HUP')
       if (PGSLib_Global_Any(USR2 /= 0)) call TLS_info ('received signal USR2')
       if (PGSLib_Global_Any(URG  /= 0)) then
          call TLS_info ('')
          call TLS_info ('received signal URG, writing timestep data and terminating')
          call TDO_write_timestep
          exit MAIN_CYCLE
       end if

       ! set current time
       t = t2

       ! Reset any time dependent conditions for the present cycle
       call Cycle_Init()

       ! perform any necessary cyclic output and check for termination
       call CYCLE_OUTPUT_DRIVER (quit, c)

       ! check for termination; exit if time to quit
       if (quit) exit MAIN_CYCLE

       ! increment time step counter
       cycle_number = cycle_number + 1
       ! set beginning cycle time (= previous cycle's end time)
       t1 = t2

       ! get the new time step
       call TIME_STEP ()

       call mem_diag_write ('Cycle ' // i_to_c(cycle_number) // ': Before output cycle:')

       ! output cycle time step
       call CYCLE_OUTPUT_PRE ()

       ! set ending cycle time (= beginning cycle time + dt)
       t2 = t1 + dt

       ! call each physics package in turn

       ! Evaluate the Joule heat source for the enthalpy calculation.
       call mem_diag_write ('Cycle ' // i_to_c(cycle_number) // ': before induction heating:')
       call INDUCTION_HEATING (t1, t2)
      
       ! move materials and associated quantities
       call mem_diag_write ('Cycle ' // i_to_c(cycle_number) // ': before advection:')
       call ADVECT_MASS ()

       ! solve heat transfer and phase change
       call mem_diag_write ('Cycle ' // i_to_c(cycle_number) // ': before heat transfer/species diffusion:')

       ! Diffusion solver: species concentration and/or heat.
       if (ds_enabled) then
         call ds_step (dt, dt_ds, errc)
         if (errc /= 0) call TLS_fatal ('CYCLE_DRIVER: Diffusion Solver step failed')
         ! The step size may have been reduced.  This assumes all other physics
         ! is off, and will need to be redone when the diffusion solver is made
         ! co-operable with the rest of the physics.
         t2 = t1 + dt
       end if

       ! calculate new velocity field
       call mem_diag_write ('Cycle ' // i_to_c(cycle_number) // ': before fluid flow:')
       call FLUID_FLOW_DRIVER (t)

       call mem_diag_write ('Cycle ' // i_to_c(cycle_number) // ': before thermomechanics:')
       call THERMO_MECHANICS ()

       ! output iteration information
       call CYCLE_OUTPUT_POST ()
       
       ! post-processing modules (no side effects)
       call mem_diag_write ('Cycle ' // i_to_c(cycle_number) // ': before microstructure:')

    end do MAIN_CYCLE
 
    ! stop the main cycle timer
    call  stop_timer ("Main Cycle")

  END SUBROUTINE CYCLE_DRIVER

  SUBROUTINE CLEANUP ()
    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !   clean up prior to termination, print reports
    !---------------------------------------------------------------------------
    use base_types_A_module,    only: BASE_TYPES_A_DEALLOCATE
    use base_types_B_module,    only: MESH_VERTEX_DEALLOCATE, BASE_TYPES_B_DEALLOCATE
    use debug_control_data
    use fluid_utilities_module, only: FLUID_DEALLOCATE
    use mesh_module,            only: Mesh, Vertex
    use time_step_module,       only: t, cycle_number
    use timing_tree
    use truchas_timing
    use diffusion_solver,       only: ds_delete
    use truchas_logging_services
    
    character(128) :: message

    !---------------------------------------------------------------------------

    !deallocate the fluidvof array, and others
    call FLUID_DEALLOCATE()

    ! deallocate the mesh
    call MESH_VERTEX_DEALLOCATE (Mesh, Vertex)

    ! deallocate the base types
    call BASE_TYPES_A_DEALLOCATE ()
    call BASE_TYPES_B_DEALLOCATE ()

    ! free the diffusion solver resources
    call ds_delete ()

    ! end of run; print out diagnostics
    Write (message, 1) t, cycle_number
1   Format(17x,'Final Time: ',1pe11.4, ' after ', i5,' steps')
    call TLS_info (message)
    call TLS_info ('')

    ! report the timing info
    call report_tree_timing()
    
    call mem_diag_write ('before termination')
    call mem_diag_close

    ! report the timing info
    call REPORT_MEMORY

  END SUBROUTINE CLEANUP

  SUBROUTINE PROGRAM_SPECIFICATIONS ()
    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !   print build time and run time information
    !---------------------------------------------------------------------------
    use code_module,          only: code_name, code_version, libraries,                      &
                                    build_date, build_architecture, build_flags, build_host, &
                                    run_architecture, run_host
    use parallel_info_module, only: p_info
    use utilities_module,     only: TIMESTAMP
    use truchas_logging_services
    use string_utilities, only: i_to_c

    character(LEN=32)  :: run_date

    run_architecture = ""
    run_host         = ""
    call getrunhostinfo (run_architecture, run_host)
    call TIMESTAMP (run_date)

    call TLS_info ('   code:                ' // code_name)
    call TLS_info ('   version:             ' // code_version)
    call TLS_info ('   libraries:           ' // libraries)
    call TLS_info ('   build architecture:  ' // build_architecture)
    call TLS_info ('   build date/time:     ' // build_date)
    call TLS_info ('   build flags:         ' // build_flags)
    call TLS_info ('   build host:          ' // build_host)
    call TLS_info ('   run architecture:    ' // run_architecture)
    call TLS_info ('   run host:            ' // run_host)
    call TLS_info ('   run date/time:       ' // run_date(5:22))
    if (p_info%nPE > 1) then
      call TLS_info ('   processors:          ' // i_to_c(p_info%nPE) // &
                     ' (processor ' // i_to_c(p_info%IO_ROOT_PE) // ' is performing I/O)')
    else
      call TLS_info ('   processors:          ' // i_to_c(p_info%nPE))
    end if

  END SUBROUTINE PROGRAM_SPECIFICATIONS

  SUBROUTINE PROCESS_COMMAND_LINE (argv)
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
    !       -d   debug, verbose  debugging output, in addition to trace output
    !       -t   debug, run all tests, mostly quiet
    !       -r   restart (not yet implemented here)
    !       -m   turn mem diagnostics off/on
    ! 
    !    We have two global integer variables, verbose and debug, that
    !    have values of 0, 1, 2.  0 is quiet, 1 is normal, 2 is
    !    verbose.  debug is 0 unless debug is explicitly turned on,
    !    then it takes the value 0, 1, or 2, according to the other
    !    flags.  -q -d is allowed, but is potentially uninteresting.
    !
    !    if you change the command line options, please change the
    !    usage string declared below
    !
    !---------------------------------------------------------------------------

    use parallel_info_module, only: p_info
    use pgslib_module,        only: PGSLib_BCAST
    use truchas_env,          only: input_dir, output_dir, prefix, input_file
    !use truchas_logging_services
    use utilities_module,     only: MAKE_DIRECTORY_HIERARCHY
    use parameter_module,     only: string_len
    use debug_control_data
    use restart_variables,    only: restart_file, restart
    use file_utility,         only: count_tokens, get_token

    ! arguments
    character(*), dimension(:), pointer :: argv

    ! local variables
    integer :: status
    logical :: file_exist
    character (LEN=1024) :: string
    character (LEN=1024) :: token
    character (LEN=8)   :: file_read
    integer :: n_tokens
    integer :: i
    integer :: indx
    integer :: j
    logical :: v
    logical :: d
    logical :: o
    logical :: r
    logical :: h
    logical :: f
    logical :: m
    character (LEN=string_len), dimension(9) :: usage = (/                       &
       'usage: truchas [options] infile                                       ', &
       '                                                                      ', &
       'options:                                                              ', &
       '  -v[:n]        verbose level (0, 1, 2)                               ', &
       '  -d[:n]        debug level (0, 1, 2)                                 ', &
       '  -o:filename   output filename root                                  ', &
       '  -r:filename   restart path/filename                                 ', &
       '  -m            turn on memory diagnostics                            ', &
       '  -h            help                                                  ' /)

    !---------------------------------------------------------------------------

    ! mark each flag as "unset"
    v = .false.                         ! verbose
    d = .false.                         ! debug
    o = .false.                         ! output path
    r = .false.                         ! restart file
    h = .false.                         ! help
    f = .false.                         ! input file
    m = .false.                         ! mem diagnostics on/off
    ! there must be at least one argument
    if (SIZE(argv) < 2) then
       call TLS_error ('insufficient arguments')
       call error_check (.true., usage, 'PROCESS_COMMAND_LINE')
    end if

    ! process the arguments
    do i = LBOUND(argv,1)+1, UBOUND(argv,1)

       ! get the next argument
       string = TRIM(argv(i))

       ! arguments must start with '-', or be a filename (last argument)
       if (string(1:1) /= '-') then
          if (i /= UBOUND(argv,1)) then
             call TLS_error ('filename must be last argument or argument must start with "-": ' // trim(string))
             call error_check (.true., usage, 'PROCESS_COMMAND_LINE')
          else
             input_file = string
             f = .true.
          end if

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
          case ('h')
             if (h) then
                call TLS_error ('repeated argument: ' // trim(string))
                call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
             end if
             h = .true.
          case ('d')
             if (d) then
                call TLS_error ('repeated argument: ' // trim(string))
                call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
             end if
             d = .true.
             if (n_tokens == 1) then
                debug = DEBUG_DEFAULT_SET
             else
                call GET_TOKEN(token,2,string,':')
                if (LEN_TRIM(token) /= 1) then
                   call TLS_error ('invalid argument: ' // trim(string))
                   call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
                end if
                select case (token(1:1))
                case ('0')
                   debug=DEBUG_NONE
                case ('1')
                   debug=DEBUG_QUIET
                case ('2')
                   debug=DEBUG_NOISY
                case default
                   call TLS_error ('invalid argument: ' // trim(string))
                   call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
                end select
             end if
          case ('v')
             !! The code path that is executed when this option is enabled is
             !! seriously broken.  So until resources are found to fix it we
             !! are going to simply disable it -- NNC, 11/13/06.

#define DISABLE_VERBOSE_OPTION
#ifdef DISABLE_VERBOSE_OPTION
             call TLS_warn ('the -v argument has been disabled in this release.')
#else
             if (v) then
                call TLS_error ('repeated argument: ' // trim(string))
                call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
             end if
             v = .true.
             if (n_tokens == 1) then
                verbose = VERBOSE_DEFAULT_SET
             else
                call GET_TOKEN(token,2,string,':')
                if (LEN_TRIM(token) /= 1) then
                   call TLS_error ('invalid argument: ' // trim(string))
                   call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
                end if
                select case (token(1:1))
                case ('0')
                   verbose=VERBOSE_QUIET
                case ('1')
                   verbose=VERBOSE_NORMAL
                case ('2')
                   verbose=VERBOSE_NOISY
                case default
                   call TLS_error ('invalid argument: ' // trim(string))
                   call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
                end select
             end if
#endif
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

    if (.not. d) then
       debug = DEBUG_DEFAULT
    end if

    if (.not. v) then
       verbose = VERBOSE_DEFAULT
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
    if (p_info%IOP) then
       inquire (FILE=Trim(input_file) , EXIST=file_exist)
    end if
    call PGSLib_BCAST (file_exist)
    if (.not. file_exist) then
       call TLS_error ('input file ' // trim(input_file) // ' does not exist')
       call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
    end if

    ! check input_file for readability
    if (p_info%IOP) then
       Inquire (FILE=TRIM(input_file) , READ=file_read)
    end if
    call PGSLib_BCAST (file_read)
    if (Trim(file_read) == 'NO') then
       call TLS_error ('input file ' // trim(input_file) // ' cannot be read')
       call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
    end if

    ! set restart flag
    if (r) then

       restart = .true.

       ! check restart_file for existence
       if (p_info%IOP) then
          inquire (FILE=Trim(restart_file) , EXIST=file_exist)
       end if
       call PGSLib_BCAST (file_exist)
       if (.not. file_exist) then
          call TLS_error ('restart file ' // trim(restart_file) // ' does not exist')
          call ERROR_CHECK (.true., usage, 'PROCESS_COMMAND_LINE')
       end if

       ! check restart_file for readability
       if (p_info%IOP) then
          Inquire (FILE=TRIM(restart_file) , READ=file_read)
       end if
       call PGSLib_BCAST (file_read)
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
      use parallel_info_module, only: p_info
      use pgslib_module, only: pgslib_finalize
      logical, intent(in) :: flag ! ignored
      integer :: n
      character(*), intent(in) :: message(:)
      character(*), intent(in) :: name  ! ignored
      if (p_info%IOP) then
        do n = 1, size(message)
          write(error_unit,'(a)') message(n)(:len_trim(message(n)))
        end do
      end if
      call pgslib_finalize
      stop
    end subroutine ERROR_CHECK
    
    subroutine TLS_error (message)
      use,intrinsic :: iso_fortran_env, only: error_unit
      use parallel_info_module, only: p_info
      character(*), intent(in) :: message
      if (p_info%IOP) write(error_unit,'(a)') 'ERROR: ' // message(:len_trim(message))
    end subroutine TLS_error
    
    subroutine TLS_warn (message)
      use,intrinsic :: iso_fortran_env, only: error_unit
      use parallel_info_module, only: p_info
      character(*), intent(in) :: message
      if (p_info%IOP) write(error_unit,'(a)') 'Warning: ' // message(:len_trim(message))
    end subroutine TLS_warn

  END SUBROUTINE PROCESS_COMMAND_LINE

  SUBROUTINE PREPROCESSOR_DEFS ()
    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !    get preprocessor definitions into global variables
    !
    !    we need this before we can punt!
    !---------------------------------------------------------------------------
    use code_module, only: code_name, code_version, libraries,                      &
                           build_date, build_architecture, build_flags, build_host

    !---------------------------------------------------------------------------

    ! Get the preprocessor definitions into variables.
    code_name           = CODE_NAME
    code_version        = VERSION
    libraries           = UBIKSOLVE
    libraries           = TRIM(ADJUSTL(libraries)) // ', ' // PGSLIB
#ifdef USE_CHACO
    libraries           = TRIM(ADJUSTL(libraries)) // ', ' // CHACO
#endif
    build_date          = BUILD_DATE
    build_architecture  = &
      ARCHITECTURE
    build_flags         = FFLAGS_TRIMMED_1 // &
                          FFLAGS_TRIMMED_2
    build_host          = HOST_NAME

  END SUBROUTINE PREPROCESSOR_DEFS


  SUBROUTINE CYCLE_INIT()
  !----------------------------------------------------------------------------------------
  !
  !  Purpose:
  !   To initialize any time varying conditions for the computational cycle at the
  !    start of each cycle
  !
  !    Jim Sicilian, CCS-2 (sicilian@lanl.gov)
  !    May 2003
  !
  !----------------------------------------------------------------------------------------
    use Update_BCS,  only: Update_Radiation_BC, Update_Dirichlet_BC, Update_HTC_External_BC
    use time_step_module, only: t

    call Update_Radiation_BC(t)
    call Update_Dirichlet_BC(t)
    call Update_HTC_External_BC(t)

  END SUBROUTINE CYCLE_INIT

END MODULE DRIVERS
