!!
!! TRUCHAS_2D_MAIN
!!
!! This program drives the two-dimensional Truchas simulations.  The
!! simulation implementation is selected with the --simulation option.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, September 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

program truchas_2d_main

  use,intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use mpi_f08
  use parallel_communication
  use fhypre, only: fhypre_initialize
  use parameter_list_type
  use parameter_list_json
  use simulation_command_line_type
  use simulation_environment_type
  use simulation_factory
  use simulation_provenance
  use simulation_type
  implicit none

  integer :: inlun, stat
  character(:), allocatable :: errmsg
  type(parameter_list), pointer :: params
  type(simulation_command_line) :: cli
  type(simulation_environment) :: env
  class(simulation), allocatable :: sim

  call MPI_Init
  call init_parallel_communication
  call fhypre_initialize

  call cli%parse(stat, errmsg)
  if (cli%help) then
    if (is_IOP) then
      call cli%write_help('Two-dimensional Truchas simulation.')
      write(output_unit,'(a)') 'Simulations: flow, thermal, flow_thermal.'
    end if
    call MPI_Finalize
    stop
  end if
  if (stat == 0 .and. .not.allocated(cli%simulation)) then
    stat = 1
    errmsg = 'the --simulation option is required'
  end if
  if (stat /= 0) then
    if (is_IOP) write(error_unit,'(2a)') trim(cli%program) // ': ', errmsg
    call MPI_Finalize
    if (is_IOP) error stop 2
    stop
  end if

  call new_simulation(cli%simulation, sim, stat, errmsg)
  if (stat /= 0) then
    if (is_IOP) write(error_unit,'(2a)') trim(cli%program) // ': ', errmsg
    call MPI_Finalize
    if (is_IOP) error stop 2
    stop
  end if

  if (is_IOP) call cli%prepare_output_dir(stat, errmsg)
  call broadcast(stat)
  if (stat /= 0) then
    call broadcast_alloc_char(errmsg)
    if (is_IOP) write(error_unit,'(a)') trim(cli%program) // ': ' // errmsg
    call MPI_Finalize
    if (is_IOP) error stop 1
    stop
  end if

  env%input_dir = cli%input_dir
  env%output_dir = cli%output_dir
  env%comm = MPI_COMM_WORLD
  call MPI_Comm_rank(env%comm, env%rank)
  call MPI_Comm_size(env%comm, env%nproc)
  call env%simlog%init(env%comm, trim(env%output_dir)//'/run.log', stat, errmsg)
  if (stat /= 0) then
    if (env%rank == 0) write(error_unit,'(a)') 'error opening log file: ' // errmsg
    call MPI_Finalize
    error stop 1
  end if

  call write_simulation_prologue(env, cli%program, cli%simulation, cli%input_file, stat, errmsg)
  if (stat /= 0) then
    call env%simlog%error('error staging input file: ' // errmsg)
    call env%simlog%close
    call MPI_Finalize
    error stop 1
  end if

  open(newunit=inlun, file=cli%input_file, action='read', access='stream')
  call parameter_list_from_json_stream(inlun, params, errmsg)
  close(inlun)
  if (.not.associated(params)) then
    call env%simlog%error('error reading input file: ' // errmsg)
    call env%simlog%close
    call MPI_Finalize
    error stop 1
  end if

  call env%timer%start('simulation')
  call env%simlog%begin_section('Initializing simulation.')
  call sim%init(env, params, stat, errmsg)
  if (stat /= 0) then
    call env%simlog%end_section('Simulation initialization failed.')
    call env%simlog%error('simulation initialization error: ' // errmsg)
    call env%timer%stop('simulation')
    call env%simlog%close
    call MPI_Finalize
    error stop 1
  else
    call env%simlog%end_section('Simulation initialization complete.')
  end if

  call sim%run(env, stat, errmsg)
  call env%timer%stop('simulation')
  if (stat /= 0) call env%simlog%error('simulation error: ' // errmsg)

  call env%simlog%info('')
  call env%simlog%info('timing-summary-begin')
  if (env%rank == 0) then
    call env%timer%write(env%simlog%unit(), indent=2)
    if (env%simlog%terminal_output_enabled()) call env%timer%write(output_unit, indent=2)
  end if
  call env%simlog%info('timing-summary-end')

  call env%simlog%close
  call MPI_Finalize
  if (stat /= 0) error stop 1

end program truchas_2d_main
