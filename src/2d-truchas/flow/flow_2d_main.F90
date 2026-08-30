!!
!! FLOW_2D_MAIN
!!
!! A basic JSON-input driver for the two-dimensional isothermal flow model.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

program flow_2d_main

  use,intrinsic :: iso_fortran_env, only: error_unit
  use mpi_f08
  use parallel_communication
  use fhypre, only: fhypre_initialize
  use parameter_list_type
  use parameter_list_json
  use simulation_command_line_type
  use simulation_environment_type
  use simulation_provenance
  use flow_2d_sim_type
  implicit none

  integer :: inlun, stat
  character(:), allocatable :: errmsg, sha256
  type(parameter_list), pointer :: params
  type(simulation_command_line) :: cli
  type(simulation_environment) :: env
  type(flow_2d_sim) :: sim

  call MPI_Init
  call init_parallel_communication

  call cli%parse(stat, errmsg)
  if (cli%help) then
    if (is_IOP) call cli%write_help('A JSON-input driver for the two-dimensional isothermal flow model.')
    call MPI_Finalize
    stop
  end if
  if (stat /= 0) then
    if (is_IOP) then
      write(error_unit,'(2a)') trim(cli%program) // ': ', errmsg
    end if
    call MPI_Finalize
    if (is_IOP) error stop 2
    stop
  end if

  call fhypre_initialize
  if (is_IOP) call cli%prepare_output_dir(stat, errmsg)
  call broadcast(stat)
  if (stat /= 0) then
    call broadcast_alloc_char(errmsg)
    if (is_IOP) then
      write(error_unit,'(a)') trim(cli%program) // ': ' // errmsg
    end if
    call MPI_Finalize
    if (is_IOP) error stop 1
    stop
  end if
  env%input_dir = cli%input_dir
  env%output_dir = cli%output_dir
  env%comm = MPI_COMM_WORLD
  call MPI_Comm_rank(env%comm, env%rank)
  call MPI_Comm_size(env%comm, env%nproc)
  allocate(env%timer)
  call env%simlog%init(env%comm, trim(env%output_dir)//'/run.log', stat, errmsg)
  if (stat /= 0) then
    if (env%rank == 0) write(error_unit,'(a)') 'error opening log file: ' // errmsg
    deallocate(env%timer)
    call MPI_Finalize
    error stop 1
  end if
  call write_program_prologue(env, cli%program)
  call stage_input_file(env, cli%input_file, 'input.json', sha256, stat, errmsg)
  if (stat /= 0) then
    call env%simlog%error('error staging input file: ' // errmsg)
    call env%simlog%close
    deallocate(env%timer)
    call MPI_Finalize
    error stop 1
  end if
  if (env%rank == 0) call env%simlog%info('initialization reading_input="input.json" sha256="' // sha256 // '"')
  open(newunit=inlun, file=cli%input_file, action='read', access='stream')
  call parameter_list_from_json_stream(inlun, params, errmsg)
  close(inlun)
  if (.not.associated(params)) then
    call env%simlog%error('error reading input file: ' // errmsg)
    call env%simlog%close
    deallocate(env%timer)
    call MPI_Finalize
    error stop 1
  end if
  call env%timer%start('simulation')
  call sim%init(env, params, stat, errmsg)
  if (stat == 0) call sim%run(env, stat, errmsg)
  call env%timer%stop('simulation')
  if (stat /= 0) then
    call env%simlog%error('flow simulation error: ' // errmsg)
    call env%simlog%close
    deallocate(env%timer)
    call MPI_Finalize
    error stop 1
  end if
  call env%simlog%info('')
  call env%simlog%info('Timing Summary:')
  call env%simlog%info('')
  if (env%rank == 0) call env%timer%write(env%simlog%unit(), indent=3)
  call env%simlog%close
  deallocate(env%timer)
  call MPI_Finalize

end program flow_2d_main
