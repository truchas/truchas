!!
!! VOF_2D_MAIN
!!
!! A JSON-input driver for the two-dimensional developer VOF simulation.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

program vof_2d_main

  use,intrinsic :: iso_fortran_env, only: error_unit
  use mpi_f08
  use parallel_communication
  use parameter_list_type
  use parameter_list_json
  use simulation_command_line_type
  use simulation_environment_type
  use simulation_provenance
  use vof_2d_sim_type
  implicit none

  integer :: inlun, stat
  character(:), allocatable :: errmsg
  type(parameter_list), pointer :: params
  type(simulation_command_line) :: cli
  type(simulation_environment) :: env
  type(vof_2d_sim) :: sim

  call MPI_Init
  call init_parallel_communication

  call cli%parse(stat, errmsg)
  if (cli%help) then
    if (is_IOP) call cli%write_help('A JSON-input driver for the two-dimensional developer VOF simulation.')
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
  call env%simlog%init(env%comm, trim(env%output_dir)//'/run.log', stat, errmsg)
  if (stat /= 0) then
    if (env%rank == 0) write(error_unit,'(a)') 'error opening log file: ' // errmsg
    call MPI_Finalize
    error stop 1
  end if
  call write_simulation_prologue(env, cli%program, cli%program, cli%input_file, stat, errmsg)
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
  call sim%init(env, params, stat, errmsg)
  if (stat /= 0) then
    call env%simlog%error('error initializing VOF simulation: ' // errmsg)
    call env%simlog%close
    call MPI_Finalize
    error stop 1
  end if
  call sim%run(env, stat, errmsg)
  if (stat /= 0) then
    call env%simlog%error('VOF simulation error: ' // errmsg)
    call env%simlog%close
    call MPI_Finalize
    error stop 1
  end if
  call env%simlog%close
  call MPI_Finalize

end program vof_2d_main
