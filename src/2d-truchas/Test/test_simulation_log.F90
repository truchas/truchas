program test_simulation_log

  use simulation_environment_type
  use simulation_log_type
  use mpi_f08
  use parallel_communication, only: init_parallel_communication, halt_parallel_communication, global_any
  implicit none

  type(simulation_environment) :: env
  character(:), allocatable :: errmsg
  integer :: ierr, stat
  logical :: valid

  call init_parallel_communication
  call MPI_Comm_dup(MPI_COMM_WORLD, env%comm, ierr)
  call require(ierr == MPI_SUCCESS, 'could not duplicate the simulation communicator')
  call MPI_Comm_rank(env%comm, env%rank, ierr)
  call require(ierr == MPI_SUCCESS, 'could not determine the simulation rank')
  call MPI_Comm_size(env%comm, env%nproc, ierr)
  call require(ierr == MPI_SUCCESS, 'could not determine the simulation size')

  call env%simlog%init(env%comm, 'simulation_log_test.log', stat, errmsg, verbosity=-1)
  call require(stat /= 0, 'invalid log verbosity was accepted')
  call require(errmsg == 'invalid log verbosity', 'invalid log verbosity message was not broadcast')

  call env%simlog%init(env%comm, 'simulation_log_test.log', stat, errmsg, LOG_NORMAL, .false.)
  call require(stat == 0, 'could not initialize normal log')
  valid = .true.
  if (env%rank == 0) valid = env%simlog%unit() /= 0
  call require(valid, 'I/O process has no log unit')
  call env%simlog%info('normal message')
  call env%simlog%begin_section('outer section')
  call env%simlog%info('nested message')
  call env%simlog%begin_section('detail section', LOG_DETAIL)
  call env%simlog%info('deeply nested message')
  call env%simlog%end_section('detail section complete', LOG_DETAIL)
  call env%simlog%warn('nested warning message')
  call env%simlog%end_section('outer section complete')
  call env%simlog%info('detail message', LOG_DETAIL)
  call env%simlog%warn('warning message')
  call env%simlog%error('error message')
  call env%simlog%close()

  call require(log_contains('simulation_log_test.log', 'normal message'), 'normal message missing')
  call require(log_contains('simulation_log_test.log', 'outer section'), 'section begin message missing')
  call require(log_contains('simulation_log_test.log', '  nested message'), 'nested message was not indented')
  call require(log_contains('simulation_log_test.log', '    deeply nested message'), 'deeply nested message was not indented')
  call require(.not.log_contains('simulation_log_test.log', 'detail section complete'), 'detail section end was not filtered')
  call require(log_contains('simulation_log_test.log', '  Warning: nested warning message'), 'nested warning was not indented')
  call require(log_contains('simulation_log_test.log', 'outer section complete'), 'section end message missing')
  call require(.not.log_contains('simulation_log_test.log', 'detail message'), 'detail message was not filtered')
  call require(log_contains('simulation_log_test.log', 'Warning: warning message'), 'warning message missing')
  call require(log_contains('simulation_log_test.log', 'ERROR: error message'), 'error message missing')

  call env%simlog%init(env%comm, 'simulation_log_test.log', stat, errmsg, LOG_DETAIL, .false.)
  call require(stat == 0, 'could not initialize detail log')
  call env%simlog%info('detail message', LOG_DETAIL, terminal=.false.)
  call env%simlog%close()
  call require(log_contains('simulation_log_test.log', 'detail message'), 'detail message missing')

  call MPI_Comm_free(env%comm, ierr)
  call require(ierr == MPI_SUCCESS, 'could not free the simulation communicator')
  call halt_parallel_communication

contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (global_any(.not.condition)) then
      if (env%rank == 0) write(*, '(2a)') 'ERROR: ', message
      error stop 1
    end if
  end subroutine


  logical function log_contains(filename, text) result(found)
    character(*), intent(in) :: filename, text

    character(256) :: line
    integer :: ios, unit

    if (env%rank == 0) then
      found = .false.
      open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
      if (ios == 0) then
        do
          read(unit, '(a)', iostat=ios) line
          if (ios /= 0) exit
          if (index(line, text) > 0) then
            found = .true.
            exit
          end if
        end do
        close(unit)
      end if
    end if
    call MPI_Bcast(found, 1, MPI_LOGICAL, 0, env%comm, ios)
  end function

end program test_simulation_log
