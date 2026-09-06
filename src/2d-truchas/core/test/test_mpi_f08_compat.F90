program test_mpi_f08_compat

  use mpi_f08
  use mpi, only: legacy_mpi_comm_world => MPI_COMM_WORLD
  use parallel_communication, only: init_parallel_communication, halt_parallel_communication, &
      this_pe, npe, is_IOP, global_any
  implicit none

  type(MPI_Comm) :: duplicated_comm
  logical :: initialized, valid
  integer :: comparison, ierr, rank, nproc, status

  call init_parallel_communication
  status = 0

  call MPI_Initialized(initialized, ierr)
  valid = ierr == MPI_SUCCESS .and. initialized
  call require(valid, 'MPI was not initialized through the legacy interface')

  valid = MPI_COMM_WORLD%mpi_val == legacy_mpi_comm_world
  call require(valid, 'mpi_f08 MPI_COMM_WORLD does not wrap the legacy communicator value')

  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  valid = ierr == MPI_SUCCESS .and. rank == this_pe - 1
  call require(valid, 'mpi_f08 rank disagrees with parallel_communication')

  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  valid = ierr == MPI_SUCCESS .and. nproc == npe
  call require(valid, 'mpi_f08 size disagrees with parallel_communication')

  call MPI_Comm_dup(MPI_COMM_WORLD, duplicated_comm, ierr)
  valid = ierr == MPI_SUCCESS
  call require(valid, 'could not duplicate mpi_f08 MPI_COMM_WORLD')

  call MPI_Comm_compare(MPI_COMM_WORLD, duplicated_comm, comparison, ierr)
  valid = ierr == MPI_SUCCESS .and. comparison == MPI_CONGRUENT
  call require(valid, 'duplicated communicator is not congruent with mpi_f08 MPI_COMM_WORLD')

  call MPI_Comm_free(duplicated_comm, ierr)
  call require(ierr == MPI_SUCCESS, 'could not free duplicated mpi_f08 communicator')

  call halt_parallel_communication
  stop status

contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (global_any(.not.condition)) then
      if (is_IOP) write(*, '(2a)') 'ERROR: ', message
      status = 1
    end if
  end subroutine

end program test_mpi_f08_compat
