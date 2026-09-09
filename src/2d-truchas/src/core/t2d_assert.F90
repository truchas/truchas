!!
!!  T2D_ASSERT -- Assertions for the 2D Truchas implementation.
!!
!!  This procedure is the companion to the t2d_assert.inc include file,
!!  which defines the ASSERT, ASSERT_MSG, INSIST, and INSIST_MSG macros.  A
!!  failing assertion reports the MPI rank and aborts the MPI job when MPI is
!!  active.
!!
!!  Neil N. Carlson <neil.carlson@gmail.com>, September 2026
!!  SPDX-License-Identifier: BSD-3-Clause
!!
!!  Note: This procedure must not acquire any collective operations.  An
!!  assertion can be triggered by one process and is not necessarily called
!!  collectively.
!!

subroutine t2d_assert(file, line, msg)

  use,intrinsic :: iso_fortran_env, only: error_unit
  use mpi_f08

  character(*), intent(in) :: file
  integer, intent(in) :: line
  character(*), intent(in) :: msg

  logical :: mpi_is_initialized, mpi_is_finalized
  integer :: rank

  call MPI_Initialized(mpi_is_initialized)
  mpi_is_finalized = .false.
  if (mpi_is_initialized) call MPI_Finalized(mpi_is_finalized)

  if (mpi_is_initialized .and. .not.mpi_is_finalized) then
    call MPI_Comm_rank(MPI_COMM_WORLD, rank)
    if (len_trim(msg) > 0) then
      write(error_unit,'(a,":",i0,": T2D assertion failed on rank ",i0,": ",a)') &
        trim(file), line, rank, trim(msg)
    else
      write(error_unit,'(a,":",i0,": T2D assertion failed on rank ",i0)') trim(file), line, rank
    end if
    flush(error_unit)
    call MPI_Abort(MPI_COMM_WORLD, 1)
  else
    if (len_trim(msg) > 0) then
      write(error_unit,'(a,":",i0,": T2D assertion failed: ",a)') trim(file), line, trim(msg)
    else
      write(error_unit,'(a,":",i0,": T2D assertion failed")') trim(file), line
    end if
    flush(error_unit)
  end if

  error stop, quiet=.true.

end subroutine t2d_assert
