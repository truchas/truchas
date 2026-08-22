!!
!! SIMULATION_ENVIRONMENT_TYPE
!!
!! This module defines SIMULATION_ENVIRONMENT, the plain collection of
!! process-wide services owned and initialized by a simulation program. The
!! program initializes and destroys these components explicitly; this type
!! intentionally has no type-bound procedures. The timer tree is allocated by
!! the owning program and is referenced, but not owned, by lower-level code.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

module simulation_environment_type

  use mpi_f08, only: MPI_Comm
  use simulation_log_type, only: simulation_log
  use timer_tree_type, only: timer_tree
  implicit none
  private

  type, public :: simulation_environment
    type(MPI_Comm) :: comm
    integer :: rank = -1
    integer :: nproc = 0
    character(:), allocatable :: output_dir
    type(simulation_log) :: simlog
    type(timer_tree), pointer :: timer => null()
  end type simulation_environment

end module simulation_environment_type
