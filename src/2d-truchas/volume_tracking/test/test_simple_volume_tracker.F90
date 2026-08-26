!!
!! Direct unit tests for the two-dimensional simple volume tracker.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

program test_simple_volume_tracker

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mpi_f08
  use parallel_communication
  use simulation_environment_type
  use unstr_2d_mesh_type
  use unstr_2d_mesh_factory, only: new_unstr_2d_quad_mesh
  use simple_volume_tracker_type
  implicit none

  type(simulation_environment) :: env
  integer :: stat
  character(:), allocatable :: errmsg

  call MPI_Init
  call init_parallel_communication
  env%comm = MPI_COMM_WORLD
  call MPI_Comm_rank(env%comm, env%rank)
  call MPI_Comm_size(env%comm, env%nproc)
  call env%simlog%init(env%comm, 'test_simple_volume_tracker.log', stat, errmsg, terminal_output=.false.)
  if (stat /= 0) error stop 'initializing simulation log: ' // errmsg

  call test_immobile_solid(env)

  call env%simlog%close
  call MPI_Finalize

contains

  subroutine test_immobile_solid(env)

    type(simulation_environment), intent(in) :: env
    type(unstr_2d_mesh), pointer :: mesh
    type(simple_volume_tracker) :: tracker
    real(r8), allocatable :: vel(:), vof_n(:,:), vof(:,:), flux_vol(:,:), int_normal(:,:,:)
    real(r8) :: q
    integer :: f, j1, j2

    mesh => new_unstr_2d_quad_mesh(env, [0.0_r8, 0.0_r8], [2.0_r8, 1.0_r8], [2, 1])
    call tracker%init(env, mesh, 1, 1, 2, .false., [1,2])
    allocate(vel(size(mesh%cface)), vof_n(2,mesh%ncell), vof(2,mesh%ncell), &
        flux_vol(2,size(mesh%cface)), int_normal(2,2,mesh%ncell))

    f = 0
    do j1 = mesh%cstart(1), mesh%cstart(2)-1
      if (mesh%cnhbr(j1) == 2) then
        f = mesh%cface(j1)
        exit
      end if
    end do
    ASSERT(f > 0)
    j2 = mesh%cstart(2) - 1 + &
        findloc(mesh%cface(mesh%cstart(2):mesh%cstart(3)-1), f, dim=1)

    vof_n = 0.5_r8
    vel = 0.0_r8
    vel(j1) = 0.1_r8
    vel(j2) = -0.1_r8
    call tracker%flux_volumes(vel, vof_n, vof, flux_vol, int_normal, 1, 0, 0.1_r8)

    q = 0.5_r8 * 0.1_r8 * 0.1_r8 * mesh%area(f) / mesh%volume(1)
    call require_close(vof(1,:2), [0.5_r8-q, 0.5_r8+q], 'fluid transport')
    call require_close(vof(2,:2), [0.5_r8, 0.5_r8], 'immobile solid')
    call require_close(flux_vol(2,:), [(0.0_r8, j1=1,size(flux_vol,2))], 'solid flux')

    deallocate(mesh)

  end subroutine test_immobile_solid


  subroutine require_close(value, expected, label)

    real(r8), intent(in) :: value(:), expected(:)
    character(*), intent(in) :: label

    if (size(value) /= size(expected) .or. maxval(abs(value-expected)) > 1.0e-14_r8) &
        error stop 'FAIL: ' // label

  end subroutine require_close

end program test_simple_volume_tracker
