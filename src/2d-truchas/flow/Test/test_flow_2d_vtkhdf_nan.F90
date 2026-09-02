!!
!! TEST_FLOW_2D_VTKHDF_NAN
!!
!! This test verifies that quiet NaNs can be passed through the two-dimensional
!! flow VTKHDF writer without raising an IEEE invalid-operation exception.
!! The Python companion test checks that VTK reads the NaNs and excludes them
!! from finite data ranges.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, September 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

program test_flow_2d_vtkhdf_nan

  use, intrinsic :: iso_fortran_env, only: r8 => real64
  use mpi_f08
  use parameter_list_type
  use simulation_environment_type
  use unstr_2d_mesh_factory, only: new_unstr_2d_quad_mesh
  use unstr_2d_mesh_type, only: unstr_2d_mesh
  use flow_2d_vtkhdf_writer_type, only: flow_2d_vtkhdf_writer
  implicit none

  type(simulation_environment) :: env
  type(unstr_2d_mesh), pointer :: mesh
  type(flow_2d_vtkhdf_writer) :: output
  type(parameter_list) :: temporal_output
  real(r8), allocatable :: pressure(:), velocity(:,:)
  logical, allocatable :: flow_active(:)
  character(:), allocatable :: errmsg
  integer :: stat

  call MPI_Init(stat)
  env%comm = MPI_COMM_WORLD
  call MPI_Comm_rank(env%comm, env%rank)
  call MPI_Comm_size(env%comm, env%nproc)
  env%output_dir = '.'

  mesh => new_unstr_2d_quad_mesh(env, [0.0_r8, 0.0_r8], [2.0_r8, 1.0_r8], [2, 1], 0.0_r8)
  if (mesh%ncell /= 2) call fail('test mesh does not contain two cells')

  allocate(pressure(mesh%ncell), velocity(2, mesh%ncell))
  pressure = [1.0_r8, 0.0_r8]
  velocity(:,1) = [0.25_r8, 0.0_r8]
  velocity(:,2) = 0.0_r8
  flow_active = [.true., .false.]

  call output%open(env, mesh, temporal_output, stat, errmsg)
  if (stat /= 0) call fail(errmsg)
  call output%write_solution(0.0_r8, pressure, velocity, temporal_output, flow_active)
  call output%close()

  deallocate(mesh)
  call MPI_Finalize()

contains

  subroutine fail(message)
    character(*), intent(in) :: message
    if (env%rank == 0) write (*, '(2a)') 'FAIL: ', message
    call MPI_Abort(env%comm, 1)
  end subroutine fail

end program test_flow_2d_vtkhdf_nan
