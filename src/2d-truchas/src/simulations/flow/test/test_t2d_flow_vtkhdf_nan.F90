!!
!! TEST_T2D_FLOW_VTKHDF_NAN
!!
!! This test verifies that quiet NaNs can be passed through the two-dimensional
!! flow VTKHDF writer without raising an IEEE invalid-operation exception.
!! The Python companion test checks that VTK reads the NaNs and excludes them
!! from finite data ranges.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, September 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

program test_t2d_flow_vtkhdf_nan

  use, intrinsic :: iso_fortran_env, only: r8 => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use mpi_f08
  use parallel_communication, only: init_parallel_communication, halt_parallel_communication
  use parameter_list_type
  use material_database_type
  use material_model_type
  use simulation_environment_type
  use t2d_unstr_mesh_factory, only: new_unstr_2d_quad_mesh
  use t2d_unstr_mesh_type, only: t2d_unstr_mesh
  use t2d_flow_vtkhdf_writer_type, only: t2d_flow_vtkhdf_writer
  implicit none

  type(simulation_environment) :: env
  type(t2d_unstr_mesh), pointer :: mesh
  type(material_database) :: database
  type(material_model) :: matl_model
  type(t2d_flow_vtkhdf_writer) :: output
  type(parameter_list) :: temporal_output
  character(1) :: no_materials(0)
  real(r8), allocatable :: pressure(:), velocity(:,:)
  real(r8), allocatable :: compliance(:), target_divergence(:)
  logical, allocatable :: flow_active(:)
  character(:), allocatable :: errmsg
  integer :: stat, c, id

  call init_parallel_communication
  env%comm = MPI_COMM_WORLD
  call MPI_Comm_rank(env%comm, env%rank)
  call MPI_Comm_size(env%comm, env%nproc)
  env%output_dir = '.'

  mesh => new_unstr_2d_quad_mesh(env, [0.0_r8, 0.0_r8], [4.0_r8, 2.0_r8], [4, 2], 0.0_r8)

  call matl_model%init(no_materials, database, stat, errmsg)
  if (stat /= 0) call fail('initializing material model: ' // errmsg)

  allocate(pressure(mesh%ncell), velocity(2, mesh%ncell), flow_active(mesh%ncell))
  allocate(compliance(mesh%ncell_onP), target_divergence(mesh%ncell_onP))
  pressure = 1.0_r8
  velocity(1,:) = 0.25_r8
  velocity(2,:) = 0.0_r8
  do c = 1, mesh%ncell
    id = mesh%cell_imap%global_index(c)
    flow_active(c) = mod(id,4) /= 0
    if (c > mesh%ncell_onP) cycle
    compliance(c) = 0.125_r8
    select case (mod(id,4))
    case (1)
      target_divergence(c) = -0.5_r8
    case (2)
      target_divergence(c) = 0.25_r8
    case (3)
      target_divergence(c) = 0.0_r8
    case (0)
      compliance(c) = 0.0_r8
      target_divergence(c) = ieee_value(0.0_r8, ieee_quiet_nan)
    end select
  end do

  call output%open(env, mesh, matl_model, temporal_output, stat, errmsg, compliance_output=.true.)
  if (stat /= 0) call fail(errmsg)
  call output%write_solution(0.0_r8, pressure, velocity, temporal_output, flow_active, &
      compliance=compliance, void_target_divergence=target_divergence)
  call output%write_solution(1.0_r8, pressure, velocity, temporal_output, flow_active)
  call output%close()

  ! The default writer omits both diagnostics when compliance is disabled.
  env%output_dir = 'disabled'
  call output%open(env, mesh, matl_model, temporal_output, stat, errmsg)
  if (stat /= 0) call fail(errmsg)
  call output%write_solution(0.0_r8, pressure, velocity, temporal_output, flow_active, &
      compliance=compliance, void_target_divergence=target_divergence)
  call output%close()

  ! Mainline collapse publishes only the signed target divergence.
  env%output_dir = 'mainline'
  call output%open(env, mesh, matl_model, temporal_output, stat, errmsg, collapse_output=.true.)
  if (stat /= 0) call fail(errmsg)
  do c = 1, mesh%ncell_onP
    if (flow_active(c)) target_divergence(c) = -abs(target_divergence(c))
  end do
  call output%write_solution(0.0_r8, pressure, velocity, temporal_output, flow_active, &
      void_target_divergence=target_divergence)
  call output%close()

  deallocate(mesh)
  call halt_parallel_communication

contains

  subroutine fail(message)
    character(*), intent(in) :: message
    if (env%rank == 0) write (*, '(2a)') 'FAIL: ', message
    call MPI_Abort(env%comm, 1)
  end subroutine fail

end program test_t2d_flow_vtkhdf_nan
