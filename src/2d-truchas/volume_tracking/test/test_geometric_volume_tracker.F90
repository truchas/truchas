!!
!! Direct unit tests for the two-dimensional geometric volume tracker.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

program test_geometric_volume_tracker

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mpi_f08
  use parallel_communication
  use simulation_environment_type
  use unstr_2d_mesh_type
  use unstr_2d_mesh_factory, only: new_unstr_2d_quad_mesh
  use geom_axisymmetric, only: mesh_axisymmetry_mod
  use geometric_volume_tracker_type
  implicit none

  type(simulation_environment) :: env
  integer :: stat
  character(:), allocatable :: errmsg

  call MPI_Init
  call init_parallel_communication
  env%comm = MPI_COMM_WORLD
  call MPI_Comm_rank(env%comm, env%rank)
  call MPI_Comm_size(env%comm, env%nproc)
  call env%simlog%init(env%comm, 'test_geometric_volume_tracker.log', stat, errmsg, &
      terminal_output=.false.)
  if (stat /= 0) error stop 'initializing simulation log: ' // errmsg

  call test_stationary(env)
  call test_planar_transport(env)
  call test_immobile_solid(env)
  call test_inflow_material(env)
  call test_proportional_inflow(env)
  call test_three_material_transport(env)
  call test_material_depletion(env)
  call test_axisymmetric_transport(env)
  call test_axisymmetric_radial_transport(env)

  call env%simlog%close
  call MPI_Finalize

contains

  subroutine test_stationary(env)

    type(simulation_environment), intent(inout) :: env
    type(unstr_2d_mesh), pointer :: mesh
    type(geometric_volume_tracker) :: tracker
    real(r8), allocatable :: vel(:), vof_n(:,:), vof(:,:), flux_vol(:,:), int_normal(:,:,:)

    mesh => new_unstr_2d_quad_mesh(env, [0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [1, 1])
    call mesh%init_face_centroid
    call tracker%init(env, mesh, 2, 2, 2, .false., [1,2])
    allocate(vel(size(mesh%cface)), vof_n(2,mesh%ncell), vof(2,mesh%ncell), &
        flux_vol(2,size(mesh%cface)), int_normal(2,2,mesh%ncell))

    vel = 0.0_r8
    vof_n(:,1) = [0.3_r8, 0.7_r8]
    call tracker%flux_volumes(env, vel, vof_n, vof, flux_vol, int_normal, 2, 0, 0.1_r8)
    call require_close_2d(vof, vof_n, 'stationary volume fractions', 1.0e-13_r8)
    call require(maxval(abs(flux_vol)) <= 1.0e-13_r8, 'stationary flux volumes')

    deallocate(mesh)

  end subroutine test_stationary


  subroutine test_planar_transport(env)

    type(simulation_environment), intent(inout) :: env
    type(unstr_2d_mesh), pointer :: mesh
    type(geometric_volume_tracker) :: tracker
    real(r8), allocatable :: vel(:), vof_n(:,:), vof(:,:), flux_vol(:,:), int_normal(:,:,:)
    real(r8) :: expected(2,2), q
    integer :: f, interface_face

    mesh => new_unstr_2d_quad_mesh(env, [0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [2, 1])
    call mesh%init_face_centroid
    call tracker%init(env, mesh, 2, 2, 2, .false., [1,2])
    allocate(vel(size(mesh%cface)), vof_n(2,mesh%ncell), vof(2,mesh%ncell), &
        flux_vol(2,size(mesh%cface)), int_normal(2,2,mesh%ncell))

    vof_n = 0.0_r8
    vof_n(:,1) = [1.0_r8, 0.0_r8]
    vof_n(:,2) = [0.0_r8, 1.0_r8]
    call uniform_face_velocity(mesh, 0.1_r8, vel)
    call tracker%flux_volumes(env, vel, vof_n, vof, flux_vol, int_normal, 2, 0, 0.1_r8)
    interface_face = 0
    do f = 1, mesh%nface
      if (mesh%fcell(2,f) /= 0) interface_face = f
    end do
    ASSERT(interface_face > 0)
    q = 0.1_r8 * 0.1_r8 * mesh%area(interface_face) / mesh%volume(1)
    expected(:,1) = [1.0_r8, 0.0_r8]
    expected(:,2) = [q, 1.0_r8-q]
    call require_close_2d(vof, expected, 'planar transport', 1.0e-12_r8)
    call require_close_1d(sum(vof, dim=1), 1.0_r8, 'planar fraction sum', 1.0e-12_r8)
    call require(all(vof >= 0.0_r8 .and. vof <= 1.0_r8), 'planar boundedness')

    deallocate(mesh)

  end subroutine test_planar_transport


  subroutine test_immobile_solid(env)

    type(simulation_environment), intent(inout) :: env
    type(unstr_2d_mesh), pointer :: mesh
    type(geometric_volume_tracker) :: tracker
    real(r8), allocatable :: vel(:), vof_n(:,:), vof(:,:), flux_vol(:,:), int_normal(:,:,:)
    integer :: f, j1, j2

    mesh => new_unstr_2d_quad_mesh(env, [0.0_r8, 0.0_r8], [2.0_r8, 1.0_r8], [2, 1])
    call mesh%init_face_centroid
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

    vof_n(:,1) = [0.4_r8, 0.6_r8]
    vof_n(:,2) = [0.6_r8, 0.4_r8]
    vel = 0.0_r8
    vel(j1) = 0.1_r8
    vel(j2) = -0.1_r8
    call tracker%flux_volumes(env, vel, vof_n, vof, flux_vol, int_normal, 1, 0, 0.1_r8)

    call require_close_2d(vof(2:2,:), vof_n(2:2,:), 'immobile solid', 1.0e-12_r8)
    call require(maxval(abs(flux_vol(2,:))) <= 1.0e-12_r8, 'immobile solid flux')

    deallocate(mesh)

  end subroutine test_immobile_solid


  subroutine test_inflow_material(env)

    type(simulation_environment), intent(inout) :: env
    type(unstr_2d_mesh), pointer :: mesh
    type(geometric_volume_tracker) :: tracker
    real(r8), allocatable :: vel(:), vof_n(:,:), vof(:,:), flux_vol(:,:), int_normal(:,:,:)
    real(r8) :: expected(2,1), q
    integer :: f, left_face

    mesh => new_unstr_2d_quad_mesh(env, [0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [1, 1])
    call mesh%init_face_centroid
    call tracker%init(env, mesh, 2, 2, 2, .false., [1,2])
    allocate(vel(size(mesh%cface)), vof_n(2,mesh%ncell), vof(2,mesh%ncell), &
        flux_vol(2,size(mesh%cface)), int_normal(2,2,mesh%ncell))

    left_face = 0
    do f = 1, mesh%nface
      if (mesh%fcell(2,f) == 0 .and. mesh%face_centroid(1,f) == 0.0_r8) left_face = f
    end do
    ASSERT(left_face > 0)
    call tracker%set_inflow_material(1, [left_face])
    call uniform_face_velocity(mesh, 0.1_r8, vel)
    vof_n(:,1) = [0.0_r8, 1.0_r8]
    call tracker%flux_volumes(env, vel, vof_n, vof, flux_vol, int_normal, 2, 0, 0.1_r8)

    q = 0.1_r8 * 0.1_r8 * mesh%area(left_face) / mesh%volume(1)
    expected(:,1) = [q, 1.0_r8-q]
    call require_close_2d(vof, expected, 'inflow material', 1.0e-12_r8)

    deallocate(mesh)

  end subroutine test_inflow_material


  subroutine test_proportional_inflow(env)

    type(simulation_environment), intent(inout) :: env
    type(unstr_2d_mesh), pointer :: mesh
    type(geometric_volume_tracker) :: tracker
    real(r8), allocatable :: vel(:), vof_n(:,:), vof(:,:), flux_vol(:,:), int_normal(:,:,:)
    real(r8) :: expected_flux(2)
    integer :: f, left_face, left_local

    mesh => new_unstr_2d_quad_mesh(env, [0.0_r8, 0.0_r8], [2.0_r8, 1.0_r8], [2, 1])
    call mesh%init_face_centroid
    call tracker%init(env, mesh, 2, 2, 2, .false., [1,2])
    allocate(vel(size(mesh%cface)), vof_n(2,mesh%ncell), vof(2,mesh%ncell), &
        flux_vol(2,size(mesh%cface)), int_normal(2,2,mesh%ncell))

    vof_n(:,1) = [0.3_r8, 0.7_r8]
    vof_n(:,2) = vof_n(:,1)
    call uniform_face_velocity(mesh, 0.1_r8, vel)
    call tracker%flux_volumes(env, vel, vof_n, vof, flux_vol, int_normal, 2, 0, 0.1_r8)
    left_face = 0
    do f = 1, mesh%nface
      if (mesh%fcell(2,f) == 0 .and. mesh%face_centroid(1,f) == 0.0_r8) left_face = f
    end do
    ASSERT(left_face > 0)
    left_local = mesh%cstart(1)-1 + &
        findloc(mesh%cface(mesh%cstart(1):mesh%cstart(2)-1), left_face, dim=1)
    expected_flux = -vof_n(:,1) * 0.1_r8 * 0.1_r8 * mesh%area(left_face)
    call require(maxval(abs(flux_vol(:,left_local)-expected_flux)) <= 1.0e-12_r8, &
        'proportional inflow flux')
    call require_close_1d(sum(vof, dim=1), 1.0_r8, 'proportional inflow fraction sum', 1.0e-12_r8)
    call require(all(vof >= 0.0_r8 .and. vof <= 1.0_r8), 'proportional inflow boundedness')
    call require(maxval(abs(flux_vol)) > 0.0_r8, 'proportional inflow flux')

    deallocate(mesh)

  end subroutine test_proportional_inflow


  subroutine test_three_material_transport(env)

    type(simulation_environment), intent(inout) :: env
    type(unstr_2d_mesh), pointer :: mesh
    type(geometric_volume_tracker) :: tracker
    real(r8), allocatable :: vel(:), vof_n(:,:), vof(:,:), flux_vol(:,:), int_normal(:,:,:)

    mesh => new_unstr_2d_quad_mesh(env, [0.0_r8, 0.0_r8], [2.0_r8, 1.0_r8], [2, 1])
    call mesh%init_face_centroid
    call tracker%init(env, mesh, 3, 3, 3, .false., [1,2,3])
    allocate(vel(size(mesh%cface)), vof_n(3,mesh%ncell), vof(3,mesh%ncell), &
        flux_vol(3,size(mesh%cface)), int_normal(2,3,mesh%ncell))

    vof_n(:,1) = [0.2_r8, 0.3_r8, 0.5_r8]
    vof_n(:,2) = vof_n(:,1)
    call uniform_face_velocity(mesh, 0.1_r8, vel)
    call tracker%flux_volumes(env, vel, vof_n, vof, flux_vol, int_normal, 3, 0, 0.1_r8)
    call require(maxval(abs(vof(:,2)-vof_n(:,2))) <= 1.0e-12_r8, &
        'three-material undisturbed cell')
    call require_close_1d(sum(vof, dim=1), 1.0_r8, 'three-material fraction sum', 1.0e-12_r8)
    call require(all(vof >= 0.0_r8 .and. vof <= 1.0_r8), 'three-material boundedness')
    call require(vof(1,1) < vof(2,1) .and. vof(2,1) < vof(3,1), 'three-material ordering')

    deallocate(mesh)

  end subroutine test_three_material_transport


  subroutine test_material_depletion(env)

    type(simulation_environment), intent(inout) :: env
    type(unstr_2d_mesh), pointer :: mesh
    type(geometric_volume_tracker) :: tracker
    real(r8), allocatable :: vel(:), vof_n(:,:), vof(:,:), flux_vol(:,:), int_normal(:,:,:)
    real(r8) :: expected(2,2)
    integer :: f, left_face

    mesh => new_unstr_2d_quad_mesh(env, [0.0_r8, 0.0_r8], [2.0_r8, 1.0_r8], [2, 1])
    call mesh%init_face_centroid
    call tracker%init(env, mesh, 2, 2, 2, .false., [1,2])
    allocate(vel(size(mesh%cface)), vof_n(2,mesh%ncell), vof(2,mesh%ncell), &
        flux_vol(2,size(mesh%cface)), int_normal(2,2,mesh%ncell))

    left_face = 0
    do f = 1, mesh%nface
      if (mesh%fcell(2,f) == 0 .and. mesh%face_centroid(1,f) == 0.0_r8) left_face = f
    end do
    ASSERT(left_face > 0)
    call tracker%set_inflow_material(2, [left_face])
    call uniform_face_velocity(mesh, 1.0_r8, vel)
    vof_n(:,1) = [0.0_r8, 1.0_r8]
    vof_n(:,2) = [0.01_r8, 0.99_r8]
    call tracker%flux_volumes(env, vel, vof_n, vof, flux_vol, int_normal, 2, 0, 0.8_r8)

    expected(:,1) = [0.0_r8, 1.0_r8]
    expected(:,2) = [0.0_r8, 1.0_r8]
    call require_close_2d(vof, expected, 'material depletion', 1.0e-12_r8)
    call require_close_1d(sum(vof, dim=1), 1.0_r8, 'depletion fraction sum', 1.0e-12_r8)
    call require(all(vof >= 0.0_r8 .and. vof <= 1.0_r8), 'depletion boundedness')
    call require(maxval(abs(sum(flux_vol, dim=2) - [0.01_r8, -0.01_r8])) <= 1.0e-12_r8, &
        'depletion material balance')

    deallocate(mesh)

  end subroutine test_material_depletion


  subroutine test_axisymmetric_transport(env)

    type(simulation_environment), intent(inout) :: env
    type(unstr_2d_mesh), pointer :: mesh
    type(geometric_volume_tracker) :: tracker
    real(r8), allocatable :: vel(:), vof_n(:,:), vof(:,:), flux_vol(:,:), int_normal(:,:,:)
    real(r8) :: expected(2,2), q
    integer :: f, interface_face

    mesh => new_unstr_2d_quad_mesh(env, [1.0_r8, 0.0_r8], [2.0_r8, 2.0_r8], [1, 2])
    call mesh%init_face_centroid
    call mesh_axisymmetry_mod(mesh)
    call tracker%init(env, mesh, 2, 2, 2, .true., [1,2])
    allocate(vel(size(mesh%cface)), vof_n(2,mesh%ncell), vof(2,mesh%ncell), &
        flux_vol(2,size(mesh%cface)), int_normal(2,2,mesh%ncell))

    vof_n = 0.0_r8
    vof_n(:,1) = [1.0_r8, 0.0_r8]
    vof_n(:,2) = [0.0_r8, 1.0_r8]
    call uniform_face_velocity(mesh, 0.1_r8, vel, component=2)
    call tracker%flux_volumes(env, vel, vof_n, vof, flux_vol, int_normal, 2, 0, 0.1_r8)
    interface_face = 0
    do f = 1, mesh%nface
      if (mesh%fcell(2,f) /= 0) interface_face = f
    end do
    ASSERT(interface_face > 0)
    q = 0.1_r8 * 0.1_r8 * mesh%area(interface_face) / mesh%volume(1)
    expected(:,1) = [1.0_r8, 0.0_r8]
    expected(:,2) = [q, 1.0_r8-q]
    call require_close_2d(vof, expected, 'axisymmetric transport', 1.0e-12_r8)

    deallocate(mesh)

  end subroutine test_axisymmetric_transport


  subroutine test_axisymmetric_radial_transport(env)

    type(simulation_environment), intent(inout) :: env
    type(unstr_2d_mesh), pointer :: mesh
    type(geometric_volume_tracker) :: tracker
    real(r8), allocatable :: vel(:), vof_n(:,:), vof(:,:), flux_vol(:,:), int_normal(:,:,:)
    real(r8) :: q
    integer :: f, interface_face

    mesh => new_unstr_2d_quad_mesh(env, [1.0_r8, 0.0_r8], [3.0_r8, 1.0_r8], [2, 1])
    call mesh%init_face_centroid
    call mesh_axisymmetry_mod(mesh)
    call tracker%init(env, mesh, 2, 2, 2, .true., [1,2])
    allocate(vel(size(mesh%cface)), vof_n(2,mesh%ncell), vof(2,mesh%ncell), &
        flux_vol(2,size(mesh%cface)), int_normal(2,2,mesh%ncell))

    vof_n(:,1) = [1.0_r8, 0.0_r8]
    vof_n(:,2) = [0.0_r8, 1.0_r8]
    call uniform_face_velocity(mesh, 0.1_r8, vel, component=1)
    call tracker%flux_volumes(env, vel, vof_n, vof, flux_vol, int_normal, 2, 0, 0.1_r8)
    interface_face = 0
    do f = 1, mesh%nface
      if (mesh%fcell(2,f) /= 0) interface_face = f
    end do
    ASSERT(interface_face > 0)
    q = 0.1_r8 * 0.1_r8 * mesh%area(interface_face) / mesh%volume(2)
    call require(maxval(abs(vof(:,1)-[1.0_r8, 0.0_r8])) <= 1.0e-12_r8, &
        'axisymmetric radial inner cell')
    call require(abs(vof(1,2)-q) <= 0.25_r8*q, 'axisymmetric radial transported fraction')
    call require_close_1d(sum(vof, dim=1), 1.0_r8, 'axisymmetric radial fraction sum', 1.0e-12_r8)
    call require(all(vof >= 0.0_r8 .and. vof <= 1.0_r8), 'axisymmetric radial boundedness')

    deallocate(mesh)

  end subroutine test_axisymmetric_radial_transport


  subroutine uniform_face_velocity(mesh, speed, vel, component)

    type(unstr_2d_mesh), intent(in) :: mesh
    real(r8), intent(in) :: speed
    real(r8), intent(out) :: vel(:)
    integer, intent(in), optional :: component

    real(r8) :: face_velocity(mesh%nface)
    integer :: c, j, f, f0, f1, d

    d = 1
    if (present(component)) d = component
    face_velocity = speed * mesh%unit_normal(d,:)
    do c = 1, mesh%ncell
      f0 = mesh%cstart(c)
      f1 = mesh%cstart(c+1)-1
      do j = f0, f1
        f = mesh%cface(j)
        if (btest(mesh%cfpar(c), 1+j-f0)) then
          vel(j) = -face_velocity(f)
        else
          vel(j) = face_velocity(f)
        end if
      end do
    end do

  end subroutine uniform_face_velocity


  subroutine require(condition, message)

    logical, intent(in) :: condition
    character(*), intent(in) :: message
    if (.not.condition) error stop 'FAIL: ' // message

  end subroutine require


  subroutine require_close_2d(actual, expected, message, tolerance)

    real(r8), intent(in) :: actual(:,:), expected(:,:)
    character(*), intent(in) :: message
    real(r8), intent(in) :: tolerance
    ASSERT(all(shape(actual) == shape(expected)))
    call require(maxval(abs(actual-expected)) <= tolerance, message)

  end subroutine require_close_2d


  subroutine require_close_1d(actual, expected, message, tolerance)

    real(r8), intent(in) :: actual(:), expected
    character(*), intent(in) :: message
    real(r8), intent(in) :: tolerance
    call require(maxval(abs(actual-expected)) <= tolerance, message)

  end subroutine require_close_1d

end program test_geometric_volume_tracker
