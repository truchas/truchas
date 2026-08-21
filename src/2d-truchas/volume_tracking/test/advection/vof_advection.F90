!!
!! This code attempts to initialize a VOF field and advect it with a constant
!! velocity using it's own time-step driver and VOF routines.
!! It also instantiates a 2D unstructured mesh (which happens to be a regular
!! Cartesian mesh).
!!
!! Aditya K. Pandare <apandare@lanl.gov>, January 2020
!! SPDX-License-Identifier: BSD-3-Clause
!!

program vof_advection

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mpi_f08, only: MPI_COMM_WORLD, MPI_Comm_rank, MPI_Comm_size
  use parallel_communication
  use simulation_environment_type
  use unstr_2d_mesh_factory
  use read_inputfile
  use gaussian_quadrature_vofinit
  use vof_2d_test_driver
  use vof_2d_vtkhdf_writer_type
  implicit none

  character(512) :: inputfile, arg
  integer  :: nx(2), tsmax, nmat, nvtrack
  real(r8) :: xmin(2), xmax(2), dxeps, ptri, dt, r
  type(unstr_2d_mesh), pointer :: mesh
  type(vof_2d_vtkhdf_writer) :: vtk_writer
  type(simulation_environment) :: env
  integer :: stat
  character(:), allocatable :: errmsg

  integer :: i, j, ngp
  real(r8) :: t_start, t_end, coord(2)
  real(r8), allocatable :: v(:,:), vof(:,:), int_normal(:,:,:)
  real(r8), allocatable :: vel_fn(:) ! fluxing velocity stored at faces
  real(r8), allocatable :: gp_coord(:,:), gp_weight(:)
  logical :: axisym

  procedure(constant_vel), pointer :: problem_vel => NULL()

  call cpu_time(t_start)

  !! Initialize MPI and other base stuff that truchas depends on
  call init_parallel_communication
  env%comm = MPI_COMM_WORLD
  call MPI_Comm_rank(env%comm, env%rank)
  call MPI_Comm_size(env%comm, env%nproc)
  call env%simlog%init(env%comm, 'mesh.log', stat, errmsg, terminal_output=.false.)
  if (stat /= 0) error stop 'initializing simulation log: ' // errmsg

  inputfile = 'input.txt'
  call get_command_argument(1, arg)
  if (len_trim(arg) > 0) inputfile = trim(arg)
  call readfile(inputfile, xmin, xmax, nx, dxeps, ptri, tsmax, dt, nmat, nvtrack)

  !! Create the mesh specified by the above input file
  mesh => new_unstr_2d_mesh(env, xmin, xmax, nx, dxeps, ptri)

  !! Cell volumes and face areas (okay, areas and lengths in 2D) are defined
  !! by default. But cell centroids and face centroids must be "requested".
  call mesh%init_cell_centroid
  call mesh%init_face_centroid

  !! Define a node-based vector field on the mesh with arbitrary value.
  allocate(v(2,mesh%nnode))
  do j = 1, mesh%nnode
    v(:,j) = mesh%x(:,j)
  end do

  !! Define a face-based normal-velocity field with arbitrary constant value.
  allocate(vel_fn(mesh%nface))
  problem_vel => constant_vel
  axisym = .false.

  !! Define a cell-based VOF field.
  allocate(vof(nmat,mesh%ncell))
  !! Initialize a "circular" VOF field
  ngp = 16
  allocate(gp_coord(2,ngp), gp_weight(ngp))

  do j = 1, mesh%ncell
    associate (cn => mesh%cnode(mesh%cstart(j):mesh%cstart(j+1)-1))

      call quadrature_qua4(ngp, gp_coord, gp_weight)
      vof(:,j) = 0.0_r8

      do i = 1, ngp
        call transform_qua4(mesh%x(:,cn), gp_coord(:,i), coord)
        r = norm2(coord(:)-0.3_r8)
        if (r<=0.125_r8) then
          vof(1,j) = vof(1,j) + gp_weight(i)*1.0_r8
        !else if (r<=0.15) then
        !  vof(2,j) = vof(2,j) + gp_weight(i)*1.0_r8
        end if
      end do !i

      vof(nmat,j) = 1.0_r8 - sum(vof(1:nmat-1,j))
    end associate
  end do

  !! Allocate a cell-based interface normal field.
  allocate(int_normal(2,nmat,mesh%ncell))
  int_normal = 0.0_r8

  call vtk_writer%open(env, mesh, nmat, stat, errmsg)
  if (stat /= 0) error stop 'opening VTKHDF output: ' // errmsg
  call vtk_writer%write_solution(0.0_r8, vof)

  !! call time-step driver
  call timestep_driver(tsmax, dt, mesh, vel_fn, nmat, nvtrack, problem_vel, vof, int_normal, axisym)

  call vtk_writer%write_solution(real(tsmax, r8)*dt, vof)
  call vtk_writer%close()

  !! Shutdown MPI
  call halt_parallel_communication

  call cpu_time(t_end)

  write(*,*) "Runtime: ", t_end-t_start

end program vof_advection
