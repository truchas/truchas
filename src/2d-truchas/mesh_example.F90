!!
!! This base example provides the boilerplate startup code needed to initialize
!! MPI and other Truchas framework stuff, and the corresponding shutdown code.
!! It also shows how to instantiate a 2D unstructured mesh (which happens to be
!! a regular Cartesian mesh.)
!!

program mesh_example

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication
  use truchas_env, only: prefix
  use truchas_logging_services
  use unstr_2d_mesh_factory
  implicit none

  integer  :: nx(2)
  real(r8) :: xmin(2), xmax(2), eps, ptri
  type(unstr_2d_mesh), pointer :: mesh

  namelist /mesh2d/ xmin, xmax, nx, eps, ptri


  !! Initialize MPI and other base stuff that truchas depends on
  call init_parallel_communication
  prefix='run'  ! TLS will write to 'run.log'
  call TLS_initialize
  call TLS_set_verbosity(TLS_VERB_NOISY)

  eps = 0.0_r8  ! do not perturb the node positions by default
  ptri = 0.0_r8 ! quad mesh by default
  read(*,nml=mesh2d)
  mesh => new_unstr_2d_mesh(xmin, xmax, nx, eps, ptri)
  !mesh => new_unstr_2d_quad_mesh(xmin, xmax, nx, eps)
  !mesh => new_unstr_2d_tri_mesh(xmin, xmax, nx, eps)

  !! Cell volumes and face areas (okay, areas and lengths in 2D) are defined
  !! by default. But cell centroids and face centroids must be "requested".
  call mesh%init_cell_centroid
  call mesh%init_face_centroid

  !! Shutdown MPI
  call halt_parallel_communication

end program
