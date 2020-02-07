!!
!! This base example provides the boilerplate startup code needed to initialize
!! MPI and other Truchas framework stuff, and the corresponding shutdown code.
!! It also shows how to instantiate a 2D unstructured mesh (which happens to be
!! a regular Cartesian mesh.) It also illustrates how to generate a graphics
!! file that Paraview can read.
!!

program mesh_example

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use pgslib_module, only: PGSLib_CL_MAX_TOKEN_LENGTH, pgslib_finalize, pgslib_barrier
  use parallel_util_module, only: parallel_init
  use parallel_communication
  use truchas_env, only: prefix
  use truchas_logging_services
  use unstr_2d_mesh_factory
  use xdmf_file_type
  implicit none

  character(PGSLib_CL_MAX_TOKEN_LENGTH), pointer :: argv(:) => null()

  integer  :: nx(2)
  real(r8) :: xmin(2), xmax(2), eps, ptri
  type(unstr_2d_mesh), pointer :: mesh

  namelist /mesh2d/ xmin, xmax, nx, eps, ptri

  integer :: j
  real(r8), allocatable :: p(:), v(:,:)
  type(xdmf_file) :: outfile

  !! Initialize MPI and other base stuff that truchas depends on
  call parallel_init(argv)
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

  !! Define a cell-based field on the mesh with arbitrary value.
  allocate(p(mesh%ncell))
  do j = 1, mesh%ncell
    p(j) = dot_product([1.0_r8, -2.0_r8], mesh%cell_centroid(:,j))
  end do

  !! Define a node-based vector field on the mesh with arbitrary value.
  allocate(v(2,mesh%nnode))
  do j = 1, mesh%nnode
    v(:,j) = mesh%x(:,j)
  end do

  !! Create XDMF-format input files for Paraview: .xmf XML metadata file and
  !! .bin binary data file. Load the .xmf file into Paraview.
  call outfile%open('example')
  call outfile%write_mesh(mesh)

  !! Write P and V fields at a time snapshot
  call outfile%begin_variables(0.0_r8)
  call outfile%write_cell_var(p, 'P')
  call outfile%write_node_var(v, 'V')
  call outfile%end_variables

  !! Write P and V fields at a later time snapshot
  p = 2*p; v = 2*v
  call outfile%begin_variables(0.1_r8)
  call outfile%write_cell_var(p, 'P')
  call outfile%write_node_var(v, 'V')
  call outfile%end_variables

  !! Close the files (and add closing tags to the .xmf XML file)
  call outfile%close

  !! Shutdown MPI
  call pgslib_finalize

end program
