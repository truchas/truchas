#include "f90_assert.fpp"

program example

  use kinds, only: r8
  use pgslib_module
  use parallel_util_module, only: parallel_init
  use parallel_communication
  use dist_mesh_type
  use mesh_broker
  use rad_driver_type
  use physical_constants, only: read_physical_constants
  use function_namelist, only: read_function_namelists
  use ER_input
  use tbrook_module, only: tbrook_initialize  ! for read_physical_constants
  use parallel_permutations
  use ER_driver_gmv
  implicit none
  
  logical :: exists
  type(rad_problem) :: prob
  type(dist_mesh), pointer :: mesh
  character(len=PGSLib_CL_MAX_TOKEN_LENGTH), pointer :: argv(:)
  real(r8), allocatable :: hc_temp(:), temp(:), qrad(:), res(:)
  real(r8) :: centroid(3), res_norm
  integer :: j, p
  character(len=31), allocatable :: encl_name(:)
  integer, parameter :: lun = 10
  
  !! Essential Truchas initialization.
  call parallel_init (argv)
  call init_parallel_communication ()
  call tbrook_initialize (istatus=p)  ! for read_physical_constants (huh?)
  
  if (is_IOP) open(lun,file='example.inp',position='rewind',action='read')
  
  !! Read essential Truchas namelists.
  call read_mesh_namelists (lun)
  call read_physical_constants (lun)
  call read_function_namelists (lun)
  
  !! Enclosure radiation namelists.
  call ERI_read_enclosure_radiation (lun)
  call ERI_read_enclosure_surface (lun)
  
  !! The mesh must be named 'main' for this example program.
  call enable_mesh ('main', exists)
  INSIST(exists)
  
  !! Read and initialize the mesh.
  call init_mesh_broker ()
  mesh => named_mesh_ptr('main')
  INSIST(associated(mesh%fnode))
  
  !! Define an interesting temperature field over the mesh.
  allocate(hc_temp(mesh%nface))
  do j = 1, mesh%nface
    centroid = sum(mesh%x(:,mesh%fnode(:,j)),dim=2) / size(mesh%fnode,dim=1)
    hc_temp(j) = 1.0_r8 + 2.0_r8 * sum(centroid**2)
  end do
  
  !! Get the list of enclosure names.
  allocate(encl_name(ERI_num_enclosures()))
  call ERI_get_names (encl_name)
  
  do p = 1, size(encl_name)
  
    if (is_IOP) write(*,'(/,a)') 'Solving radiosity system "' // trim(encl_name(p)) // '" ...'
  
    call prob%init (mesh, encl_name(p))
  
    allocate(temp(prob%nface_hc), qrad(prob%nface_hc), res(prob%nface_hc))
    
    call ERD_gmv_open (trim(encl_name(p)) // '.gmv')
    call ERD_gmv_write_enclosure (prob)
    call ERD_gmv_begin_variables ()

    !! The given enclosure surface temperature.
    temp = hc_temp(prob%faces)
    call ERD_gmv_write_var (prob, temp, 'temp')
    
    !! Compute the surface radiosity.
    qrad = 0.0_r8
    call prob%solve_radiosity (0.0_r8, temp, qrad)
    call ERD_gmv_write_var (prob, qrad, 'radiosity')
  
    !! Compute the residual of the radiosity system (should be small).
    call prob%residual  (0.0_r8, qrad, temp, res)
    res_norm = sqrt(global_sum(res**2))
    if (is_IOP) write(*,'(2x,a,es13.5)') 'residual norm =', res_norm
  
    !! Compute the heat flux through the surface.
    call prob%heat_flux (0.0_r8, qrad, res)
    call ERD_gmv_write_var (prob, res,  'heatflux')
    
    call ERD_gmv_end_variables ()
    call ERD_gmv_close ()
    
    if (is_IOP) write(*,'(2x,a)') 'graphics written to "' // trim(encl_name(p)) // '.gmv"'
    
    deallocate(temp, qrad, res)
    
    call ERD_problem_destroy (prob)
    
  end do
  
  call halt_parallel_communication ()

end program
