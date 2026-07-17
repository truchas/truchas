program main

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication
  use functions
  use face_boundary_values
  use cell_based_function
  use properties
  use material_properties
  use mesh_manager
  use spec_diff
  use bdf2_dae
  implicit none
  
  integer :: n, stat, status
  logical :: exists
  real(r8) :: tout, t, h
  real(r8), allocatable :: u(:), udot(:)
  type(sf_s) :: g
  type(prop), pointer :: diff_array(:) => null()
  type(state) :: s
  character(len=64) :: errmsg
  real :: cpu0, cpu1
  character(len=PGSLib_CL_MAX_TOKEN_LENGTH), pointer :: argv(:) => null()
  
  call init_parallel_communication(argv)
  
  if (is_IOP) call cpu_time (cpu0)

  open(unit=10,file='test1.inp',position='rewind',action='read',status='old')
  call read_mesh_namelists (10)
  call enable_mesh ('testmesh', exists)
  call init_mesh_manager ()
  
  if (is_IOP) then
    call cpu_time (cpu1)
    write(*,*) 'Mesh initialization time =', cpu1-cpu0
    cpu0 = cpu1
  end if
  
  mesh => named_mesh_ptr('testmesh')
  
  !! Nonlinear diffusion problem -Div((0.01+u)*Grad(u)) = 1.
  !! Homogeneous Neumann conditions on x,y,z=0 (symmetry planes).
  !! 0-Dirichlet conditions on x,y,z=1.

  !! Setup the Dirichlet BC.
  call face_bv_prep (dir, mesh)
  call face_bv_add  (dir, 0.0_r8, (/2,4,6/), stat, errmsg)
  call check_stat (stat, 'Error setting up Dirichlet BC: ' // trim(errmsg))
  call face_bv_done (dir)

  !! Setup the homogeneous Neuman BC.
  call face_bv_prep (neu, mesh)
  call face_bv_add  (neu, 0.0_r8, (/1,3,5/), stat, errmsg)
  call check_stat (stat, 'Error setting up Neumann BC: ' // trim(errmsg))
  call face_bv_done (neu)

  !! Setup the diffusion coefficient.
  !call cbf_prep (d, mesh)
  !call create_sf_v_user (g, index=1)
  !call cbf_add_u_fun (d, g, (/1/), stat, errmsg)
  !call destroy (g)
  !call cbf_done (d, stat, errmsg)
  !call check_stat (stat, 'Error setting up diffusion coefficient: ' // trim(errmsg))
  
  !! Single material; VFRAC is unity.
  allocate(vfrac(1,mesh%ncell))
  vfrac = 1.0_r8
  
  !! Species diffusivity for material 1.
  allocate(diff_array(1))
  !call create_sf_s_user (g, index=1)
  call create_sf_s_poly (g, (/0.01_r8, 1.0_r8/), (/0,1/))
  call create_prop_c_fun (diff_array(1), g)
  call destroy (g)
  
  !! (Total) species diffusivity object.
  call create_mat_prop (diffusivity, diff_array, vfrac)
  
  !! Setup the source function.
  call cbf_prep (q, mesh)
  call cbf_add_const (q, 1.0_r8, (/1/), stat, errmsg)
  call cbf_done (q, stat, errmsg)
  call check_stat (stat, 'Error setting up source function: ' // trim(errmsg))
  
  call spec_diff_init ()

  !! Initial value.
  n = mesh%ncell + mesh%nface_on_PE
  allocate(u(n), udot(n))
  u = 0.0_r8
  !udot = 0.0_r8 ! Completely wrong -- start with an extra small time step!
  t = 0.0d0
  call eval_udot (t, u, udot)
  call user (t, u)
  
  if (is_IOP) then
    call cpu_time (cpu1)
    write(*,*) 'Problem initialization time =', cpu1-cpu0
    cpu0 = cpu1
  end if

 !!
 !! INITIALIZE THE INTEGRATION PROCEDURE
 !!
  
  rtol = 0.0d0
  atol = 1.0d-5
  
  nsweep = 4
  omega = 1.4_r8
  
  call create_state (s, size(u), mitr=5, mvec=4, ntol=0.1d0)
  call set_initial_state (s, t, u, udot)
  
  open(unit=1,file='bdf2.out')
  call verbose_stepping (s, 1)
  
  tout = 1.0d0
  h = 1.0d-5
  
  call bdf2_step_driver (s, pcfun, updpc, enorm, h, status, tout=tout)
  select case (status)
  case (SOLVED_TO_TOUT)
    if (is_IOP) then
      call cpu_time (cpu1)
      write(*,*) 'Problem integration time =', cpu1-cpu0
      cpu0 = cpu1
    end if
    call user (tout, interpolated_solution (s, tout))
  case (SOLVED_TO_NSTEP)
    print *, 'Failed to reach final time within specified time steps.'
  case default
    print *, 'BDF2_STEP_DRIVER returned an error condition: status=', status
    print *, t
  end select
  
  if (is_IOP) call write_bdf2_stepping_statistics (s, 6)
  
  call halt_parallel_communication ()

contains

  subroutine check_stat (stat, errmsg)
    integer, intent(in) :: stat
    character(len=*), intent(in) :: errmsg
    if (global_any(stat /= 0)) then
      write(0,*) errmsg
      stop
    end if
  end subroutine check_stat

  subroutine user (t, u)
    use gmv_writer
    real(r8), intent(in) :: t, u(:)
    integer, save :: n=0
    character(len=12) :: fname
    write(fname,'(a,i4.4)') 'out.gmv.', n
    call gmv_open (fname)
    call gmv_write_dist_mesh (mesh)
    call gmv_begin_variables (t)
    call gmv_write_dist_cell_var (mesh, u(:mesh%ncell), name='u')
    call gmv_end_variables ()
    call gmv_close ()
    n = n + 1
  end subroutine user

end program main

function user_sf_s (index, s) result (f)
  use,intrinsic :: iso_fortran_env, only: r8 => real64
  integer,  intent(in) :: index
  real(r8), intent(in) :: s
  real(r8) :: f
  select case (index)
  case (1)
    if (s < 0.0_r8) print *, 'FOOBAR', s
    f = 0.01_r8 + max(s, 0.0_r8)
  case default
    print *, 'bad user_sf_s index:', index
    stop
  end select
end function user_sf_s
