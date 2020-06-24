!!
!! This code attempts to initialize a VOF field and advect it with a constant
!! velocity using it's own time-step driver and VOF routines.
!! It also instantiates a 2D unstructured mesh (which happens to be a regular
!! Cartesian mesh) and generates a graphics file that Paraview can read.
!!

program vof_vortex

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use pgslib_module, only: PGSLib_CL_MAX_TOKEN_LENGTH, pgslib_finalize, pgslib_barrier
  use parallel_util_module, only: parallel_init
  use parallel_communication
  use truchas_env, only: prefix
  use truchas_logging_services
  use unstr_2d_mesh_factory
  use xdmf_file_type
  use read_inputfile
  use gaussian_quadrature_vofinit
  use vof_2d_test_driver
  implicit none

  character(PGSLib_CL_MAX_TOKEN_LENGTH), pointer :: argv(:) => null()

  character(len=100) :: inputfile
  integer  :: nx(2), tsmax, nmat, nvtrack
  real(r8) :: xmin(2), xmax(2), dxeps, ptri, dt, r
  type(unstr_2d_mesh), pointer :: mesh

  integer :: i, j, ngp, test_run, nelem, gncell
  integer, allocatable :: global_xcell(:)
  real(r8) :: t_start, t_end, coord(2), vof_err
  real(r8), allocatable :: vof(:,:), int_normal(:,:,:), vof_std(:), global_vof(:), myproc(:)
  real(r8), allocatable :: vel_fn(:) ! fluxing velocity stored at faces
  real(r8), allocatable :: gp_coord(:,:), gp_weight(:)
  logical :: test_failure, axisym
  type(xdmf_file) :: outfile

  procedure(constant_vel), pointer :: problem_vel => NULL()

  call cpu_time(t_start)

  !! Initialize MPI and other base stuff that truchas depends on
  call parallel_init(argv)
  call init_parallel_communication
  prefix='run'  ! TLS will write to 'run.log'
  call TLS_initialize
  call TLS_set_verbosity(TLS_VERB_NOISY)

  !! Read input file "input_vortex.txt"
  inputfile = 'input_vortex.txt'
  call readfile(inputfile, xmin, xmax, nx, dxeps, ptri, tsmax, dt, nmat, nvtrack, test_run)

  !! Create the mesh specified by the above input file
  mesh => new_unstr_2d_mesh(xmin, xmax, nx, dxeps, ptri)

  !! Cell volumes and face areas (okay, areas and lengths in 2D) are defined
  !! by default. But cell centroids and face centroids must be "requested".
  call mesh%init_cell_centroid
  call mesh%init_face_centroid

  !! Define a face-based normal-velocity field with arbitrary constant value.
  allocate(vel_fn(mesh%nface))
  problem_vel => vortex_vel
  axisym = .false.

  !! Define a cell-based VOF field.
  allocate(vof(nmat,mesh%ncell), myproc(mesh%ncell))
  !! Initialize a "circular" VOF field
  ngp = 16
  allocate(gp_coord(2,ngp), gp_weight(ngp))

  do j = 1, mesh%ncell
    associate (cn => mesh%cnode(mesh%cstart(j):mesh%cstart(j+1)-1))

      call quadrature_qua4(ngp, gp_coord, gp_weight)
      vof(:,j) = 0.0_r8

      myproc(j) = this_pe

      do i = 1, ngp
        call transform_qua4(mesh%x(:,cn), gp_coord(:,i), coord)
        r = norm2(coord(1:2)-[0.5_r8, 0.75_r8])
        if (r<=0.15_r8) then
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

  !! Create XDMF-format input files for Paraview: .xmf XML metadata file and
  !! .bin binary data file. Load the .xmf file into Paraview.
  call outfile%open('vof')
  call outfile%write_mesh(mesh)

  !! Write P and V fields at a time snapshot
  call outfile%begin_variables(0.0_r8)
  call outfile%write_cell_var(vof(1,:), 'VOF1')
  call outfile%write_cell_var(vof(2,:), 'VOF2')
  call outfile%write_cell_var(myproc(:), 'mype')
  call outfile%write_cell_var(int_normal(1,1,:), 'x-normal')
  call outfile%write_cell_var(int_normal(2,1,:), 'y-normal')
  call outfile%end_variables

  !! call time-step driver
  call timestep_driver(tsmax, dt, mesh, vel_fn, nmat, nvtrack, problem_vel, vof, outfile, &
    int_normal, axisym, myproc)

  !! Close the files (and add closing tags to the .xmf XML file)
  call outfile%close

  if (test_run == 0) then
    open(3, file='circlevort_vof.txt')

    write(3,'(I8)') mesh%ncell
    do j = 1, mesh%ncell
      write(3,'(I8, E20.10)') j, vof(1,j)
    end do

    close(3)
  end if

  !! Collect local VOF arrays into a global VOF array
  gncell = global_sum(mesh%ncell_onP)

  allocate(global_vof(merge(gncell,0,is_IOP)))
  allocate(global_xcell(merge(gncell,0,is_IOP)))

  call collate(global_vof, vof(1,:mesh%ncell_onP))
  call collate(global_xcell, mesh%xcell(:mesh%ncell_onP))

  if (is_iop) global_vof(global_xcell) = global_vof

  ! Testing
  if (test_run == 1 .and. is_iop) then
    test_failure = .false.
    write(*,*) 'Comparing output to circlevort_vof.txt'
    open(3, file='circlevort_vof.txt', action='read', status='old')

    read(3,*) nelem
    if (nelem /= gncell) then
      call TLS_fatal('Number of mesh cells in test standard and current mesh do not match')
    end if

    allocate(vof_std(nelem))

    do j = 1, nelem
      read(3,*) i, vof_std(j)
    end do

    vof_err = 0.0_r8
    do j = 1, nelem
      vof_err = max(vof_err, abs(vof_std(j)-global_vof(j)))
    end do
    test_failure = vof_err > 1e-07_r8

    close(3)
  end if

  !! Shutdown MPI
  call pgslib_finalize

  call cpu_time(t_end)

  if (test_run == 1 .and. is_iop) then
    if (test_failure) then
      write(*,*) "FAIL: VOF linf error=", vof_err, " (tol=1.00e-7)"
      stop 1
    else
      write(*,*) "PASS: VOF linf error=", vof_err, " (tol=1.00e-7)"
    end if
  end if

  write(*,*) "Runtime: ", t_end-t_start

end program vof_vortex
