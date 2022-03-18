program test_pressure_poisson
#ifdef NAGFOR
  use,intrinsic :: f90_unix, only: exit
#endif
  use kinds
  use parallel_communication
  use truchas_env, only: prefix
  use truchas_logging_services
  use parameter_list_type
  use flow_operators
  use flow_props_type
  use flow_mesh_type
  use flow_bc_type
  use flow_type
  implicit none

  type(parameter_list), pointer :: p, pp
  type(flow_props) :: props
  type(flow_bc), pointer :: bc
  type(flow_mesh), pointer :: mesh
  type(flow) :: f
  integer :: i, in
  real(r8), allocatable :: flux_volumes(:,:)

  call init_parallel_communication
  prefix='run'  ! TLS will write to 'run.log'
  call TLS_initialize()
  call TLS_set_verbosity(TLS_VERB_NOISY)
  call init_mesh(mesh)
  call flow_operators_init(mesh)
  call init_flow_props(mesh, props)
  call init_flow_bcs(mesh, bc)

  allocate(p)
  pp => p%sublist("projection")
  pp => pp%sublist("solver")
  call pp%set('krylov-method', 'cg')
  call pp%set('rel-tol', 1.0e-8_r8)
  call pp%set('conv-rate-tol', 1.0e-8_r8)
  call pp%set('max-ds-iter', 50)
  call pp%set('max-amg-iter', 20)
  call pp%set('cg-use-two-norm', .true.)
  call pp%set('logging-level', 1)
  call pp%set('print-level', 0)

  call f%read_params(p)
  call f%init(mesh, bc, P_cc=1.0_r8)
  allocate(flux_volumes(1,size(mesh%mesh%cface)))
  flux_volumes = 0.0_r8
  call f%step(0.0_r8, 1.0_r8, props, flux_volumes)

  open(file='p_cc', newunit=in)
  associate (m => mesh%mesh, cc => mesh%cell_centroid, p => f%P_cc)
    do i = 1, m%ncell_onP
      write(in, '(4es15.5)') cc(:,i), p(i)
    end do
  end associate
  close(in)
  print *, "done writing file"
  call halt_parallel_communication()

contains

  subroutine init_mesh(mesh)

    use mesh_manager
    use unstr_mesh_type

    type(flow_mesh), pointer, intent(out) :: mesh
    character(:), allocatable :: indir
    type(parameter_list) :: pmesh
    type(parameter_list), pointer :: p
    type(unstr_mesh), pointer :: umesh
    integer :: i

    call process_command_line(indir)

    p => pmesh%sublist("mesh")
    call p%set('mesh-file', indir//"/square-32x32.g")
    call init_mesh_manager(pmesh)

    umesh => unstr_mesh_ptr("MESH")
    allocate(mesh)
    call mesh%init(umesh)
    print *, "INTIALIZED MESH"

    print *, "Information on cell 2"
    print *, "cell centroid: ", mesh%cell_centroid(:,2)
    print *, "faces of cell 2: ", umesh%cface(umesh%xcface(2):umesh%xcface(3)-1)
    print *, "cnhbr of cell 2: ", umesh%cnhbr(umesh%xcnhbr(2):umesh%xcnhbr(3)-1)
    print *, "face centroids of cell 2: "
    associate (fn => umesh%cface(umesh%xcface(2):umesh%xcface(3)-1))
      do i = 1, size(fn)
        print *, "centroid: ", fn(i), mesh%face_centroid(:,fn(i))
        if (btest(umesh%cfpar(2),pos=i)) then
          print *, fn(i), " IS INWARD"
          print *, "normal: ", fn(i), -umesh%normal(:,fn(i))
        else
          print *, fn(i), " IS OUTWARD"
          print *, "normal: ", fn(i), umesh%normal(:,fn(i))
        end if
      end do
    end associate
  end subroutine init_mesh

  subroutine init_flow_props(mesh, props)
    use scalar_func_containers
    use scalar_func_factories
    type(flow_mesh), pointer, intent(in) :: mesh
    type(flow_props), intent(inout) :: props
    !-
    type(scalar_func_box) :: density_delta(1), viscosity(1)
    real(r8) :: density(1)
    logical :: contains_void
    real(r8), allocatable :: vof(:,:), temperature(:)
    type(parameter_list) :: p

    density = [1.0_r8]
    contains_void = .false.
    call alloc_const_scalar_func(density_delta(1)%f, 0.0_r8)
    call alloc_const_scalar_func(viscosity(1)%f, 1.0_r8)

    call props%read_params(p)
    call props%init(mesh, density, density_delta, viscosity, contains_void)

    associate (nc => mesh%mesh%ncell)
      allocate(vof(1,nc))
      allocate(temperature(nc))
      vof = 1.0_r8
      temperature = 1.0_r8

      call props%update_cc(vof, temperature, initial=.true.)
      call props%update_fc(initial=.true.)
    end associate

    print *, "INTIALIZED FLOW PROPS"
  end subroutine init_flow_props


  subroutine init_flow_bcs(mesh, bc)
    type(flow_mesh), pointer, intent(in) :: mesh
    type(flow_bc), pointer, intent(out) :: bc
    !-
    type(parameter_list), pointer :: p, pp

    allocate(p)
!!$    pp => p%sublist("bc")
!!$    pp => pp%sublist('yz')
    pp => p%sublist('yz')
    call pp%set('condition', 'velocity-dirichlet')
    call pp%set('face-sets', [3,4,5])
    call pp%set('data', [0.0_r8, 0.0_r8, 0.0_r8])
!!$    pp => p%sublist("bc")
!!$    pp => pp%sublist('yz-p')
    pp => p%sublist('yz-p')
    call pp%set('condition', 'pressure-neumann')
    call pp%set('face-sets', [3,4,5])
    call pp%set('data', 0.0_r8)
!!$    pp => p%sublist("bc")
!!$    pp => pp%sublist('xm')
    pp => p%sublist('xm')
    call pp%set('condition', 'pressure-dirichlet')
    call pp%set('face-sets', [1])
    call pp%set('data', 0.0_r8)
!!$    pp => p%sublist("bc")
!!$    pp => pp%sublist('xp')
    pp => p%sublist('xp')
    call pp%set('condition', 'pressure-dirichlet')
    call pp%set('face-sets', [2])
    call pp%set('data', 1.0_r8)

    allocate(bc)
    call bc%init(mesh, p)

    print *, "INTIALIZED FLOW BCS"
  end subroutine init_flow_bcs


!!$  subroutine pressure_poisson_test(paramstr, dir)
!!$
!!$    use material_properties_type
!!$
!!$    character(*), intent(in) :: paramstr
!!$    integer, intent(in) :: dir
!!$
!!$    type(projection_solver) :: projection
!!$    type(parameter_list), pointer :: projection_params => null()
!!$    character(:), allocatable :: errmsg
!!$    type(amrex_multifab) :: velocity, pressure, fluxing_velocity(3), gradp_dyn_rho, fluid_rho, vof, &
!!$        temperature
!!$    type(amrex_string) :: label(2)
!!$    real(r8) :: linf, linf_global, xcen(3)
!!$    real(r8), pointer :: pressurep(:,:,:,:)
!!$    type(amrex_mfiter) :: mfi
!!$    type(amrex_box) :: bx
!!$    integer :: ix,iy,iz
!!$    type(material_properties) :: matl_prop
!!$
!!$    call parameter_list_from_json_string(paramstr, projection_params, errmsg)
!!$    INSIST(associated(projection_params))
!!$    call projection%init(mesh, projection_params)
!!$
!!$    call mesh%multifab_build(fluxing_velocity(1), nc=1, ng=1, nodal=xface_nodal)
!!$    call mesh%multifab_build(fluxing_velocity(2), nc=1, ng=1, nodal=yface_nodal)
!!$    call mesh%multifab_build(fluxing_velocity(3), nc=1, ng=1, nodal=zface_nodal)
!!$    call mesh%multifab_build(velocity, nc=3, ng=1)
!!$    call mesh%multifab_build(pressure, nc=1, ng=1)
!!$    call mesh%multifab_build(gradp_dyn_rho, nc=3, ng=1)
!!$    call mesh%multifab_build(fluid_rho, nc=1, ng=1)
!!$    call mesh%multifab_build(vof, nc=1, ng=1)
!!$    call mesh%multifab_build(temperature, nc=1, ng=1)
!!$
!!$    call fluxing_velocity(1)%setval(0.0_r8)
!!$    call fluxing_velocity(2)%setval(0.0_r8)
!!$    call fluxing_velocity(3)%setval(0.0_r8)
!!$    call velocity%setval(0.0_r8)
!!$    call pressure%setval(0.0_r8)
!!$    call gradp_dyn_rho%setval(0.0_r8)
!!$    call vof%setval(1.0_r8)
!!$    call fluid_rho%setval(1.0_r8)
!!$    call temperature%setval(0.0_r8)
!!$    matl_prop%density = 1
!!$    matl_prop%viscosity = 0
!!$
!!$    call projection%solve(1.0_r8, 0.0_r8, [0.0_r8, 0.0_r8, 0.0_r8], 1.0_r8, matl_prop, &
!!$        fluid_rho, vof, temperature, velocity, pressure, fluxing_velocity, gradp_dyn_rho, status)
!!$
!!$    ! check error
!!$    linf = 0
!!$    call mesh%mfiter_build(mfi)
!!$    do while (mfi%next())
!!$      bx = mfi%tilebox()
!!$      pressurep => pressure%dataptr(mfi)
!!$
!!$      do iz = bx%lo(3),bx%hi(3)
!!$        do iy = bx%lo(2),bx%hi(2)
!!$          do ix = bx%lo(1),bx%hi(1)
!!$            xcen = mesh%geom%get_physical_location([ix,iy,iz]) + mesh%geom%dx / 2
!!$            linf = max(linf, abs(pressurep(ix,iy,iz,1) - xcen(dir)))
!!$          end do
!!$        end do
!!$      end do
!!$    end do
!!$
!!$    call MPI_Reduce(linf, linf_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
!!$
!!$    if (this_rank==0) then
!!$      print *, 'linf: ', linf_global
!!$      if (linf_global > 1e-11_r8) status = 1
!!$    end if
!!$    call MPI_Bcast(status, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
!!$
!!$    ! call amrex_string_build(label(1), 'pressure')
!!$    ! call amrex_write_plotfile('plt0', 1, [pressure], label, [mesh%geom], 0.0_r8, [1], [1])
!!$
!!$  end subroutine pressure_poisson_test


    !! Handle a single optional argument which is the directory where the mesh
  !! files are located.  This is needed for cmake/ctest which puts and runs
  !! the executable from a different directory that this source code and mesh
  !! files.  Also handle '-h' or '--help' as an option to write out the usage.

  subroutine process_command_line (indir)

    use,intrinsic :: iso_fortran_env, only: output_unit

    character(:), allocatable, intent(out) :: indir

    character(:), allocatable :: prog
    character(256) :: arg
    integer :: n, stat

    call get_command_argument (0, arg)
    n = scan(arg, '/', back=.true.) ! remove the leading path component, if any
    prog = trim(arg(n+1:))

    stat = 0
    select case (command_argument_count())
    case (0)
      indir = '.'
    case (1)
      call get_command_argument (1, arg)
      select case (arg)
      case ('-h','--help')
        stat = 1
      case default
        indir = trim(arg)
      end select
    case default
      stat = 1
    end select

    if (stat /= 0) then
      if (is_IOP) then
        write(output_unit,'(a)') 'Usage: ' // prog // ' [INDIR]'
        write(output_unit,'(a)') 'INDIR is the directory containing the input files (default ".")'
      end if
      call exit (stat) ! do not want this to count as a successful test
    end if

  end subroutine process_command_line

end program test_pressure_poisson
