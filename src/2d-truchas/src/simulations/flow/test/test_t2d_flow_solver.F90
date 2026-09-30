program test_t2d_flow_solver
  use,intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mpi_f08, only: MPI_COMM_WORLD, MPI_Comm_rank, MPI_Comm_size
  use parallel_communication
  use fhypre, only: fhypre_initialize
  use truchas_env, only: prefix, overwrite_output
  use truchas_logging_services
  use parameter_list_type
  use parameter_list_json
  use simulation_environment_type
  use t2d_unstr_mesh_type
  use t2d_unstr_mesh_factory
  use material_database_type
  use material_model_type
  use material_factory, only: load_material_database
  use t2d_flow_model_type
  use flow_domain_types, only: regular_void_t
  use t2d_flow_solver_type
  implicit none

  integer :: status, stat
  character(:), allocatable :: errmsg
  type(simulation_environment) :: env

  call init_parallel_communication
  call fhypre_initialize
  prefix = 'run'
  overwrite_output = .true.
  call TLS_initialize
  call TLS_set_verbosity(TLS_VERB_NORMAL)
  env%comm = MPI_COMM_WORLD
  call MPI_Comm_rank(env%comm, env%rank)
  call MPI_Comm_size(env%comm, env%nproc)
  call env%simlog%init(env%comm, 'test_t2d_flow_solver.log', stat, errmsg, terminal_output=.false.)
  if (stat /= 0) call TLS_fatal('initializing simulation log: ' // errmsg)

  status = 0
  call test_step
  call test_collapse(0.0_r8)
  call test_collapse(0.5_r8)
  call test_collapse(0.0_r8, mainline=.true.)
  call test_mainline_history

  call halt_parallel_communication
  stop status

contains

  subroutine test_step
    type(t2d_unstr_mesh), pointer :: mesh
    type(t2d_flow_model), target :: model
    type(t2d_flow_solver), target :: solver
    type(material_database) :: database
    type(material_model) :: matl_model
    type(parameter_list), pointer :: matl_params, plist, momentum_params, projection_params, tracking_params
    type(parameter_list), target :: bc_params, solver_params
    real(r8), allocatable :: velocity(:,:), flux(:), vfrac(:,:), temperature(:)
    real(r8), pointer :: pressure(:), velocity_state(:,:), velocity_face(:)
    character(:), allocatable :: errmsg
    integer :: stat

    mesh => new_unstr_2d_mesh(env, [0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [8, 8], 0.0_r8, 0.0_r8)
    plist => bc_params%sublist('wall')
    call plist%set('type', 'no-slip')
    call plist%set('face-set-ids', [1,2,3,4])
    call model%init_core(env, mesh, bc_params, stat, errmsg)
    call require(stat == 0, 'Navier--Stokes model initialization failed')
    if (stat /= 0) return
    call parameter_list_from_json_string( &
        '{"liquid":{"properties":{"fluid":true,"density":1.0,"viscosity":0.1}}}', matl_params, errmsg)
    call require(associated(matl_params), 'parsing material database failed')
    if (.not.associated(matl_params)) return
    call load_material_database(database, matl_params, stat, errmsg)
    call require(stat == 0, 'loading material database failed')
    if (stat /= 0) return
    call matl_model%init(['liquid'], database, stat, errmsg)
    call require(stat == 0, 'material model initialization failed')
    if (stat /= 0) return
    momentum_params => solver_params%sublist('momentum-solver')
    call set_solver_params(momentum_params)
    projection_params => solver_params%sublist('projection-solver')
    call set_solver_params(projection_params)
    tracking_params => solver_params%sublist('volume-tracking')
    call tracking_params%set('algorithm', 'simple')
    call solver%init(env, model, matl_model, solver_params, stat, errmsg)
    call require(stat == 0, 'Navier--Stokes solver initialization failed')
    if (stat /= 0) return
    allocate(velocity(2,mesh%ncell_onP), flux(mesh%ncell_onP), vfrac(1,mesh%ncell), temperature(mesh%ncell_onP))
    vfrac = 1.0_r8
    temperature = 0.0_r8
    velocity = spread([1.0_r8, -0.5_r8], dim=2, ncopies=mesh%ncell_onP)
    call solver%set_initial_material_state(vfrac, temperature)
    call solver%set_initial_state(env, 0.0_r8, 0.01_r8, velocity, stat)
    call require(stat == 0, 'Navier--Stokes initial-condition solve did not converge')
    if (stat /= 0) return
    call solver%step(env, 0.0_r8, 0.01_r8, stat, errmsg)
    call require(stat == 0, 'Navier--Stokes solver step did not converge')
    if (stat /= 0) return
    call solver%get_cell_flow_soln(pressure, velocity_state)
    call solver%get_face_velocity(velocity_face)
    call model%operators%divergence(velocity_face, flux)
    call require(maxval(abs(flux)) < 1.0e-8_r8, 'Navier--Stokes step did not make face velocity solenoidal')
    call require(solver%courant_time_step() > 0.0_r8, 'Navier--Stokes Courant time step is not positive')
  end subroutine


  subroutine test_collapse(sigma, mainline)
    real(r8), intent(in) :: sigma
    logical, optional, intent(in) :: mainline
    type(t2d_unstr_mesh), pointer :: mesh
    type(t2d_flow_model), target :: model
    type(t2d_flow_solver), target :: solver
    type(material_database) :: database
    type(material_model) :: matl_model
    type(parameter_list), pointer :: matl_params, plist, bc_params, tracking_params
    type(parameter_list), target :: model_params, solver_params
    real(r8), allocatable :: velocity(:,:), vfrac(:,:), temperature(:), bias(:), flux(:), compliance(:)
    real(r8), pointer :: velocity_face(:), pressure(:), velocity_state(:,:), target_divergence(:)
    real(r8) :: before, after, initial_volume, inflow, dt, time
    character(:), allocatable :: errmsg
    integer :: stat, c, f, step
    logical :: mainline_, diagnostic_ok

    mainline_ = .false.
    if (present(mainline)) mainline_ = mainline
    mesh => new_unstr_2d_mesh(env, [0.0_r8,0.0_r8], [1.0_r8,1.0_r8], [8,8], 0.0_r8, 0.0_r8)
    call model_params%set('inviscid', .true.)
    plist => model_params%sublist('void-collapse')
    if (mainline_) then
      call plist%set('model', 'mainline')
    else
      call plist%set('pressure-time-scale', 1.0_r8)
    end if
    if (sigma > 0.0_r8) call plist%set('capillary-coefficient', sigma)
    bc_params => model_params%sublist('bc')
    plist => bc_params%sublist('wall')
    call plist%set('type', 'free-slip')
    call plist%set('face-set-ids', [2,3,4])
    plist => bc_params%sublist('feed')
    call plist%set('type', 'pressure')
    call plist%set('face-set-ids', [1])
    call plist%set('pressure', 1.0_r8)
    call model%init(env, mesh, model_params, stat, errmsg)
    call require(stat == 0, 'collapse model initialization failed')
    if (stat /= 0) return
    call parameter_list_from_json_string( &
        '{"liquid":{"properties":{"fluid":true,"density":1.0}}}', matl_params, errmsg)
    call load_material_database(database, matl_params, stat, errmsg)
    call require(stat == 0, 'collapse material database failed')
    if (stat /= 0) return
    call matl_model%init(['liquid','VOID  '], database, stat, errmsg)
    call require(stat == 0, 'collapse material model failed')
    if (stat /= 0) return
    plist => solver_params%sublist('projection-solver')
    call set_solver_params(plist)
    tracking_params => solver_params%sublist('volume-tracking')
    call tracking_params%set('algorithm', 'simple')
    call solver%init(env, model, matl_model, solver_params, stat, errmsg)
    call require(stat /= 0, 'collapse accepted simple tracking')
    call tracking_params%set('algorithm', 'geometric')
    call tracking_params%set('cutoff', 1.0e-10_r8)
    call solver%init(env, model, matl_model, solver_params, stat, errmsg)
    call require(stat == 0, 'collapse solver initialization failed')
    if (stat /= 0) return
    allocate(velocity(2,mesh%ncell_onP), vfrac(2,mesh%ncell), temperature(mesh%ncell_onP))
    vfrac(1,:) = 1.0_r8
    do c = 1, mesh%ncell
      if (mesh%cell_centroid(1,c) > 0.875_r8) vfrac(1,c) = 0.5_r8
    end do
    vfrac(2,:) = 1.0_r8-vfrac(1,:)
    temperature = 0.0_r8
    velocity = 0.0_r8
    dt = 0.005_r8
    call solver%set_initial_material_state(vfrac, temperature)
    bias = model%collapse_pressure_bias()
    do c = 1, mesh%ncell_onP
      if (mesh%cell_centroid(1,c) > 0.875_r8) then
        bias(c) = bias(c)-sigma*0.125_r8/sqrt(mesh%volume(c))
      end if
    end do
    call require(maxval(abs(bias)) < 1.0e-12_r8, 'wrong capillary pressure bias')
    initial_volume = global_sum(sum(mesh%volume(:mesh%ncell_onP)*model%matl_props%vof_novoid(:mesh%ncell_onP)))
    call solver%set_initial_state(env, 0.0_r8, dt, velocity, stat)
    call require(stat == 0, 'collapse initial solve failed')
    if (stat /= 0) return
    call solver%get_cell_flow_soln(pressure, velocity_state, target_divergence)
    compliance = model%collapse_compliance()
    diagnostic_ok = .true.
    do c = 1, mesh%ncell_onP
      if (compliance(c) > 0.0_r8 .or. (mainline_ .and. model%matl_props%cell_t(c) == regular_void_t)) then
        diagnostic_ok = diagnostic_ok .and. ieee_is_finite(target_divergence(c))
      else
        diagnostic_ok = diagnostic_ok .and. ieee_is_nan(target_divergence(c))
      end if
    end do
    call require(diagnostic_ok, 'incorrect initial compliance target mask')
    allocate(flux(mesh%ncell_onP))
    do step = 1, 20
      before = global_sum(sum(mesh%volume(:mesh%ncell_onP)*model%matl_props%vof_novoid(:mesh%ncell_onP)))
      call solver%get_face_velocity(velocity_face)
      inflow = 0.0_r8
      do f = 1, mesh%nface_onP
        if (mesh%fcell(2,f) == 0) inflow = inflow-dt*mesh%area(f)*velocity_face(f)
      end do
      inflow = global_sum(inflow)
      time = (step-1)*dt
      call solver%step(env, time, time+dt, stat, errmsg)
      call require(stat == 0, 'coupled collapse step failed')
      if (stat /= 0) return
      after = global_sum(sum(mesh%volume(:mesh%ncell_onP)*model%matl_props%vof_novoid(:mesh%ncell_onP)))
      ! Initial velocity is zero; the first projection supplies the next transport velocity.
      if (step > 1 .and. .not.mainline_) &
          call require(after > before, 'stationary pocket stopped collapsing under pressure loading')
      call require(abs(after-before-inflow) < 1.0e-8_r8, 'collapse liquid gain differs from boundary inflow')
      call require(all(model%matl_props%vof_novoid >= 0.0_r8 .and. &
          model%matl_props%vof_novoid <= 1.0_r8), 'collapse produced unbounded liquid fractions')
      call solver%get_cell_flow_soln(pressure, velocity_state, target_divergence)
      call solver%get_face_velocity(velocity_face)
      call model%operators%divergence(velocity_face, flux)
      compliance = model%collapse_compliance()
      diagnostic_ok = .true.
      do c = 1, mesh%ncell_onP
        if (compliance(c) > 0.0_r8 .or. (mainline_ .and. model%matl_props%cell_t(c) == regular_void_t)) then
          if (ieee_is_finite(target_divergence(c))) then
            diagnostic_ok = diagnostic_ok .and. abs(target_divergence(c)-flux(c)/mesh%volume(c)) < 1.0e-8_r8
          else
            diagnostic_ok = .false.
          end if
        else
          diagnostic_ok = diagnostic_ok .and. ieee_is_nan(target_divergence(c))
        end if
      end do
      call require(diagnostic_ok, 'incorrect accepted compliance target divergence or mask')
    end do
    if (mainline_) then
      ! A pressure-loaded but unchanged pocket supplies no history-based source.
      call require(abs(after-initial_volume) < 1.0e-8_r8, 'mainline collapsed a stationary unchanged pocket')
    else
      call require(after > initial_volume+1.0e-4_r8, 'collapse did not replace measurable void with liquid')
    end if
  end subroutine


  subroutine test_mainline_history
    type(t2d_unstr_mesh), pointer :: mesh
    type(t2d_flow_model) :: model
    type(parameter_list), target :: params
    type(parameter_list), pointer :: plist, bc_params
    real(r8), allocatable :: vfrac(:,:), temperature(:), expected(:)
    character(:), allocatable :: errmsg
    integer :: stat, c

    mesh => new_unstr_2d_mesh(env, [0.0_r8,0.0_r8], [1.0_r8,1.0_r8], [8,8], 0.0_r8, 0.0_r8)
    call params%set('inviscid', .true.)
    bc_params => params%sublist('bc')
    plist => bc_params%sublist('walls')
    call plist%set('type', 'free-slip')
    call plist%set('face-set-ids', [1,2,3,4])
    plist => params%sublist('void-collapse')
    call plist%set('model', 'unknown')
    call model%init(env, mesh, params, stat, errmsg)
    call require(stat /= 0, 'unknown collapse model accepted')
    call plist%set('model', 'mainline')
    call model%init(env, mesh, params, stat, errmsg)
    call require(stat == 0, 'mainline collapse without compliance parameters failed')
    call require(model%collapse_enabled(), 'mainline collapse not enabled')
    if (stat /= 0) return
    call model%matl_props%init(mesh, [1.0_r8], .true., stat, errmsg, nfluid=2)
    call require(stat == 0, 'mainline material properties initialization failed')
    if (stat /= 0) return
    allocate(vfrac(2,mesh%ncell), temperature(mesh%ncell_onP), expected(mesh%ncell_onP))
    temperature = 0.0_r8
    vfrac(1,:) = 0.6_r8
    vfrac(2,:) = 0.4_r8
    call model%matl_props%set_initial_state(vfrac, temperature)
    call require(all(model%collapse_fraction() == 0.0_r8), 'nonzero initial collapse source')
    vfrac(1,:) = 0.65_r8
    vfrac(2,:) = 0.35_r8
    call model%matl_props%set_volume_fractions(vfrac)
    call require(maxval(abs(model%collapse_fraction()-0.035_r8)) < 1.0e-14_r8, &
        'wrong mainline source for decreasing VOID; default relaxation cap not applied')
    call require(maxval(abs(model%matl_props%void_old-0.4_r8)) < 1.0e-14_r8, &
        'trial property update overwrote committed VOID')
    vfrac(1,:) = 0.59_r8
    vfrac(2,:) = 0.41_r8
    call model%matl_props%set_volume_fractions(vfrac)
    call require(maxval(abs(model%collapse_fraction()-0.01_r8)) < 1.0e-14_r8, &
        'wrong mainline source for increasing VOID')
    ! Rejection restores current fractions without advancing their history.
    vfrac(1,:) = 0.6_r8
    vfrac(2,:) = 0.4_r8
    call model%matl_props%set_volume_fractions(vfrac)
    call require(maxval(abs(model%collapse_fraction())) < 1.0e-14_r8, 'nonzero source after restoring fractions')
    vfrac(1,:) = 0.59_r8
    vfrac(2,:) = 0.41_r8
    call model%matl_props%set_volume_fractions(vfrac)
    call model%accept_material_state()
    call require(all(model%collapse_fraction() == 0.0_r8), 'accept failed to advance VOID history')

    ! Mainline includes solid/liquid/VOID cells, unlike the compliance prototype.
    vfrac(1,:) = 0.3_r8
    vfrac(2,:) = 0.2_r8
    call model%matl_props%set_initial_state(vfrac, temperature)
    vfrac(1,:) = 0.29_r8
    vfrac(2,:) = 0.21_r8
    call model%matl_props%set_volume_fractions(vfrac)
    call require(maxval(abs(model%collapse_fraction()-0.01_r8)) < 1.0e-14_r8, &
        'mainline excluded solid/liquid/VOID cells')
    call require(all(model%collapse_compliance() == 0.0_r8), 'mainline enabled pressure compliance')
    call require(all(model%collapse_pressure_bias() == 0.0_r8), 'mainline enabled pressure bias')
    do c = 1, mesh%ncell
      if (mesh%cell_centroid(1,c) < 0.125_r8) then
        vfrac(:,c) = [1.0_r8, 0.0_r8]
      else if (mesh%cell_centroid(1,c) < 0.25_r8) then
        vfrac(:,c) = 0.0_r8
      else if (mesh%cell_centroid(1,c) > 0.875_r8) then
        vfrac(:,c) = [0.0_r8, 0.5_r8] ! solid/VOID is classified as VOID
      end if
    end do
    call model%matl_props%set_volume_fractions(vfrac)
    expected = 0.01_r8
    do c = 1, mesh%ncell_onP
      if (mesh%cell_centroid(1,c) < 0.25_r8 .or. mesh%cell_centroid(1,c) > 0.75_r8) expected(c) = 0.0_r8
    end do
    call require(maxval(abs(model%collapse_fraction()-expected)) < 1.0e-14_r8, &
        'wrong mainline eligibility: pure fluid, solid, VOID, or VOID neighbor')

    call plist%set('relaxation', -0.1_r8)
    call model%init(env, mesh, params, stat, errmsg)
    call require(stat /= 0, 'negative mainline relaxation accepted')
    call plist%set('relaxation', 1.1_r8)
    call model%init(env, mesh, params, stat, errmsg)
    call require(stat /= 0, 'mainline relaxation above one accepted')
    call plist%set('relaxation', 0.0_r8)
    call model%init(env, mesh, params, stat, errmsg)
    call require(stat == 0, 'zero mainline relaxation rejected')
    call model%matl_props%init(mesh, [1.0_r8], .true., stat, errmsg, nfluid=2)
    call model%matl_props%set_volume_fractions(vfrac)
    call require(all(model%collapse_fraction() == 0.0_r8), 'zero relaxation produced a source')
  end subroutine


  subroutine set_solver_params(params)
    type(parameter_list), intent(inout) :: params

    call params%set('rel-tol', 1.0e-10_r8)
    call params%set('max-ds-iter', 100)
    call params%set('max-amg-iter', 100)
  end subroutine


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (global_any(.not.condition)) then
      if (is_IOP) print '("ERROR: ",a)', message
      status = 1
    end if
  end subroutine

end program test_t2d_flow_solver
