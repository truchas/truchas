program test_flow_2d_ic_solver

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication
  use fhypre, only: fhypre_initialize
  use truchas_env, only: prefix, overwrite_output
  use truchas_logging_services
  use parameter_list_type
  use unstr_2d_mesh_type
  use unstr_2d_mesh_factory
  use flow_2d_model_type
  use flow_2d_state_type
  use flow_2d_solver_type
  implicit none

  integer :: status

  call init_parallel_communication
  call fhypre_initialize
  prefix = 'run'
  overwrite_output = .true.
  call TLS_initialize
  call TLS_set_verbosity(TLS_VERB_NORMAL)

  status = 0
  call test_uniform_velocity
  call test_incompatible_velocity

  call halt_parallel_communication
  stop status

contains

  subroutine test_uniform_velocity
    type(unstr_2d_mesh), pointer :: mesh
    type(flow_2d_model), target :: model
    type(flow_2d_state), target :: state
    type(flow_2d_solver) :: solver
    type(parameter_list), target :: bc_params, momentum_params, projection_params
    type(parameter_list), pointer :: plist
    real(r8), allocatable :: velocity(:,:), flux(:)
    character(:), allocatable :: errmsg
    integer :: stat

    mesh => new_unstr_2d_mesh([0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [8, 8], 0.0_r8, 0.0_r8)
    plist => bc_params%sublist('flow')
    call plist%set('type', 'velocity')
    call plist%set('face-set-ids', [1,2,3,4])
    call plist%set('velocity', [1.0_r8, 0.0_r8])
    call model%init(mesh, bc_params, 1.0_r8, 1.0_r8, stat, errmsg)
    call require(stat == 0, 'flow model initialization failed')
    if (stat /= 0) return
    call state%init(mesh)
    call set_solver_params(momentum_params)
    call set_solver_params(projection_params)
    call solver%init(model, state, momentum_params, projection_params)
    allocate(velocity(2,mesh%ncell), flux(mesh%ncell_onP))
    velocity = spread([1.0_r8, 0.0_r8], dim=2, ncopies=mesh%ncell)
    call solver%set_initial_state(0.0_r8, 1.0_r8, velocity, stat)
    call require(stat == 0, 'initial-condition solve did not converge')
    if (stat /= 0) return
    call model%operators%divergence(state%vel_fn, flux)
    call require(maxval(abs(flux)) < 1.0e-8_r8, 'initial face velocity is not solenoidal')
    call require(maxval(abs(state%vel_cc(:,1:mesh%ncell_onP) - velocity(:,1:mesh%ncell_onP))) < 1.0e-8_r8, &
        'compatible uniform initial velocity was changed')
  end subroutine


  subroutine test_incompatible_velocity
    type(unstr_2d_mesh), pointer :: mesh
    type(flow_2d_model), target :: model
    type(flow_2d_state), target :: state
    type(flow_2d_solver) :: solver
    type(parameter_list), target :: bc_params, momentum_params, projection_params
    type(parameter_list), pointer :: plist
    real(r8), allocatable :: velocity(:,:), flux(:)
    real(r8) :: boundary_velocity
    character(:), allocatable :: errmsg
    integer :: stat, f

    mesh => new_unstr_2d_mesh([0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [8, 8], 0.0_r8, 0.0_r8)
    plist => bc_params%sublist('wall')
    call plist%set('type', 'no-slip')
    call plist%set('face-set-ids', [1,2,3,4])
    call model%init(mesh, bc_params, 1.0_r8, 1.0_r8, stat, errmsg)
    call require(stat == 0, 'incompatible-velocity model initialization failed')
    if (stat /= 0) return
    call state%init(mesh)
    call set_solver_params(momentum_params)
    call set_solver_params(projection_params)
    call solver%init(model, state, momentum_params, projection_params)
    allocate(velocity(2,mesh%ncell), flux(mesh%ncell_onP))
    velocity = spread([1.0_r8, -0.5_r8], dim=2, ncopies=mesh%ncell)
    call solver%set_initial_state(0.0_r8, 1.0_r8, velocity, stat)
    call require(stat == 0, 'incompatible initial-condition solve did not converge')
    if (stat /= 0) return
    call model%operators%divergence(state%vel_fn, flux)
    call require(maxval(abs(flux)) < 1.0e-8_r8, &
        'repaired initial face velocity is not solenoidal')
    boundary_velocity = 0.0_r8
    do f = 1, mesh%nface_onP
      if (mesh%fcell(2,f) == 0) boundary_velocity = max(boundary_velocity, abs(state%vel_fn(f)))
    end do
    call require(boundary_velocity < 1.0e-8_r8, &
        'repaired initial velocity does not satisfy no-slip boundaries')
    call require(maxval(abs(state%vel_cc(:,1:mesh%ncell_onP) - velocity(:,1:mesh%ncell_onP))) > 1.0e-6_r8, &
        'incompatible initial velocity was not repaired')
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

end program test_flow_2d_ic_solver
