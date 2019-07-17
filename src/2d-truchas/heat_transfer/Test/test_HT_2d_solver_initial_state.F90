!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_HT_2d_solver_initial_state

#ifdef NAGFOR
  use,intrinsic :: f90_unix, only: exit
#endif
  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication
  use truchas_env, only: prefix, overwrite_output
  use truchas_logging_services
  use parameter_list_type
  use parameter_list_json
  use unstr_2d_mesh_factory
  use matl_mesh_func_type
  use material_database_type
  use source_mesh_function
  use scalar_func_factories
  use mfd_2d_disc_type
  use HT_2d_model_type
  use HT_2d_solver_type
  use bitfield_type
  use test_ht_2d_common
  implicit none

  type(unstr_2d_mesh), pointer :: mesh
  type(mfd_2d_disc), target :: mfd_disc
  type(material_database), target :: matl_db
  type(matl_mesh_func), target :: mmf
  real(r8) :: xmin(2), xmax(2), tol, eps
  integer  :: nx(2)
  integer :: status = 0

  !! Initialize MPI and other base stuff that Truchas depends on
  call init_parallel_communication
  prefix='run'  ! TLS will write to 'run.log'
  overwrite_output = .true.
  call TLS_initialize
  call TLS_set_verbosity(TLS_VERB_NORMAL)

  TOL = 1E-9_r8
  eps = 0.0_r8  ! mesh distortion
  xmin = [0.0_r8, 0.0_r8]
  xmax = [1.0_r8, 1.0_r8]
  nx  = [64, 64]

  !! Create the mesh specified by the above input file
  mesh => new_unstr_2d_mesh(xmin, xmax, nx, eps)

  !! Initialize state needed by all tests
  call mfd_disc%init(mesh)
  call init_materials(mesh, matl_db, mmf)

  !! Run test problems
  call test_linear_dir(mfd_disc, mmf, tol)
  call test_linear_flux(mfd_disc, mmf, tol)

  !! Wrap up
  call halt_parallel_communication
  call exit (status)

contains

  subroutine error_exit(errmsg)
    character(len=*), intent(in) :: errmsg
    status = 1
    if (is_IOP) print '("FATAL: ",a)', errmsg
    call halt_parallel_communication
    call exit(status)
  end subroutine error_exit

  !! Initialize parameter list for all tests
  subroutine init_params(params)

    type(parameter_list), intent(out) :: params
    type(parameter_list), pointer :: sublist

    call params%set_name('solver')

    sublist => params%sublist('preconditioner')
    call sublist%set('method','BoomerAMG')
    sublist => sublist%sublist('params')
    call sublist%set('num-cycles', 1)

    sublist => params%sublist('error-norm')
    call sublist%set('temp-rel-tol', 1.0d-4)
    call sublist%set('enth-rel-tol', 1.0d-4)

    sublist => params%sublist('integrator')
    call sublist%set('nlk-max-iter', 5)
    call sublist%set('nlk-tol', 0.01_r8)

    call params%set('hmin', tiny(1.0_r8))
    call params%set('max_step_tries', 10)

  end subroutine init_params

  !! Tests the HT_2d_solver on a linear problem with Dirichlet boundary conditions
  subroutine test_linear_dir(disc, mmf, tol)

    type(mfd_2d_disc), target, intent(in) :: disc
    type(matl_mesh_func), target, intent(in) :: mmf
    real(r8), intent(in) :: tol

    type(HT_2d_solver) :: HT_solver
    type(HT_2d_model), target :: HT_model
    type(parameter_list) :: solver_params
    type(parameter_list), pointer :: model_params
    class(scalar_func), allocatable :: f
    integer :: exps(3,2) = reshape([0,1,0,0,0,1],[3,2])  ! exponents of u(x,t)
    real(r8) :: lcoef(2) = [1.0_r8, 2.0_r8]  ! coefficients of u(x,t)
    real(r8), allocatable, target :: u(:), udot(:)
    real(r8), allocatable :: state(:,:), Hcell(:), Tcell(:), Tface(:)
    real(r8), pointer :: view(:)
    character(:), allocatable :: errmsg, string
    integer :: n, stat, max_itr
    real(r8) :: t, dt, rel_tol

    if (is_IOP) print '(/,"Testing linear problem with Dirichlet BCs")'

    t = 0.0_r8
    dt = 1E-3_r8

    !! 2D HT model parameters
    string = &
    '{"bc": { &
        "all-sides": { &
          "type": "temperature", &
          "face-set-ids": [1,2,3,4], &
          "temp": { &
            "type": "polynomial", &
            "poly-coef": [1.0, 2.0], &
            "poly-powers": [[0,1,0],[0,0,1]]}}}}'
    call parameter_list_from_json_string(string, model_params, errmsg)
    if (.not. associated(model_params)) call error_exit(errmsg)

    !! Initialize 2D HT model
    call HT_model%init(disc, mmf, model_params, stat, errmsg)
    if (stat /= 0) call error_exit(errmsg)

    !! Define state variables
    n = HT_model%num_dof()
    allocate(u(n), udot(n))

    !! Define the linear function ax+by
    call alloc_mpoly_scalar_func(f, lcoef, exps)

    !! Compute cell-based and face-based fields on the mesh.
    allocate(Tcell(disc%mesh%ncell_onP), Tface(disc%mesh%nface_onP))
    call average_integral(disc, f, Tcell, Tface)

    !! Initialize 2D HT solver
    call init_params(solver_params)
    call HT_solver%init(HT_model, solver_params)

    !! Expected cell enthalpy
    allocate(Hcell(disc%mesh%ncell))
    call HT_model%set_cell_temp(Tcell, u)
    call HT_model%new_state_array(u, state)
    call HT_model%H_of_T%compute_value(state, Hcell)
    deallocate(state)

    !! Compute consistent u and udot
    max_itr = 100
    rel_tol = tol * 1E-5_r8 !TODO: good practice?
    call HT_solver%compute_initial_state(t, dt, Tcell, u, udot, rel_tol, max_itr, stat, errmsg)
    if (stat/=0) call error_exit(errmsg)

    !! u must match expected values
    call HT_model%get_cell_heat_view(u, view)
    if (global_any(abs(Hcell(:size(view))-view) > tol)) then
      if (is_IOP) print '("ERROR: cell enthalpy exceeds expected value; tol=",es9.2)', tol
      status = 1
    end if
    call HT_model%get_cell_temp_view(u, view)
    if (global_any(abs(Tcell-view) /= 0.0_r8)) then
      if (is_IOP) print '("ERROR: cell temp exceeds expected value; tol=",es9.2)', tol
      status = 1
    end if
    call HT_model%get_face_temp_view(u, view)
    if (global_any(abs(Tface-view) > tol)) then
      if (is_IOP) print '("ERROR: face temp exceeds expected value; tol=",es9.2)', tol
      status = 1
    end if

    !! udot must be 0
    call HT_model%get_cell_heat_view(udot, view)
    if (global_any(abs(view) > tol)) then
      if (is_IOP) print '("ERROR: Hdot is nonzero; tol=",es9.2)', tol
      status = 1
    end if
    call HT_model%get_cell_temp_view(udot, view)
    if (global_any(abs(view) > tol)) then
      if (is_IOP) print '("ERROR: cell temp derivative is nonzero; tol=",es9.2)', tol
      status = 1
    end if
    call HT_model%get_face_temp_view(udot, view)
    if (global_any(abs(view) > tol)) then
      if (is_IOP) print '("ERROR: face temp derivative is nonzero; tol=",es9.2)', tol
      status = 1
    end if

  end subroutine test_linear_dir


  !! Tests the HT_2d_solver on a linear problem with Neumann boundary conditions
  subroutine test_linear_flux(disc, mmf, tol)

    type(mfd_2d_disc), target, intent(in) :: disc
    type(matl_mesh_func), target, intent(in) :: mmf
    real(r8), intent(in) :: tol

    type(HT_2d_solver) :: HT_solver
    type(HT_2d_model), target :: HT_model
    type(parameter_list) :: solver_params
    type(parameter_list), pointer :: model_params
    class(scalar_func), allocatable :: f
    integer :: exps(3,2) = reshape([0,1,0,0,0,1],[3,2])  ! exponents of u(x,t)
    real(r8) :: lcoef(2) = [1.0_r8, 2.0_r8]  ! coefficients of u(x,t)
    real(r8) :: dcoef = 1.0_r8  ! diffusion coefficient
    real(r8), allocatable, target :: u(:), udot(:)
    real(r8), allocatable :: state(:,:), Hcell(:), Tcell(:), Tface(:)
    real(r8), pointer :: view(:)
    character(:), allocatable :: errmsg, string
    real(r8) :: t, dt, rel_tol
    integer :: n, stat, max_itr

    if (is_IOP) print '(/,"Testing linear problem with Neumann BCs")'

    t = 0.0_r8
    dt = 1E-3_r8

    !! 2D HT model parameters
    string = &
    '{"bc": { &
        "left": { &
          "type": "flux", &
          "face-set-ids": [1], &
          "flux": 1.0 }, &
        "right": { &
          "type": "flux", &
          "face-set-ids": [2], &
          "flux": -1.0 }, &
        "bottom": { &
          "type": "flux", &
          "face-set-ids": [3], &
          "flux": 2.0 }, &
        "top": { &
          "type": "flux", &
          "face-set-ids": [4], &
          "flux": -2.0 }}, &
      "source": { &
        "all-cells": { &
          "cell-set-ids": [1], &
          "source": 1.0 }}}'
    call parameter_list_from_json_string(string, model_params, errmsg)
    if (.not. associated(model_params)) call error_exit(errmsg)

    !! Initialize 2D HT model
    call HT_model%init(disc, mmf, model_params, stat, errmsg)
    if (stat /= 0) call error_exit(errmsg)

    !! Define state variables
    n = HT_model%num_dof()
    allocate(u(n), udot(n))

    !! Define the linear function ax+by
    call alloc_mpoly_scalar_func(f, lcoef, exps)

    !! Compute cell-based and face-based fields on the mesh.
    allocate(Tcell(disc%mesh%ncell_onP), Tface(disc%mesh%nface_onP))
    call average_integral(disc, f, Tcell, Tface)

    !! Initialize 2D HT solver
    call init_params(solver_params)
    call HT_solver%init(HT_model, solver_params)

    !! Expected cell enthalpy
    allocate(Hcell(disc%mesh%ncell))
    call HT_model%set_cell_temp(Tcell, u)
    call HT_model%new_state_array(u, state)
    call HT_model%H_of_T%compute_value(state, Hcell)
    deallocate(state)

    !! Compute consistent u and udot
    max_itr = 100
    rel_tol = tol * 1E-5_r8 !TODO: good practice?
    call HT_solver%compute_initial_state(t, dt, Tcell, u, udot, rel_tol, max_itr, stat, errmsg)
    if (stat/=0) call error_exit(errmsg)

    !! u must match expected values
    call HT_model%get_cell_heat_view(u, view)
    if (global_any(abs(Hcell(:size(view))-view) > tol)) then
      if (is_IOP) print '("ERROR: cell enthalpy exceeds expected value; tol=",es9.2)', tol
      status = 1
    end if
    call HT_model%get_cell_temp_view(u, view)
    if (global_any(abs(Tcell-view) > tol)) then
      if (is_IOP) print '("ERROR: cell temp exceeds expected value; tol=",es9.2)', tol
      status = 1
    end if
    call HT_model%get_face_temp_view(u, view)
    if (global_any(abs(Tface-view) > tol)) then
      if (is_IOP) print '("ERROR: face temp exceeds expected value; tol=",es9.2)', tol
      status = 1
    end if

    !! udot must match expected values
    call HT_model%get_cell_heat_view(udot, view)
    if (global_any(abs(view-dcoef) > tol)) then
      if (is_IOP) print '("ERROR: Hdot exceeds expected value; tol=",es9.2)', tol
      status = 1
    end if
    call HT_model%get_cell_temp_view(udot, view)
    if (global_any(abs(view) > tol)) then
      if (is_IOP) print '("ERROR: cell temp derivative is nonzero; tol=",es9.2)', tol
      status = 1
    end if
    call HT_model%get_face_temp_view(udot, view)
    if (global_any(abs(view) > tol)) then
      if (is_IOP) print '("ERROR: face temp derivative is nonzero; tol=",es9.2)', tol
      status = 1
    end if

  end subroutine test_linear_flux

end program test_HT_2d_solver_initial_state
