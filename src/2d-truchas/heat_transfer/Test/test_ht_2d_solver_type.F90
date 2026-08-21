!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_HT_2d_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication
  use fhypre, only: fhypre_initialize
  use parameter_list_type
  use parameter_list_json
  use unstr_2d_mesh_factory
  use material_model_type
  !use source_mesh_function
  use scalar_func_factories
  use mfd_2d_disc_type
  use ht_2d_model_type
  use ht_2d_solver_type
  use bitfield_type
  use test_ht_2d_common
  implicit none

  type(unstr_2d_mesh), pointer :: mesh
  type(mfd_2d_disc), target :: mfd_disc
  type(material_model), target :: matl_model
  real(r8) :: xmin(2), xmax(2), tol, eps
  integer  :: nx(2)
  integer :: status = 0
  real(r8), parameter :: PI = 3.1415926535897932_r8


  !! Initialize MPI and other base stuff that Truchas depends on
  call init_parallel_communication
  call fhypre_initialize
  call init_test_environment('test_ht_2d_solver_type.log')

  TOL = 2E-4_r8
  eps = 0.0_r8  ! mesh distortion
  xmin = [0.0_r8, 0.0_r8]
  xmax = [1.0_r8, 1.0_r8]
  nx  = [64, 64]

  !! Create the mesh specified by the above input file
  mesh => new_unstr_2d_mesh(test_env, xmin, xmax, nx, eps)

  !! Initialize state needed by all tests
  call mfd_disc%init(mesh)
  call init_materials(mesh, matl_model)

  !! Run test problems
  call test1(mfd_disc, mesh, matl_model, tol)

  !! Wrap up
  call halt_parallel_communication
  stop status

contains

  subroutine error_exit(errmsg)
    character(len=*), intent(in) :: errmsg
    status = 1
    if (is_IOP) print '("FATAL: ",a)', errmsg
    call halt_parallel_communication
    stop status
  end subroutine error_exit

  !! Initial temperature distribution. Note that x(1) is time.
  function sinusoid(x, p) result (fx)
    real(r8), intent(in) :: x(*), p(*)
    real(r8) :: fx
    fx = 1 + sin(p(2)*x(2)) * sin(p(3)*x(3))
  end function sinusoid

  !! Initialize solver parameter list
  subroutine init_params(params)

    type(parameter_list), intent(out) :: params
    type(parameter_list), pointer :: sublist

    call params%set_path('solver')

    sublist => params%sublist('preconditioner')
    call sublist%set('method','BoomerAMG')
    sublist => sublist%sublist('params')
    call sublist%set('num-cycles', 1)

    sublist => params%sublist('error-norm')
    call sublist%set('rel-t-tol', 1.0d-4)
    call sublist%set('rel-h-tol', 1.0d-4)

    sublist => params%sublist('integrator')
    call sublist%set('nlk-max-iter', 5)
    call sublist%set('nlk-tol', 0.01_r8)

  end subroutine init_params


  !! Tests the HT_2d_solver on a linear problem with Dirichlet boundary conditions
  subroutine test1(disc, mesh, matl_model, tol)

    use idaesol_type, only: SOLVED_TO_TOUT

    type(mfd_2d_disc), target, intent(in) :: disc
    type(unstr_2d_mesh), target, intent(in) :: mesh
    type(material_model), target, intent(in) :: matl_model
    real(r8), intent(in) :: tol

    type(ht_2d_solver), target :: HT_solver, bad_solver
    type(ht_2d_model), target :: HT_model
    type(parameter_list) :: solver_params
    type(parameter_list), pointer :: model_params, sublist
    class(scalar_func), allocatable :: f
    real(r8), allocatable :: Tface(:), Tcell0(:), Tcell1(:)
    character(:), allocatable :: errmsg, string
    character(80) :: metrics(2)
    integer :: stat
    real(r8) :: t, dt, h, max_error, l2_error

    if (is_IOP) print '(/,"Testing sinusoidal problem with mixed BCs")'

    !! 2D HT model parameters
    string = &
    '{"bc": { &
        "left-bottom": { &
          "type": "temperature", &
          "face-set-ids": [1,3], &
          "temp": 1.0}, &
        "right-top": { &
          "type": "flux", &
          "face-set-ids": [2,4], &
          "flux": 0.0}}}'
    call parameter_list_from_json_string(string, model_params, errmsg)
    if (.not. associated(model_params)) call error_exit(errmsg)

    !! Initialize 2D HT model
    call HT_model%init(disc%mesh, matl_model, material_composition_ref(), test_env, model_params, stat, errmsg)
    if (stat /= 0) call error_exit(errmsg)

    !! Define the function f=1+sin(PI/2 x)*sin(PI/2 y)
    call alloc_fptr_scalar_func(f, sinusoid, [0.0_r8, PI/2, PI/2])

    !! Compute cell-based and face-based fields on the mesh.
    allocate(Tcell0(mesh%ncell_onP), Tface(mesh%nface_onP))
    call average_integral(disc, f, Tcell0, Tface)

    !! Initialize 2D HT solver
    call init_params(solver_params)
    call HT_solver%init(test_env, HT_model, solver_params, stat, errmsg)
    if (stat /= 0) call error_exit(errmsg)

    !! Invalid initial-condition parameters are reported to the caller.
    sublist => solver_params%sublist('initial-condition')
    call sublist%set('rel-tol', 0.0_r8)
    call bad_solver%init(test_env, HT_model, solver_params, stat, errmsg)
    if (stat == 0 .or. index(errmsg, '"rel-tol"') == 0) then
      if (is_IOP) print '("ERROR: invalid initial-condition tolerance was accepted")'
      status = 1
      return
    end if

    t = 0.0_r8
    dt = 1E-3_r8
    call HT_solver%set_initial_state(test_env, t, dt, Tcell0, stat, errmsg)
    if (stat /= 0) call error_exit(errmsg)

    !! Run solver
    h = 1E-7_r8
    call HT_solver%integrate(h, stat, tout=0.1_r8)
    if (is_IOP) then
      print '("stat=",i3)', stat
      call HT_solver%write_metrics(metrics)
      print '(a)', metrics
    end if

    if (stat /= SOLVED_TO_TOUT) then
      if (is_IOP) print '("ERROR: failed to integrate to final time")'
      status = 1
      return
    end if

    t = HT_solver%last_time()

    allocate(Tcell1(mesh%ncell_onP))
    call HT_solver%get_cell_temp_soln(Tcell1)

    !! result must match analytical solution
    !!   u(x,y,t) = 1 + sin(PI/2 x)*sin(PI/2 y)*exp(-PI^2/2 t)
    Tcell0 = exp(-0.5_r8*PI**2*t)*(Tcell0-1) - (Tcell1-1)

    max_error = global_maxval(abs(Tcell0))
    Tcell0 = Tcell0 * mesh%volume(1:mesh%ncell_onP)
    l2_error = sqrt(global_dot_product(Tcell0, Tcell0))

    if (is_IOP) print '("l2 error=",es10.3,/,"max error=",es10.3)', l2_error, max_error

    if (max_error > tol) then
      if (is_IOP) print '("ERROR: max error exceeds tolerance; tol=",es9.2)', tol
      status = 1
    end if

  end subroutine test1

end program test_HT_2d_solver_type
