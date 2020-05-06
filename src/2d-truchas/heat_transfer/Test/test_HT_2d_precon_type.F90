!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_HT_2d_precon_type

#ifdef NAGFOR
  use,intrinsic :: f90_unix, only: exit
#endif
  use kinds, only: r8
  use pgslib_module, only: PGSLib_CL_MAX_TOKEN_LENGTH
  use parallel_util_module, only: parallel_init
  use parallel_communication
  use truchas_env, only: prefix
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
  use HT_2d_precon_type
  use bitfield_type
  use test_ht_2d_common
  implicit none

  character(PGSLib_CL_MAX_TOKEN_LENGTH), pointer :: argv(:) => null()

  type(unstr_2d_mesh), pointer :: mesh
  type(mfd_2d_disc), target :: mfd_disc
  type(material_database), target :: matl_db
  type(matl_mesh_func), target :: mmf
  real(r8) :: xmin(2), xmax(2), tol, eps
  integer  :: nx(2)
  integer :: status = 0

  !! Initialize MPI and other base stuff that Truchas depends on
  call parallel_init(argv)
  call init_parallel_communication
  prefix='run'  ! TLS will write to 'run.log'
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

  !! Tests the HT_2d_precon on a linear problem with Dirichlet boundary conditions
  subroutine test_linear_dir(disc, mmf, tol)

    type(mfd_2d_disc), target, intent(in) :: disc
    type(matl_mesh_func), target, intent(in) :: mmf
    real(r8), intent(in) :: tol

    type(HT_2d_precon) :: HT_precon
    type(HT_2d_model), target :: HT_model
    type(parameter_list), pointer :: params, sublist
    class(scalar_func), allocatable :: f
    integer :: exps(3,2) = reshape([0,1,0,0,0,1],[3,2])  ! exponents of u(x,t)
    real(r8) :: lcoef(2) = [1.0_r8, 2.0_r8]  ! coefficients of u(x,t)
    real(r8), allocatable, target :: u(:), udot(:), r(:)
    real(r8), allocatable :: state(:,:), Hcell(:)
    real(r8), pointer :: Tcell(:), Tface(:), view(:)
    character(:), allocatable :: errmsg, string
    integer :: n, stat
    real(r8) :: t, dt

    if (is_IOP) print '(/,"Testing preconditioner with Dirichlet BCs")'

    t = 0.0_r8

    !! Define problem parameters
    string = &
    '{"model": { &
        "bc": { &
          "all-sides": { &
            "type": "temperature", &
            "face-set-ids": [1,2,3,4], &
            "temp": { &
              "type": "polynomial", &
              "poly-coef": [1.0, 2.0], &
              "poly-powers": [[0,1,0],[0,0,1]]}}}}, &
      "preconditioner": { &
        "method": "BoomerAMG", &
        "parameters": { &
          "num-cycles": 25}}}'
    call parameter_list_from_json_string(string, params, errmsg)
    if (.not. associated(params)) call error_exit(errmsg)

    !! Initialize 2D HT model
    sublist => params%sublist('model')
    call HT_model%init(disc, mmf, sublist, stat, errmsg)
    if (stat /= 0) call error_exit(errmsg)

    !! Initialize 2D HT preconditioner
    sublist => params%sublist('preconditioner')
    call HT_precon%init(HT_model, sublist)

    !! Define state variables
    n = HT_model%num_dof()
    allocate(u(n), udot(n), r(n))

    !! Compute RHS of Jacobian system
    u = 0.0_r8
    udot = 0.0_r8
    call HT_model%compute_f(t, u, udot, r)

    !! Preconditioner fully solves for steady state solution
    dt = huge(0.0_r8)  !TODO: test finite dt
    call HT_precon%compute(t, u, dt)
    u = -r
    call HT_precon%apply(u)

    !! Check residual
    call HT_model%compute_f(t, u, udot, r)

    call HT_model%get_cell_heat_view(r, view)
    if (global_any(view > tol)) then
      if (is_IOP) print '("ERROR: cell enthalpy residual is nonzero; tol=",es9.2)', tol
      status = 1
    end if
    call HT_model%get_cell_temp_view(r, view)
    if (global_any(view > tol)) then
      if (is_IOP) print '("ERROR: cell temp residual is nonzero; tol=",es9.2)', tol
      status = 1
    end if
    call HT_model%get_face_temp_view(r, view)
    if (global_any(view > tol)) then
      if (is_IOP) print '("ERROR: face temp residual is nonzero; tol=",es9.2)', tol
      status = 1
    end if

    !! Expected cell and face temperature fields.
    call alloc_mpoly_scalar_func(f, lcoef, exps)
    call HT_model%get_cell_temp_view(r, Tcell)
    call HT_model%get_face_temp_view(r, Tface)
    call average_integral(disc, f, Tcell, Tface)

    !! Expected cell enthalpy field.
    allocate(Hcell(disc%mesh%ncell))
    call HT_model%new_state_array(r, state)
    call HT_model%H_of_T%compute_value(state, Hcell)
    deallocate(state)

    !! Check solution
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

  end subroutine test_linear_dir


  !! Tests the HT_2d_precon on a linear problem with Neumann boundary conditions
  subroutine test_linear_flux(disc, mmf, tol)

    type(mfd_2d_disc), target, intent(in) :: disc
    type(matl_mesh_func), target, intent(in) :: mmf
    real(r8), intent(in) :: tol

    type(HT_2d_precon) :: HT_precon
    type(HT_2d_model), target :: HT_model
    type(parameter_list), pointer :: params, sublist
    class(scalar_func), allocatable :: f
    integer :: exps(3,2) = reshape([0,1,0,0,0,1],[3,2])  ! exponents of u(x,t)
    real(r8) :: lcoef(2) = [1.0_r8, 2.0_r8]  ! coefficients of u(x,t)
    real(r8), allocatable, target :: u(:), udot(:), r(:)
    real(r8), allocatable :: state(:,:), Hcell(:)
    real(r8), pointer :: Tcell(:), Tface(:), view(:)
    character(:), allocatable :: errmsg, string
    integer :: n, stat
    real(r8) :: t, dt, shift

    if (is_IOP) print '(/,"Testing preconditioner with Neumann BCs")'

    t = 0.0_r8

    !! Define problem parameters
    string = &
    '{"model": { &
        "bc": { &
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
             "flux": -2.0 }}}, &
      "preconditioner": { &
        "method": "BoomerAMG", &
        "parameters": { &
          "num-cycles": 25}}}'
    call parameter_list_from_json_string(string, params, errmsg)
    if (.not. associated(params)) call error_exit(errmsg)

    !! Initialize 2D HT model
    sublist => params%sublist('model')
    call HT_model%init(disc, mmf, sublist, stat, errmsg)
    if (stat /= 0) call error_exit(errmsg)

    !! Initialize 2D HT preconditioner
    sublist => params%sublist('preconditioner')
    call HT_precon%init(HT_model, sublist)

    !! Define state variables
    n = HT_model%num_dof()
    allocate(u(n), udot(n), r(n))

    !! Compute RHS of Jacobian system
    u = 0.0_r8
    udot = 0.0_r8
    call HT_model%compute_f(t, u, udot, r)

    !! Preconditioner fully solves for steady state solution
    dt = huge(0.0_r8)  !TODO: test finite dt
    call HT_precon%compute(t, u, dt)
    u = -r
    call HT_precon%apply(u)

    !! Check residual
    call HT_model%compute_f(t, u, udot, r)

    call HT_model%get_cell_heat_view(r, view)
    if (global_any(view > tol)) then
      if (is_IOP) print '("ERROR: cell enthalpy residual is nonzero; tol=",es9.2)', tol
      status = 1
    end if
    call HT_model%get_cell_temp_view(r, view)
    if (global_any(view > tol)) then
      if (is_IOP) print '("ERROR: cell temp residual is nonzero; tol=",es9.2)', tol
      status = 1
    end if
    call HT_model%get_face_temp_view(r, view)
    if (global_any(view > tol)) then
      if (is_IOP) print '("ERROR: face temp residual is nonzero; tol=",es9.2)', tol
      status = 1
    end if

    !! Expected cell and face temperature fields.
    call alloc_mpoly_scalar_func(f, lcoef, exps)
    call HT_model%get_cell_temp_view(r, Tcell)
    call HT_model%get_face_temp_view(r, Tface)
    call average_integral(disc, f, Tcell, Tface)

    !! The solution of this steady-state problem with flux BC is determined
    !! only up to an additive constant.  Shift the solution error so that the
    !! cell temperature error is zero on the first cell of the rank 0 process.
    if (this_PE == 1) then
      call HT_model%get_cell_temp_view(u, view)
      shift = Tcell(1) - view(1)
    end if
    call broadcast(shift)

    Tcell = Tcell - shift
    Tface = Tface - shift

    !! Expected cell enthalpy field.
    allocate(Hcell(disc%mesh%ncell))
    call HT_model%new_state_array(r, state)
    call HT_model%H_of_T%compute_value(state, Hcell)
    deallocate(state)

    !! Check solution
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

  end subroutine test_linear_flux

end program test_HT_2d_precon_type
