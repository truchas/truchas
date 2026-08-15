!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_HT_2d_precon_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication
  use fhypre, only: fhypre_initialize
  use truchas_env, only: prefix, overwrite_output
  use truchas_logging_services
  use parameter_list_type
  use parameter_list_json
  use unstr_2d_mesh_factory
  use matl_mesh_func_type
  use material_model_type
  !use source_mesh_function
  use scalar_func_factories
  use mfd_2d_disc_type
  use ht_2d_model_type
  use ht_2d_precon_type
  use ht_2d_vector_type
  use bitfield_type
  use test_ht_2d_common
  implicit none

  type(unstr_2d_mesh), pointer :: mesh
  type(mfd_2d_disc), target :: mfd_disc
  type(material_model), target :: matl_model
  type(matl_mesh_func), target :: mmf
  real(r8) :: xmin(2), xmax(2), tol, eps
  integer  :: nx(2)
  integer :: status = 0

  !! Initialize MPI and other base stuff that Truchas depends on
  call init_parallel_communication
  call fhypre_initialize
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
  call init_materials(mesh, matl_model, mmf)

  !! Run test problems
  call test_linear_dir(mfd_disc, matl_model, tol)
  call test_linear_flux(mfd_disc, matl_model, tol)

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

  !! Tests the HT_2d_precon on a linear problem with Dirichlet boundary conditions
  subroutine test_linear_dir(disc, matl_model, tol)

    type(mfd_2d_disc), target, intent(in) :: disc
    type(material_model), target, intent(in) :: matl_model
    real(r8), intent(in) :: tol

    type(ht_2d_precon) :: HT_precon
    type(ht_2d_model), target :: HT_model
    type(parameter_list), pointer :: params, sublist
    class(scalar_func), allocatable :: f
    integer :: exps(3,2) = reshape([0,1,0,0,0,1],[3,2])  ! exponents of u(x,t)
    real(r8) :: lcoef(2) = [1.0_r8, 2.0_r8]  ! coefficients of u(x,t)
    type(ht_2d_vector) :: state, u, udot, r
    real(r8), allocatable :: func_state(:,:), Hcell(:)
    real(r8), allocatable :: Tcell(:), Tface(:)
    character(:), allocatable :: errmsg, string
    integer :: stat
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
        "params": { &
          "num-cycles": 25}}}'
    call parameter_list_from_json_string(string, params, errmsg)
    if (.not. associated(params)) call error_exit(errmsg)

    !! Initialize 2D HT model
    sublist => params%sublist('model')
    call HT_model%init(disc, matl_model, material_composition_ref(), sublist, stat, errmsg)
    if (stat /= 0) call error_exit(errmsg)

    !! Initialize 2D HT preconditioner
    sublist => params%sublist('preconditioner')
    call HT_precon%init(HT_model, sublist)

    !! Define state variables
    call state%init(disc%mesh)
    call u%init(state)
    call udot%init(state)
    call r%init(state)

    !! Compute RHS of Jacobian system
    call state%setval(0.0_r8)
    call udot%setval(0.0_r8)
    call r%setval(0.0_r8)
    call HT_model%residual(t, state, udot, r)

    !! Preconditioner fully solves for steady state solution
    dt = huge(0.0_r8)  !TODO: test finite dt
    call HT_precon%compute(t, state, dt)
    call u%copy(r)
    call u%scale(-1.0_r8)
    call HT_precon%apply(t, state, u)

    !! Check residual
    call HT_model%residual(t, u, udot, r)

    if (global_any(r%hc(:disc%mesh%ncell_onP) > tol)) then
      if (is_IOP) print '("ERROR: cell enthalpy residual is nonzero; tol=",es9.2)', tol
      status = 1
    end if
    if (global_any(r%tc(:disc%mesh%ncell_onP) > tol)) then
      if (is_IOP) print '("ERROR: cell temp residual is nonzero; tol=",es9.2)', tol
      status = 1
    end if
    if (global_any(r%tf(:disc%mesh%nface_onP) > tol)) then
      if (is_IOP) print '("ERROR: face temp residual is nonzero; tol=",es9.2)', tol
      status = 1
    end if

    !! Expected cell and face temperature fields.
    call alloc_mpoly_scalar_func(f, lcoef, exps)
    allocate(Tcell(disc%mesh%ncell_onP), Tface(disc%mesh%nface_onP))
    call average_integral(disc, f, Tcell, Tface)

    !! Expected cell enthalpy field.
    allocate(Hcell(disc%mesh%ncell))
    allocate(func_state(disc%mesh%ncell,0:0))
    func_state(:disc%mesh%ncell_onP,0) = Tcell
    call disc%mesh%cell_imap%gather_offp(func_state(:,0))
    call HT_model%H_of_T%compute_value(func_state, Hcell)
    deallocate(func_state)

    !! Check solution
    if (global_any(abs(Hcell(:disc%mesh%ncell_onP)-u%hc(:disc%mesh%ncell_onP)) > tol)) then
      if (is_IOP) print '("ERROR: cell enthalpy exceeds expected value; tol=",es9.2)', tol
      status = 1
    end if
    if (global_any(abs(Tcell-u%tc(:disc%mesh%ncell_onP)) > tol)) then
      if (is_IOP) print '("ERROR: cell temp exceeds expected value; tol=",es9.2)', tol
      status = 1
    end if
    if (global_any(abs(Tface-u%tf(:disc%mesh%nface_onP)) > tol)) then
      if (is_IOP) print '("ERROR: face temp exceeds expected value; tol=",es9.2)', tol
      status = 1
    end if

  end subroutine test_linear_dir


  !! Tests the HT_2d_precon on a linear problem with Neumann boundary conditions
  subroutine test_linear_flux(disc, matl_model, tol)

    type(mfd_2d_disc), target, intent(in) :: disc
    type(material_model), target, intent(in) :: matl_model
    real(r8), intent(in) :: tol

    type(ht_2d_precon) :: HT_precon
    type(ht_2d_model), target :: HT_model
    type(parameter_list), pointer :: params, sublist
    class(scalar_func), allocatable :: f
    integer :: exps(3,2) = reshape([0,1,0,0,0,1],[3,2])  ! exponents of u(x,t)
    real(r8) :: lcoef(2) = [1.0_r8, 2.0_r8]  ! coefficients of u(x,t)
    type(ht_2d_vector) :: state, u, udot, r
    real(r8), allocatable :: func_state(:,:), Hcell(:)
    real(r8), allocatable :: Tcell(:), Tface(:)
    character(:), allocatable :: errmsg, string
    integer :: stat
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
        "params": { &
          "num-cycles": 25}}}'
    call parameter_list_from_json_string(string, params, errmsg)
    if (.not. associated(params)) call error_exit(errmsg)

    !! Initialize 2D HT model
    sublist => params%sublist('model')
    call HT_model%init(disc, matl_model, material_composition_ref(), sublist, stat, errmsg)
    if (stat /= 0) call error_exit(errmsg)

    !! Initialize 2D HT preconditioner
    sublist => params%sublist('preconditioner')
    call HT_precon%init(HT_model, sublist)

    !! Define state variables
    call state%init(disc%mesh)
    call u%init(state)
    call udot%init(state)
    call r%init(state)

    !! Compute RHS of Jacobian system
    call state%setval(0.0_r8)
    call udot%setval(0.0_r8)
    call r%setval(0.0_r8)
    call HT_model%residual(t, state, udot, r)

    !! Preconditioner fully solves for steady state solution
    dt = huge(0.0_r8)  !TODO: test finite dt
    call HT_precon%compute(t, state, dt)
    call u%copy(r)
    call u%scale(-1.0_r8)
    call HT_precon%apply(t, state, u)

    !! Check residual
    call HT_model%residual(t, u, udot, r)

    if (global_any(r%hc(:disc%mesh%ncell_onP) > tol)) then
      if (is_IOP) print '("ERROR: cell enthalpy residual is nonzero; tol=",es9.2)', tol
      status = 1
    end if
    if (global_any(r%tc(:disc%mesh%ncell_onP) > tol)) then
      if (is_IOP) print '("ERROR: cell temp residual is nonzero; tol=",es9.2)', tol
      status = 1
    end if
    if (global_any(r%tf(:disc%mesh%nface_onP) > tol)) then
      if (is_IOP) print '("ERROR: face temp residual is nonzero; tol=",es9.2)', tol
      status = 1
    end if

    !! Expected cell and face temperature fields.
    call alloc_mpoly_scalar_func(f, lcoef, exps)
    allocate(Tcell(disc%mesh%ncell_onP), Tface(disc%mesh%nface_onP))
    call average_integral(disc, f, Tcell, Tface)

    !! The solution of this steady-state problem with flux BC is determined
    !! only up to an additive constant.  Shift the solution error so that the
    !! cell temperature error is zero on the first cell of the rank 0 process.
    if (this_PE == 1) then
      shift = Tcell(1) - u%tc(1)
    end if
    call broadcast(shift)

    Tcell = Tcell - shift
    Tface = Tface - shift

    !! Expected cell enthalpy field.
    allocate(Hcell(disc%mesh%ncell))
    allocate(func_state(disc%mesh%ncell,0:0))
    func_state(:disc%mesh%ncell_onP,0) = Tcell
    call disc%mesh%cell_imap%gather_offp(func_state(:,0))
    call HT_model%H_of_T%compute_value(func_state, Hcell)
    deallocate(func_state)

    !! Check solution
    if (global_any(abs(Hcell(:disc%mesh%ncell_onP)-u%hc(:disc%mesh%ncell_onP)) > tol)) then
      if (is_IOP) print '("ERROR: cell enthalpy exceeds expected value; tol=",es9.2)', tol
      status = 1
    end if
    if (global_any(abs(Tcell-u%tc(:disc%mesh%ncell_onP)) > tol)) then
      if (is_IOP) print '("ERROR: cell temp exceeds expected value; tol=",es9.2)', tol
      status = 1
    end if
    if (global_any(abs(Tface-u%tf(:disc%mesh%nface_onP)) > tol)) then
      if (is_IOP) print '("ERROR: face temp exceeds expected value; tol=",es9.2)', tol
      status = 1
    end if

  end subroutine test_linear_flux

end program test_HT_2d_precon_type
