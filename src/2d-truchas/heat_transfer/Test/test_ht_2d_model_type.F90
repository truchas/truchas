!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_HT_2d_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication
  use fhypre, only: fhypre_initialize
  use parameter_list_type
  use parameter_list_json
  use unstr_2d_mesh_factory
  use material_model_type
  use scalar_func_factories
  use cell_geometry, only: normalized
  use mfd_2d_disc_type
  use ht_2d_model_type
  use ht_2d_ic_solver_type
  use ht_2d_vector_type
  use bitfield_type
  use test_ht_2d_common
  implicit none

  type(unstr_2d_mesh), pointer :: mesh
  type(mfd_2d_disc), target :: mfd_disc
  type(material_model), target :: matl_model
  real(r8) :: xmin(2), xmax(2), tol, eps
  integer  :: nx(2)
  integer :: status = 0

  !! Initialize MPI and other base stuff that Truchas depends on
  call init_parallel_communication
  call fhypre_initialize
  call init_test_environment('test_ht_2d_model_type.log')

  TOL = 1E-10_r8
  eps = 0.0_r8  ! mesh distortion
  xmin = [0.0_r8, 0.0_r8]
  xmax = [1.0_r8, 1.0_r8]
  nx  = [64, 64]

  !! Create the mesh specified by the above input file
  !TODO: breaks if nproc < ncell
  mesh => new_unstr_2d_mesh(test_env, xmin, xmax, nx, eps)

  !! Initialize state needed by all tests
  call mfd_disc%init(mesh)
  call init_materials(mesh, matl_model)

  !! Run test problems
  call test_linear_dir(mfd_disc, matl_model, tol)
  call test_linear_flux(mfd_disc, matl_model, tol)
  !! Quadratic boundary-condition verification remains deferred until its
  !! expected discrete convergence behavior is established.

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

  subroutine check_prop_extent(model, tol)
    type(ht_2d_model), intent(in) :: model
    real(r8), intent(in) :: tol

    real(r8), allocatable :: state(:), value(:), value0(:), deriv(:), deriv0(:)
    integer :: ncell

    ncell = model%mesh%ncell_onP
    allocate(state(ncell+1), value(ncell), value0(ncell), deriv(ncell), deriv0(ncell))
    state = 1.0_r8

    call model%H_of_T%compute_value(state(:ncell), value0)
    call model%H_of_T%compute_value(state, value)
    if (global_any(abs(value - value0) > tol)) then
      call error_exit('property value changed when state included an extra cell')
    end if

    call model%H_of_T%compute_deriv(state(:ncell), 1, deriv0)
    call model%H_of_T%compute_deriv(state, 1, deriv)
    if (global_any(abs(deriv - deriv0) > tol)) then
      call error_exit('property derivative changed when state included an extra cell')
    end if
  end subroutine check_prop_extent

  !! Tests the HT_2d_model on a linear problem with Dirichlet boundary conditions
  subroutine test_linear_dir(disc, matl_model, tol)

    type(mfd_2d_disc), target, intent(in) :: disc
    type(material_model), target, intent(in) :: matl_model
    real(r8), intent(in) :: tol

    type(ht_2d_model), target :: HT_model
    type(ht_2d_ic_solver) :: ic
    type(parameter_list), pointer :: params
    type(parameter_list), target :: ic_params
    class(scalar_func), allocatable :: f
    integer :: exps(3,2) = reshape([0,1,0,0,0,1],[3,2])  ! exponents of u(x,t)
    real(r8) :: lcoef(2) = [1.0_r8, 2.0_r8]  ! coefficients of u(x,t)
    type(ht_2d_vector) :: u, udot, r
    real(r8), allocatable :: Tcell(:), Tface(:)
    character(:), allocatable :: errmsg, string
    integer :: n, stat, max_itr
    real(r8) :: t, dt, rel_tol

    if (is_IOP) print '(/,"Testing linear problem with Dirichlet BCs")'

    t = 0.0_r8

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
    call parameter_list_from_json_string(string, params, errmsg)
    if (.not. associated(params)) call error_exit(errmsg)

    !! Initialize 2D HT model
    call HT_model%init(disc%mesh, matl_model, material_composition_ref(), test_env, params, stat, errmsg)
    if (stat /= 0) call error_exit(errmsg)
    call check_prop_extent(HT_model, tol)

    call u%init(disc%mesh)
    call udot%init(u)
    call r%init(u)

    !! Define the linear function ax+by
    call alloc_mpoly_scalar_func(f, lcoef, exps)

    !! Compute cell-based and face-based fields on the mesh.
    allocate(Tcell(disc%mesh%ncell_onP), Tface(disc%mesh%nface_onP))
    call average_integral(disc, f, Tcell, Tface)

    !! Check boundary conditions match face temperatures
    call HT_model%bc_dir%compute(t)
    n = count(HT_model%bc_dir%index <= disc%mesh%nface_onP)
    associate (value => HT_model%bc_dir%value(:n), index => HT_model%bc_dir%index(:n))
      if (global_any(abs(value - Tface(index)) > tol)) call error_exit('incorrect BC values')
    end associate

    !! Compute a consistent vector state, including face temperatures.
    dt = 1.0e-3_r8
    max_itr = 100
    rel_tol = tol
    call ic_params%set('dt', dt)
    call ic_params%set('rel-tol', rel_tol)
    call ic_params%set('max-iter', max_itr)
    call ic%init(HT_model, ic_params)
    call ic%compute(test_env, t, Tcell, u, udot, stat, errmsg)
    if (stat/=0) call error_exit(errmsg)

    !! Compute heat transfer residuals
    call udot%setval(0.0_r8)
    call r%setval(0.0_r8)
    call HT_model%residual(t, u, udot, r)

    !! Face and cell residuals must be 0
    if (global_any(abs(r%tf(:disc%mesh%nface_onP)) > 20.0_r8*tol)) then
      if (is_IOP) print '("ERROR: face residuals are nonzero; tol=",es9.2)', tol
      if (is_IOP) print '("  max |rface| = ",es12.4)', global_maxval(maxval(abs(r%tf(:disc%mesh%nface_onP))))
      status = 1
    end if

    if (global_any(abs(r%tc(:disc%mesh%ncell_onP)) > 20.0_r8*tol)) then
      if (is_IOP) print '("ERROR: cell residuals are nonzero; tol=",es9.2)', tol
      if (is_IOP) print '("  max |rcell| = ",es12.4)', global_maxval(maxval(abs(r%tc(:disc%mesh%ncell_onP))))
      status = 1
    end if

  end subroutine test_linear_dir


  !! Tests the HT_2d_model on a linear problem with Neumann boundary conditions
  subroutine test_linear_flux(disc, matl_model, tol)

    type(mfd_2d_disc), target, intent(in) :: disc
    type(material_model), target, intent(in) :: matl_model
    real(r8), intent(in) :: tol

    type(ht_2d_model), target :: HT_model
    type(ht_2d_ic_solver) :: ic
    type(parameter_list), pointer :: params
    type(parameter_list), target :: ic_params
    class(scalar_func), allocatable :: f
    integer :: exps(3,2) = reshape([0,1,0,0,0,1],[3,2])  ! exponents of u(x,t)
    real(r8) :: lcoef(2) = [1.0_r8, 2.0_r8]  ! coefficients of u(x,t)
    real(r8) :: dcoef = 1.0_r8  ! diffusion coefficient
    type(ht_2d_vector) :: u, udot, r
    real(r8), allocatable :: Tcell(:), Tface(:)
    character(:), allocatable :: errmsg, string
    real(r8) :: t, dt, rel_tol
    integer :: j, stat, max_itr

    if (is_IOP) print '(/,"Testing linear problem with Neumann BCs")'

    t = 0.0_r8

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
          "flux": -2.0 }}}'
    call parameter_list_from_json_string(string, params, errmsg)
    if (.not. associated(params)) call error_exit(errmsg)

    !! Initialize 2D HT model
    call HT_model%init(disc%mesh, matl_model, material_composition_ref(), test_env, params, stat, errmsg)
    if (stat /= 0) call error_exit(errmsg)
    call check_prop_extent(HT_model, tol)

    call u%init(disc%mesh)
    call udot%init(u)
    call r%init(u)

    !! Define the linear function ax+by
    call alloc_mpoly_scalar_func(f, lcoef, exps)

    !! Compute cell-based and face-based fields on the mesh.
    allocate(Tcell(disc%mesh%ncell_onP), Tface(disc%mesh%nface_onP))
    call average_integral(disc, f, Tcell, Tface)

    !! Check boundary conditions match expected value
    block
      logical, allocatable :: mask(:)
      real(r8) :: expected, normal(2)
      call HT_model%bc_flux%compute(t)
      do j = 1, size(disc%mesh%face_set_id)
        mask = btest(disc%mesh%face_set_mask(HT_model%bc_flux%index), disc%mesh%face_set_id(j))
        associate (index => pack(HT_model%bc_flux%index, mask), &
                   value => pack(HT_model%bc_flux%value, mask))
          if (size(index) < 1) cycle
          normal = normalized(disc%mesh%normal(:,index(1)))
          !! The boundary face fluxes must equal -k*(GRAD u).(face_normal)
          expected = -dcoef*dot_product(lcoef, normal)
          if (global_any(abs(value - expected) > tol)) call error_exit('incorrect BC values')
        end associate
      end do
    end block

    !! Compute a consistent vector state, including face temperatures.
    dt = 1.0e-3_r8
    max_itr = 100
    rel_tol = tol
    call ic_params%set('dt', dt)
    call ic_params%set('rel-tol', rel_tol)
    call ic_params%set('max-iter', max_itr)
    call ic%init(HT_model, ic_params)
    call ic%compute(test_env, t, Tcell, u, udot, stat, errmsg)
    if (stat/=0) call error_exit(errmsg)

    !! Compute heat transfer residuals
    call udot%setval(0.0_r8)
    call r%setval(0.0_r8)
    call HT_model%residual(t, u, udot, r)

    !! Face and cell residuals must be 0
    if (global_any(abs(r%tf(:disc%mesh%nface_onP)) > 20.0_r8*tol)) then
      if (is_IOP) print '("ERROR: face residuals are nonzero; tol=",es9.2)', tol
      if (is_IOP) print '("  max |rface| = ",es12.4)', global_maxval(maxval(abs(r%tf(:disc%mesh%nface_onP))))
      status = 1
    end if

    if (global_any(abs(r%tc(:disc%mesh%ncell_onP)) > 20.0_r8*tol)) then
      if (is_IOP) print '("ERROR: cell residuals are nonzero; tol=",es9.2)', tol
      if (is_IOP) print '("  max |rcell| = ",es12.4)', global_maxval(maxval(abs(r%tc(:disc%mesh%ncell_onP))))
      status = 1
    end if

  end subroutine test_linear_flux


  !! Tests the HT_2d_model on a quadratic problem with Dirichlet boundary conditions
  !! The quadratic has the form ax^2+by^2, where a and b are constant

end program test_HT_2d_model_type
