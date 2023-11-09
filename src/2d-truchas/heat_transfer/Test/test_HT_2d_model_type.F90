!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_HT_2d_model_type

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
  use cell_geometry, only: normalized
  use mfd_2d_disc_type
  use HT_2d_model_type
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

  TOL = 1E-10_r8
  eps = 0.0_r8  ! mesh distortion
  xmin = [0.0_r8, 0.0_r8]
  xmax = [1.0_r8, 1.0_r8]
  nx  = [64, 64]

  !! Create the mesh specified by the above input file
  !TODO: breaks if nproc < ncell
  mesh => new_unstr_2d_mesh(xmin, xmax, nx, eps)

  !! Initialize state needed by all tests
  call mfd_disc%init(mesh)
  call init_materials(mesh, matl_db, mmf)

  !! Run test problems
  call test_linear_dir(mfd_disc, mmf, tol)
  call test_linear_flux(mfd_disc, mmf, tol)
  call test_quadratic_dir(mfd_disc, mmf, tol)
  call test_quadratic_flux(mfd_disc, mmf, tol)

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

  !! Tests the HT_2d_model on a linear problem with Dirichlet boundary conditions
  subroutine test_linear_dir(disc, mmf, tol)

    type(mfd_2d_disc), target, intent(in) :: disc
    type(matl_mesh_func), target, intent(in) :: mmf
    real(r8), intent(in) :: tol

    type(HT_2d_model) :: HT_model
    type(parameter_list), pointer :: params
    class(scalar_func), allocatable :: f
    integer :: exps(3,2) = reshape([0,1,0,0,0,1],[3,2])  ! exponents of u(x,t)
    real(r8) :: lcoef(2) = [1.0_r8, 2.0_r8]  ! coefficients of u(x,t)
    real(r8), allocatable, target :: u(:), udot(:), r(:)
    real(r8), pointer :: Hcell(:), ucell(:), uface(:), rcell(:), rface(:)
    character(:), allocatable :: errmsg, string
    integer :: n, stat, max_itr
    real(r8) :: t, rel_tol

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
    call HT_model%init(disc, mmf, params, stat, errmsg)
    if (stat /= 0) call error_exit(errmsg)

    !! Define state variables
    n = HT_model%num_dof()
    allocate(u(n), udot(n), r(n))
    call HT_model%get_cell_heat_view(u, Hcell)
    call HT_model%get_cell_temp_view(u, ucell)
    call HT_model%get_face_temp_view(u, uface)

    !! Define the linear function ax+by
    call alloc_mpoly_scalar_func(f, lcoef, exps)

    !! Compute cell-based and face-based fields on the mesh.
    call average_integral(disc, f, ucell, uface)

    !! Check boundary conditions match face temperatures
    call HT_model%bc_dir%compute(t)
    n = count(HT_model%bc_dir%index <= disc%mesh%nface_onP)
    associate (value => HT_model%bc_dir%value(:n), index => HT_model%bc_dir%index(:n))
      if (global_any(abs(value - uface(index)) > tol)) call error_exit('incorrect BC values')
    end associate

    !! Compute consistent face temperatures
    Hcell = 0.0_r8
    max_itr = 100
    rel_tol = tol
    call HT_model%compute_face_temp(t, u, rel_tol, max_itr, stat, errmsg)
    if (stat/=0) call error_exit(errmsg)

    !! Compute heat transfer residuals
    udot = 0.0_r8
    call HT_model%compute_f(t, u, udot, r)
    call HT_model%get_cell_temp_view(r, rcell)
    call HT_model%get_face_temp_view(r, rface)

    !! Face and cell residuals must be 0
    if (global_any(abs(rface) > tol)) then
      if (is_IOP) print '("ERROR: face residuals are nonzero; tol=",es9.2)', tol
      status = 1
    end if

    if (global_any(abs(rcell) > tol)) then
      if (is_IOP) print '("ERROR: cell residuals are nonzero; tol=",es9.2)', tol
      status = 1
    end if

  end subroutine test_linear_dir


  !! Tests the HT_2d_model on a linear problem with Neumann boundary conditions
  subroutine test_linear_flux(disc, mmf, tol)

    type(mfd_2d_disc), target, intent(in) :: disc
    type(matl_mesh_func), target, intent(in) :: mmf
    real(r8), intent(in) :: tol

    type(HT_2d_model) :: HT_model
    type(parameter_list), pointer :: params
    class(scalar_func), allocatable :: f
    integer :: exps(3,2) = reshape([0,1,0,0,0,1],[3,2])  ! exponents of u(x,t)
    real(r8) :: lcoef(2) = [1.0_r8, 2.0_r8]  ! coefficients of u(x,t)
    real(r8) :: dcoef = 1.0_r8  ! diffusion coefficient
    real(r8), allocatable, target :: u(:), udot(:), r(:)
    real(r8), pointer :: Hcell(:), ucell(:), uface(:), rcell(:), rface(:)
    character(:), allocatable :: errmsg, string
    real(r8) :: t, rel_tol
    integer :: n, j, stat, max_itr

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
    call HT_model%init(disc, mmf, params, stat, errmsg)
    if (stat /= 0) call error_exit(errmsg)

    !! Define state variables
    n = HT_model%num_dof()
    allocate(u(n), udot(n), r(n))
    call HT_model%get_cell_heat_view(u, Hcell)
    call HT_model%get_cell_temp_view(u, ucell)
    call HT_model%get_face_temp_view(u, uface)

    !! Define the linear function ax+by
    call alloc_mpoly_scalar_func(f, lcoef, exps)

    !! Compute cell-based and face-based fields on the mesh.
    call average_integral(disc, f, ucell, uface)

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

    !! Compute consistent face temperatures
    Hcell = 0.0_r8
    max_itr = 100
    rel_tol = tol
    call HT_model%compute_face_temp(t, u, rel_tol, max_itr, stat, errmsg)
    if (stat/=0) call error_exit(errmsg)

    !! Compute heat transfer residuals
    udot = 0.0_r8
    call HT_model%compute_f(t, u, udot, r)
    call HT_model%get_cell_temp_view(r, rcell)
    call HT_model%get_face_temp_view(r, rface)

    !! Face and cell residuals must be 0
    if (global_any(abs(rface) > tol)) then
      if (is_IOP) print '("ERROR: face residuals are nonzero; tol=",es9.2)', tol
      status = 1
    end if

    if (global_any(abs(rcell) > tol)) then
      if (is_IOP) print '("ERROR: cell residuals are nonzero; tol=",es9.2)', tol
      status = 1
    end if

  end subroutine test_linear_flux


  !! Tests the HT_2d_model on a quadratic problem with Dirichlet boundary conditions
  !! The quadratic has the form ax^2+by^2, where a and b are constant
  subroutine test_quadratic_dir(disc, mmf, tol)

    type(mfd_2d_disc), target, intent(in) :: disc
    type(matl_mesh_func), target, intent(in) :: mmf
    real(r8), intent(in) :: tol

    type(HT_2d_model) :: HT_model
    type(parameter_list), pointer :: params
    class(scalar_func), allocatable :: f
    integer :: exps(3,2) = reshape([0,2,0,0,0,2],[3,2])  ! exponents of u(x,t)
    real(r8) :: qcoef(2) = [1.0_r8, 2.0_r8]  ! coefficient of u(x,t)
    real(r8) :: dcoef = 1.0_r8  ! diffusion coefficient
    real(r8), allocatable, target :: u(:), udot(:), r(:)
    real(r8), allocatable :: uface_all(:)
    real(r8), pointer :: Hcell(:), ucell(:), uface(:), Hdot(:), rcell(:), rface(:)
    character(:), allocatable :: errmsg, string
    real(r8) :: t, rel_tol
    integer :: n, max_itr, stat

    if (is_IOP) print '(/,"Testing quadratic problem with Dirichlet BCs")'

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
            "poly-powers": [[0,2,0],[0,0,2]]}}}}'
    call parameter_list_from_json_string(string, params, errmsg)
    if (.not. associated(params)) call error_exit(errmsg)

    !! Initialize 2D HT model
    call HT_model%init(disc, mmf, params, stat, errmsg)
    if (stat/=0) call error_exit(errmsg)

    !! Define state variables
    n = HT_model%num_dof()
    allocate(u(n), udot(n), r(n))
    call HT_model%get_cell_heat_view(u, Hcell)
    call HT_model%get_cell_temp_view(u, ucell)
    call HT_model%get_face_temp_view(u, uface)
    call HT_model%get_cell_heat_view(udot, Hdot)

    !! Define the quadratic function ax^2+by^2
    call alloc_mpoly_scalar_func(f, qcoef, exps)

    !! Compute cell-based and face-based fields on the mesh.
    call average_integral(disc, f, ucell, uface)

    ! TODO: is there a better way?
    ! NOTE: bc_dir%compute(t) computes the BC functions at the center of each face
    !       but average_integral computes their integral over the face divided
    !       by the face length.  For linear functions, these match exactly.
    ! TODO: which is correct for this problem, BC functions correct, or average integral?
    allocate(uface_all(mesh%nface))
    call HT_model%get_face_temp_copy(u, uface_all)
    call mesh%face_imap%gather_offp(uface_all)
    call HT_model%bc_dir%compute(t)
    HT_model%bc_dir%value = uface_all(HT_model%bc_dir%index)
    call HT_model%bc_dir%compute(t)
    associate (value => HT_model%bc_dir%value, index => HT_model%bc_dir%index)
      if (global_any(abs(value - uface_all(index)) > tol)) call error_exit('incorrect BC values')
    end associate

    !! Compute consistent face temperatures
    Hcell = 0.0_r8
    max_itr = 100
    rel_tol = tol
    call HT_model%compute_face_temp(t, u, rel_tol, max_itr, stat, errmsg)
    if (stat/=0) call error_exit(errmsg)

    !! Compute heat transfer residuals
    udot = 0.0_r8
    Hdot = dcoef*sum(2*qcoef)  ! negative of expected cell residuals
    call HT_model%compute_f(t, u, udot, r)
    call HT_model%get_cell_temp_view(r, rcell)
    call HT_model%get_face_temp_view(r, rface)

    !TODO: SWITCH TO TEST LIKE THIS?
    ! INSIST(global_sum(abs(rcell-val*disc%mesh%volume))/disc%mesh%cell_ip%global_size() <= 1E-3)

    !! Face and cell residuals must be 0
    if (global_any(abs(rface) > tol)) then
      if (is_IOP) print '("ERROR: face residuals are nonzero; tol=",es9.2)', tol
      status = 1
    end if

    if (global_any(abs(rcell) > tol)) then
      if (is_IOP) print '("ERROR: cell residuals are nonzero; tol=",es9.2)', tol
      status = 1
    end if

  end subroutine test_quadratic_dir


  !! Tests the HT_2d_model on a quadratic problem with Neumann boundary conditions
  !! The quadratic has the form ax^2+by^2, where a and b are constant
  subroutine test_quadratic_flux(disc, mmf, tol)

    type(mfd_2d_disc), target, intent(in) :: disc
    type(matl_mesh_func), target, intent(in) :: mmf
    real(r8), intent(in) :: tol

    type(HT_2d_model) :: HT_model
    type(parameter_list), pointer :: params
    class(scalar_func), allocatable :: f
    integer :: exps(3,2) = reshape([0,2,0,0,0,2],[3,2])  ! exponents of u(x,t)
    real(r8) :: qcoef(2) = [1.0_r8, 2.0_r8]  ! coefficients of u(x,t)
    real(r8) :: dcoef = 1.0_r8  ! diffusion coefficient
    real(r8), allocatable, target :: u(:), udot(:), r(:)
    real(r8), pointer :: Hcell(:), ucell(:), uface(:), Hdot(:), rcell(:), rface(:)
    character(:), allocatable :: errmsg, string
    real(r8) :: t, rel_tol
    integer :: n, max_itr, stat

    if (is_IOP) print '(/,"Testing quadratic problem with Neumann BCs")'

    t = 0.0_r8

    !! 2D HT model parameters
    string = &
    '{"bc": { &
        "left": { &
          "type": "flux", &
          "face-set-ids": [1], &
          "flux": { &
            "type": "polynomial", &
            "poly-coef": [2.0], &
            "poly-powers": [[0,1,0]]}}, &
        "right": { &
          "type": "flux", &
          "face-set-ids": [2], &
          "flux": { &
            "type": "polynomial", &
            "poly-coef": [-2.0], &
            "poly-powers": [[0,1,0]]}}, &
        "bottom": { &
          "type": "flux", &
          "face-set-ids": [3], &
          "flux": { &
            "type": "polynomial", &
            "poly-coef": [4.0], &
            "poly-powers": [[0,0,1]]}}, &
        "top": { &
          "type": "flux", &
          "face-set-ids": [4], &
          "flux": { &
            "type": "polynomial", &
            "poly-coef": [-4.0], &
            "poly-powers": [[0,0,1]]}}}}'
    call parameter_list_from_json_string(string, params, errmsg)
    if (.not. associated(params)) call error_exit(errmsg)

    !! Initialize 2D HT model
    call HT_model%init(disc, mmf, params, stat, errmsg)
    if (stat /= 0) call error_exit(errmsg)

    !! Define state variables
    n = HT_model%num_dof()
    allocate(u(n), udot(n), r(n))
    call HT_model%get_cell_heat_view(u, Hcell)
    call HT_model%get_cell_temp_view(u, ucell)
    call HT_model%get_face_temp_view(u, uface)
    call HT_model%get_cell_heat_view(udot, Hdot)

    !! Define the quadratic function ax^2+by^2
    call alloc_mpoly_scalar_func(f, qcoef, exps)

    !! Compute cell-based and face-based fields on the mesh.
    call average_integral(disc, f, ucell, uface)

    !! Check boundary conditions match expected value
    block
      logical, allocatable :: mask(:)
      real(r8) :: expected, normal(2), face_centroid(2)
      logical :: failed = .false.
      integer :: i, j, f
      call HT_model%bc_flux%compute(t)
      do i = 1, size(disc%mesh%face_set_id)
        mask = btest(disc%mesh%face_set_mask(HT_model%bc_flux%index), disc%mesh%face_set_id(i))
        associate (index => pack(HT_model%bc_flux%index, mask), &
                   value => pack(HT_model%bc_flux%value, mask))
          if (size(index) < 1) cycle
          normal = normalized(disc%mesh%normal(:,index(1)))
          !! The boundary face fluxes must equal -k*(GRAD u).(face_normal)
          do j = 1, size(index)
            f = index(j)
            face_centroid = 0.5_r8 * sum(disc%mesh%x(:,disc%mesh%fnode(:,f)), dim=2)
            expected = -dcoef*dot_product((2*qcoef)*face_centroid, normal)
            if (abs(value(j) - expected) > tol) failed = .true.
          end do
        end associate
        if (global_any(failed)) call error_exit('incorrect BC values')
      end do
    end block

    !! Compute consistent face temperatures
    Hcell = 0.0_r8
    max_itr = 100
    rel_tol = 1E-10_r8
    call HT_model%compute_face_temp(t, u, rel_tol, max_itr, stat, errmsg)
    if (stat/=0) call error_exit(errmsg)

    !! Compute heat transfer residuals
    udot = 0.0_r8
    Hdot = dcoef*sum(2*qcoef)  ! negative of expected cell residuals
    call HT_model%compute_f(t, u, udot, r)
    call HT_model%get_cell_temp_view(r, rcell)
    call HT_model%get_face_temp_view(r, rface)

    !! Face and cell residuals must be 0
    if (global_any(abs(rface) > tol)) then
      if (is_IOP) print '("ERROR: face residuals are nonzero; tol=",es9.2)', tol
      status = 1
    end if

    if (global_any(abs(rcell) > tol)) then
      if (is_IOP) print '("ERROR: cell residuals are nonzero; tol=",es9.2)', tol
      status = 1
    end if

  end subroutine test_quadratic_flux

end program test_HT_2d_model_type
