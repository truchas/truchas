!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_HT_2d_solver_initial_state

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication
  use fhypre, only: fhypre_initialize
  use truchas_env, only: prefix, overwrite_output
  use truchas_logging_services
  use parameter_list_type
  use parameter_list_json
  use unstr_2d_mesh_factory
  use material_model_type
  use material_database_type
  use material_factory, only: load_material_database
  use material_utilities, only: add_enthalpy_prop
  use material_composition_type
  !use source_mesh_function
  use scalar_func_factories
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
  call init_materials(mesh, matl_model)

  !! Run test problems
  call test_linear_dir(mfd_disc, matl_model, tol)
  call test_linear_flux(mfd_disc, matl_model, tol)
  call diagnose_multimaterial_dirichlet(mfd_disc, tol)

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

  !! Tests consistent initial-state construction for a linear Dirichlet problem.
  subroutine test_linear_dir(disc, matl_model, tol)

    type(mfd_2d_disc), target, intent(in) :: disc
    type(material_model), target, intent(in) :: matl_model
    real(r8), intent(in) :: tol

    type(ht_2d_ic_solver) :: ic
    type(ht_2d_model), target :: HT_model
    type(parameter_list), pointer :: model_params
    class(scalar_func), allocatable :: f
    integer :: exps(3,2) = reshape([0,1,0,0,0,1],[3,2])  ! exponents of u(x,t)
    real(r8) :: lcoef(2) = [1.0_r8, 2.0_r8]  ! coefficients of u(x,t)
    type(ht_2d_vector) :: u, udot
    real(r8), allocatable :: state(:,:), Hcell(:), Tcell(:), Tface(:)
    character(:), allocatable :: errmsg, string
    integer :: stat, max_itr
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
    call HT_model%init(disc, matl_model, material_composition_ref(), model_params, stat, errmsg)
    if (stat /= 0) call error_exit(errmsg)

    call ic%init(HT_model)
    call u%init(disc%mesh)
    call udot%init(u)

    !! Define the linear function ax+by
    call alloc_mpoly_scalar_func(f, lcoef, exps)

    !! Compute cell-based and face-based fields on the mesh.
    allocate(Tcell(disc%mesh%ncell_onP), Tface(disc%mesh%nface_onP))
    call average_integral(disc, f, Tcell, Tface)

    !! Expected cell enthalpy
    allocate(Hcell(disc%mesh%ncell_onP))
    allocate(state(disc%mesh%ncell,0:0))
    state(:disc%mesh%ncell_onP,0) = Tcell
    call disc%mesh%cell_imap%gather_offp(state(:,0))
    call HT_model%H_of_T%compute_value(state, Hcell)
    deallocate(state)

    !! Compute consistent u and udot
    max_itr = 100
    rel_tol = tol * 1E-5_r8 !TODO: good practice?
    call ic%compute(t, dt, Tcell, u, udot, rel_tol, max_itr, stat, errmsg)
    if (stat/=0) call error_exit(errmsg)

    !! u must match expected values
    if (global_any(abs(Hcell(:disc%mesh%ncell_onP)-u%hc(:disc%mesh%ncell_onP)) > tol)) then
      if (is_IOP) print '("ERROR: cell enthalpy exceeds expected value; tol=",es9.2)', tol
      status = 1
    end if
    if (global_any(abs(Tcell-u%tc(:disc%mesh%ncell_onP)) /= 0.0_r8)) then
      if (is_IOP) print '("ERROR: cell temp exceeds expected value; tol=",es9.2)', tol
      status = 1
    end if
    if (global_any(abs(Tface-u%tf(:disc%mesh%nface_onP)) > tol)) then
      if (is_IOP) print '("ERROR: face temp exceeds expected value; tol=",es9.2)', tol
      status = 1
    end if

    !! udot must be 0
    if (global_any(abs(udot%hc(:disc%mesh%ncell_onP)) > tol)) then
      if (is_IOP) print '("ERROR: Hdot is nonzero; tol=",es9.2)', tol
      status = 1
    end if
    if (global_any(abs(udot%tc(:disc%mesh%ncell_onP)) > tol)) then
      if (is_IOP) print '("ERROR: cell temp derivative is nonzero; tol=",es9.2)', tol
      status = 1
    end if
    if (global_any(abs(udot%tf(:disc%mesh%nface_onP)) > tol)) then
      if (is_IOP) print '("ERROR: face temp derivative is nonzero; tol=",es9.2)', tol
      status = 1
    end if

  end subroutine test_linear_dir


  !! Tests consistent initial-state construction for a linear flux problem.
  subroutine test_linear_flux(disc, matl_model, tol)

    type(mfd_2d_disc), target, intent(in) :: disc
    type(material_model), target, intent(in) :: matl_model
    real(r8), intent(in) :: tol

    type(ht_2d_ic_solver) :: ic
    type(ht_2d_model), target :: HT_model
    type(parameter_list), pointer :: model_params
    class(scalar_func), allocatable :: f
    integer :: exps(3,2) = reshape([0,1,0,0,0,1],[3,2])  ! exponents of u(x,t)
    real(r8) :: lcoef(2) = [1.0_r8, 2.0_r8]  ! coefficients of u(x,t)
    real(r8) :: dcoef = 1.0_r8  ! diffusion coefficient
    type(ht_2d_vector) :: u, udot
    real(r8), allocatable :: state(:,:), Hcell(:), Tcell(:), Tface(:)
    character(:), allocatable :: errmsg, string
    real(r8) :: t, dt, rel_tol
    integer :: stat, max_itr

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
    call HT_model%init(disc, matl_model, material_composition_ref(), model_params, stat, errmsg)
    if (stat /= 0) call error_exit(errmsg)

    call ic%init(HT_model)
    call u%init(disc%mesh)
    call udot%init(u)

    !! Define the linear function ax+by
    call alloc_mpoly_scalar_func(f, lcoef, exps)

    !! Compute cell-based and face-based fields on the mesh.
    allocate(Tcell(disc%mesh%ncell_onP), Tface(disc%mesh%nface_onP))
    call average_integral(disc, f, Tcell, Tface)

    !! Expected cell enthalpy
    allocate(Hcell(disc%mesh%ncell_onP))
    allocate(state(disc%mesh%ncell,0:0))
    state(:disc%mesh%ncell_onP,0) = Tcell
    call disc%mesh%cell_imap%gather_offp(state(:,0))
    call HT_model%H_of_T%compute_value(state, Hcell)
    deallocate(state)

    !! Compute consistent u and udot
    max_itr = 100
    rel_tol = tol * 1E-5_r8 !TODO: good practice?
    call ic%compute(t, dt, Tcell, u, udot, rel_tol, max_itr, stat, errmsg)
    if (stat/=0) call error_exit(errmsg)

    !! u must match expected values
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

    !! udot must match expected values
    if (global_any(abs(udot%hc(:disc%mesh%ncell_onP)-dcoef) > tol)) then
      if (is_IOP) print '("ERROR: Hdot exceeds expected value; tol=",es9.2)', tol
      status = 1
    end if
    if (global_any(abs(udot%tc(:disc%mesh%ncell_onP)-dcoef) > tol)) then
      if (is_IOP) print '("ERROR: cell temp derivative exceeds expected value; tol=",es9.2)', tol
      status = 1
    end if
    if (global_any(abs(udot%tf(:disc%mesh%nface_onP)-dcoef) > tol)) then
      if (is_IOP) print '("ERROR: face temp derivative exceeds expected value; tol=",es9.2)', tol
      status = 1
    end if

  end subroutine test_linear_flux


  !! Diagnose the discrete residual of the exact one-dimensional two-material
  !! steady conduction solution.  The conductivity interface is aligned with
  !! a mesh face at x=0.5; density and specific heat are one in both regions.
  subroutine diagnose_multimaterial_dirichlet(disc, tol)

    type(mfd_2d_disc), target, intent(in) :: disc
    real(r8), intent(in) :: tol

    type(material_database) :: matl_db
    type(material_model), target :: matl_model
    type(material_composition), target :: composition
    type(ht_2d_model), target :: model
    type(ht_2d_ic_solver) :: ic
    type(parameter_list), pointer :: matl_params, region_params, model_params
    type(ht_2d_vector) :: u, udot, zero_udot, residual
    real(r8), allocatable :: temp(:), state(:,:), conductivity(:)
    real(r8) :: flux, cell_norm, face_norm, hdot_norm, t, dt, rel_tol
    character(:), allocatable :: errmsg, string
    integer :: j, max_itr, stat

    if (is_IOP) print '(/,"Diagnosing two-material Dirichlet problem")'

    string = '{"low":{"properties":{"conductivity":1.0,"density":1.0,"specific-heat":1.0}},&
              &"high":{"properties":{"conductivity":10.0,"density":1.0,"specific-heat":1.0}}}'
    call parameter_list_from_json_string(string, matl_params, errmsg)
    if (.not.associated(matl_params)) call error_exit(errmsg)
    call load_material_database(matl_db, matl_params, stat, errmsg)
    if (stat /= 0) call error_exit(errmsg)
    call matl_model%init(['low ', 'high'], matl_db, stat, errmsg)
    if (stat /= 0) call error_exit(errmsg)
    call add_enthalpy_prop(matl_model, stat, errmsg)
    if (stat /= 0) call error_exit(errmsg)

    string = '{"low-half":{"material":"low","type":"half-plane","point":[0.5,0.0],"normal":[1.0,0.0]},&
              &"high-background":{"material":"high","type":"background"}}'
    call parameter_list_from_json_string(string, region_params, errmsg)
    if (.not.associated(region_params)) call error_exit(errmsg)
    call composition%init(disc%mesh, matl_model, region_params, 6, stat, errmsg)
    if (stat /= 0) call error_exit(errmsg)

    string = '{"bc":{"left":{"type":"temperature","face-set-ids":[1],"temp":1.0},&
              &"right":{"type":"temperature","face-set-ids":[2],"temp":0.0},&
              &"bottom-top":{"type":"flux","face-set-ids":[3,4],"flux":0.0}}}'
    call parameter_list_from_json_string(string, model_params, errmsg)
    if (.not.associated(model_params)) call error_exit(errmsg)
    call model%init(disc, matl_model, composition, model_params, stat, errmsg)
    if (stat /= 0) call error_exit(errmsg)

    flux = 1.0_r8 / (0.5_r8 + 0.5_r8 / 10.0_r8)
    call disc%mesh%init_cell_centroid
    allocate(temp(disc%mesh%ncell_onP))
    do j = 1, size(temp)
      if (disc%mesh%cell_centroid(1,j) <= 0.5_r8) then
        temp(j) = 1.0_r8 - flux*disc%mesh%cell_centroid(1,j)
      else
        temp(j) = 1.0_r8 - 0.5_r8*flux - flux/10.0_r8 * (disc%mesh%cell_centroid(1,j)-0.5_r8)
      end if
    end do
    allocate(state(disc%mesh%ncell_onP,0:0), conductivity(disc%mesh%ncell_onP))
    state(:,0) = temp
    call model%conductivity%compute_value(state, conductivity)

    t = 0.0_r8
    dt = 1.0e-3_r8
    rel_tol = tol*1.0e-5_r8
    max_itr = 100
    call ic%init(model)
    call u%init(disc%mesh)
    call udot%init(u)
    call ic%compute(t, dt, temp, u, udot, rel_tol, max_itr, stat, errmsg)
    if (stat /= 0) call error_exit(errmsg)

    call zero_udot%init(u)
    call zero_udot%setval(0.0_r8)
    call residual%init(u)
    call residual%setval(0.0_r8)
    call model%residual(t, u, zero_udot, residual)
    cell_norm = sqrt(global_dot_product(residual%tc(:disc%mesh%ncell_onP), &
                                        residual%tc(:disc%mesh%ncell_onP)))
    face_norm = sqrt(global_dot_product(residual%tf(:disc%mesh%nface_onP), &
                                        residual%tf(:disc%mesh%nface_onP)))
    hdot_norm = sqrt(global_dot_product(udot%hc(:disc%mesh%ncell_onP), &
                                        udot%hc(:disc%mesh%ncell_onP)))
    if (is_IOP) then
      print '("  conductivity left / right       = ",2(es12.4,1x))', &
        minval(pack(conductivity, disc%mesh%cell_centroid(1,:disc%mesh%ncell_onP) < 0.5_r8)), &
        maxval(pack(conductivity, disc%mesh%cell_centroid(1,:disc%mesh%ncell_onP) > 0.5_r8))
      print '("  ||Rcell(T,Tface; Hdot=0)||_2 = ",es12.4)', cell_norm
      print '("  ||Rface(T,Tface; Hdot=0)||_2 = ",es12.4)', face_norm
      print '("  ||initial Hdot||_2             = ",es12.4)', hdot_norm
    end if
    if (cell_norm > tol .or. face_norm > tol .or. hdot_norm > 10.0_r8*tol) then
      if (is_IOP) print '("ERROR: two-material steady state is not consistent; tol=",es9.2)', tol
      status = 1
    end if

  end subroutine diagnose_multimaterial_dirichlet

end program test_HT_2d_solver_initial_state
