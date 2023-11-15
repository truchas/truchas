!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_HT_2d_norm_type

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
  use material_model_type
  use source_mesh_function
  use scalar_func_factories
  use mfd_2d_disc_type
  use HT_2d_model_type
  use HT_2d_norm_type
  use bitfield_type
  use test_ht_2d_common
  implicit none

  type(unstr_2d_mesh), pointer :: mesh
  type(mfd_2d_disc), target :: mfd_disc
  type(material_model), target :: matl_model
  type(matl_mesh_func), target :: mmf
  real(r8) :: xmin(2), xmax(2), eps
  integer  :: nx(2)
  integer :: status = 0

  !! Initialize MPI and other base stuff that Truchas depends on
  call init_parallel_communication
  prefix='run'  ! TLS will write to 'run.log'
  overwrite_output = .true.
  call TLS_initialize
  call TLS_set_verbosity(TLS_VERB_NORMAL)

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
  call test_abs(mfd_disc, matl_model)
  call test_rel(mfd_disc, matl_model)
  call test_mixed(mfd_disc, matl_model)
  call test_global(mfd_disc, matl_model)

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


  !! Test the HT_2d_norm type with only the absolute tolerance set
  subroutine test_abs(disc, matl_model)

    type(mfd_2d_disc), target, intent(in) :: disc
    type(material_model), target, intent(in) :: matl_model

    type(HT_2d_model), target :: HT_model
    type(HT_2d_norm), target :: HT_norm
    type(parameter_list), pointer :: params, sublist
    real(r8), allocatable :: u(:), du(:)
    character(:), allocatable :: errmsg, string
    real(r8) :: du_norm
    integer :: n, stat

    if (is_IOP) print '(/,"Testing absolute tolerance")'

    !! Define problem parameters
    string = &
    '{"model": { &
        "bc": { &
          "all-sides": { &
            "type": "temperature", &
            "face-set-ids": [1,2,3,4], &
            "temp": 0.0}}}, &
      "error-norm": { &
        "temp-abs-tol": 0.5, &
        "enth-abs-tol": 0.5}}'
    call parameter_list_from_json_string(string, params, errmsg)
    if (.not. associated(params)) call error_exit(errmsg)

    !! Initialize 2D HT model
    sublist => params%sublist('model')
    call HT_model%init(disc, matl_model, sublist, stat, errmsg)
    if (stat /= 0) call error_exit(errmsg)

    !! Initialize 2D HT norm
    sublist => params%sublist('error-norm')
    call HT_norm%init(HT_model, sublist)

    !! Compute norm
    n = HT_model%num_dof()
    allocate(u(n), du(n))
    u = -10.0_r8
    du = -0.5_r8
    call HT_norm%compute(u, du, du_norm)

    !! Check result
    if (global_maxval(abs(du_norm-1.0_r8)) /= 0.0_r8) then
      if (is_IOP) print '("ERROR: norm is nonzero")'
      status = 1
    end if

  end subroutine test_abs


  !! Test the HT_2d_norm type with only the relative tolerance set
  subroutine test_rel(disc, matl_model)

    type(mfd_2d_disc), target, intent(in) :: disc
    type(material_model), target, intent(in) :: matl_model

    type(HT_2d_model), target :: HT_model
    type(HT_2d_norm), target :: HT_norm
    type(parameter_list), pointer :: params, sublist
    real(r8), allocatable :: u(:), du(:)
    character(:), allocatable :: errmsg, string
    real(r8) :: du_norm
    integer :: n, stat

    if (is_IOP) print '(/,"Testing relative tolerance")'

    !! Define problem parameters
    string = &
    '{"model": { &
        "bc": { &
          "all-sides": { &
            "type": "temperature", &
            "face-set-ids": [1,2,3,4], &
            "temp": 0.0}}}, &
      "error-norm": { &
        "temp-rel-tol": 0.25, &
        "enth-rel-tol": 0.25}}'
    call parameter_list_from_json_string(string, params, errmsg)
    if (.not. associated(params)) call error_exit(errmsg)

    !! Initialize 2D HT model
    sublist => params%sublist('model')
    call HT_model%init(disc, matl_model, sublist, stat, errmsg)
    if (stat /= 0) call error_exit(errmsg)

    !! Initialize 2D HT norm
    sublist => params%sublist('error-norm')
    call HT_norm%init(HT_model, sublist)

    !! Compute norm
    n = HT_model%num_dof()
    allocate(u(n), du(n))
    u = -2.0_r8
    du = -1.0_r8
    call HT_norm%compute(u, du, du_norm)

    !! Check result
    if (global_maxval(abs(du_norm-2.0_r8)) /= 0.0_r8) then
      if (is_IOP) print '("ERROR: norm is nonzero")'
      status = 1
    end if

  end subroutine test_rel


  !! Test the HT_2d_norm type with mixed tolerances
  subroutine test_mixed(disc, matl_model)

    type(mfd_2d_disc), target, intent(in) :: disc
    type(material_model), target, intent(in) :: matl_model

    type(HT_2d_model), target :: HT_model
    type(HT_2d_norm), target :: HT_norm
    type(parameter_list), pointer :: params, sublist
    real(r8), allocatable :: u(:), du(:)
    character(:), allocatable :: errmsg, string
    real(r8) :: du_norm
    integer :: n, stat

    if (is_IOP) print '(/,"Testing mixed tolerances")'

    !! Define problem parameters
    string = &
    '{"model": { &
        "bc": { &
          "all-sides": { &
            "type": "temperature", &
            "face-set-ids": [1,2,3,4], &
            "temp": 0.0}}}, &
      "error-norm": { &
        "temp-abs-tol": 1.0, &
        "temp-rel-tol": 2.0, &
        "enth-abs-tol": 1.0, &
        "enth-rel-tol": 2.0}}'
    call parameter_list_from_json_string(string, params, errmsg)
    if (.not. associated(params)) call error_exit(errmsg)

    !! Initialize 2D HT model
    sublist => params%sublist('model')
    call HT_model%init(disc, matl_model, sublist, stat, errmsg)
    if (stat /= 0) call error_exit(errmsg)

    !! Initialize 2D HT norm
    sublist => params%sublist('error-norm')
    call HT_norm%init(HT_model, sublist)

    !! Compute norm
    n = HT_model%num_dof()
    allocate(u(n), du(n))
    u = -1.5_r8
    du = -4.0_r8
    call HT_norm%compute(u, du, du_norm)

    !! Check result
    if (global_maxval(abs(du_norm-1.0_r8)) /= 0.0_r8) then
      if (is_IOP) print '("ERROR: norm is nonzero")'
      status = 1
    end if

  end subroutine test_mixed


  !! Test the global consistency of the HT_2d_norm type
  subroutine test_global(disc, matl_model)

    type(mfd_2d_disc), target, intent(in) :: disc
    type(material_model), target, intent(in) :: matl_model

    type(HT_2d_model), target :: HT_model
    type(HT_2d_norm), target :: HT_norm
    type(parameter_list), pointer :: params, sublist
    real(r8), allocatable :: u(:), du(:)
    character(:), allocatable :: errmsg, string
    real(r8) :: du_norm, maxnorm, minnorm
    integer :: n, stat

    if (is_IOP) print '(/,"Testing global consistency")'

    !! Define problem parameters
    string = &
    '{"model": { &
        "bc": { &
          "all-sides": { &
            "type": "temperature", &
            "face-set-ids": [1,2,3,4], &
            "temp": 0.0}}}, &
      "error-norm": { &
        "temp-abs-tol": 1.0, &
        "enth-abs-tol": 1.0}}'
    call parameter_list_from_json_string(string, params, errmsg)
    if (.not. associated(params)) call error_exit(errmsg)

    !! Initialize 2D HT model
    sublist => params%sublist('model')
    call HT_model%init(disc, matl_model, sublist, stat, errmsg)
    if (stat /= 0) call error_exit(errmsg)

    !! Initialize 2D HT norm
    sublist => params%sublist('error-norm')
    call HT_norm%init(HT_model, sublist)

    !! Compute norm
    n = HT_model%num_dof()
    allocate(u(n), du(n))
    u = 0.0_r8
    du = this_PE
    call HT_norm%compute(u, du, du_norm)

    !! Check result
    maxnorm = global_maxval(du_norm)
    minnorm = global_minval(du_norm)
    if (maxnorm /= minnorm) then
      if (is_IOP) print '("ERROR: max and min norms do not match")'
      status = 1
    end if
    if (du_norm /= nPE) then
      if (is_IOP) print '("ERROR: du_norm is not the maximum rank; du_norm = ",es10.3)', du_norm
      status = 1
    end if

  end subroutine test_global

end program test_HT_2d_norm_type
