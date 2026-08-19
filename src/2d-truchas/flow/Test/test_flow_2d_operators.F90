module test_flow_2d_bc_types

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use bndry_func1_class
  use bndry_vfunc_class
  implicit none

  type, extends(bndry_func1) :: test_bndry_func1
  contains
    procedure :: compute => compute_func1
  end type

  type, extends(bndry_vfunc) :: test_bndry_vfunc
  contains
    procedure :: compute => compute_vfunc
  end type

contains

  subroutine compute_func1(this, t)
    class(test_bndry_func1), intent(inout) :: this
    real(r8), intent(in) :: t
  end subroutine


  subroutine compute_vfunc(this, t)
    class(test_bndry_vfunc), intent(inout) :: this
    real(r8), intent(in) :: t
  end subroutine

end module test_flow_2d_bc_types


program test_flow_2d_operators

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication
  use truchas_env, only: prefix, overwrite_output
  use truchas_logging_services
  use unstr_2d_mesh_type
  use unstr_2d_mesh_factory
  use parameter_list_type
  use flow_2d_state_type
  use flow_2d_operators_type
  use test_flow_2d_bc_types
  implicit none

  integer :: status

  call init_parallel_communication
  prefix = 'run'
  overwrite_output = .true.
  call TLS_initialize
  call TLS_set_verbosity(TLS_VERB_NORMAL)

  status = 0
  if (command_argument_count() == 0) then
    call test_mesh(0.0_r8, 'quad mesh')
    call test_mesh(0.0_r8, 'rotated quad mesh', 30.0_r8)
    call test_mesh(0.5_r8, 'mixed triangle/quad mesh')
  else
    call test_external_mesh()
  end if

  call halt_parallel_communication
  stop status

contains

  subroutine test_mesh(triangle_probability, name, rotation_angle)
    real(r8), intent(in) :: triangle_probability
    character(*), intent(in) :: name
    real(r8), optional, intent(in) :: rotation_angle

    type(unstr_2d_mesh), pointer :: mesh
    type(flow_2d_state) :: state
    type(flow_2d_operators) :: ops

    mesh => new_unstr_2d_mesh([0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [8, 8], &
        0.0_r8, triangle_probability)
    if (present(rotation_angle)) call rotate_mesh(mesh, rotation_angle)
    call ops%init(mesh)
    call state%init(mesh)
    call test_divergence(mesh, ops, state, name)
    call test_derivative(mesh, ops, state, name, triangle_probability == 0.0_r8)
    call test_interpolation(mesh, ops, state, name)
    call test_boundary_operators(mesh, ops, state, name, triangle_probability == 0.0_r8)
  end subroutine


  subroutine test_external_mesh()
    type(parameter_list) :: params
    type(unstr_2d_mesh), pointer :: mesh
    type(flow_2d_state) :: state
    type(flow_2d_operators) :: ops
    character(512) :: path
    character(:), allocatable :: errmsg
    integer :: stat

    call get_command_argument(1, path)
    call params%set('mesh-file', trim(path))
    mesh => new_unstr_2d_mesh(params, stat, errmsg)
    if (stat /= 0) return
    call require(associated(mesh), 'external QUAD4 mesh initialization failed')
    if (.not.associated(mesh)) return

    call ops%init(mesh)
    call state%init(mesh)
    call test_divergence(mesh, ops, state, 'external quad mesh')
    call test_derivative(mesh, ops, state, 'external quad mesh', .false.)
    call test_interpolation(mesh, ops, state, 'external quad mesh')
    call test_boundary_operators(mesh, ops, state, 'external quad mesh', .false.)
    call report_derivative_accuracy(mesh, ops, state)
  end subroutine


  subroutine test_boundary_operators(mesh, ops, state, name, check_exact)
    type(unstr_2d_mesh), target, intent(in) :: mesh
    type(flow_2d_operators), intent(in) :: ops
    type(flow_2d_state), intent(inout) :: state
    character(*), intent(in) :: name
    logical, intent(in) :: check_exact

    type(test_bndry_func1) :: pdir, pneu, zero_normal
    type(test_bndry_vfunc) :: vdir
    real(r8), allocatable :: derivative(:), velocity_f(:), gradient_c(:,:)
    real(r8), parameter :: grad(2) = [1.0_r8, 2.0_r8]
    real(r8), parameter :: velocity(2) = [1.5_r8, -0.75_r8]
    integer :: c, f, i, nboundary
    logical :: dirichlet_ok, neumann_ok, velocity_ok, gradient_ok

    do c = 1, mesh%ncell
      state%p_cc(c) = dot_product(grad, mesh%cell_centroid(:,c))
      state%vel_cc(:,c) = velocity
    end do
    nboundary = count(mesh%fcell(2,1:mesh%nface_onP) == 0)
    allocate(pdir%index(nboundary), pdir%value(nboundary), pneu%index(nboundary), pneu%value(nboundary), &
        zero_normal%index(nboundary), zero_normal%value(nboundary), vdir%index(nboundary), &
        vdir%value(2,nboundary), derivative(mesh%nface), velocity_f(mesh%nface), gradient_c(2,mesh%ncell))
    i = 0
    do f = 1, mesh%nface_onP
      if (mesh%fcell(2,f) /= 0) cycle
      i = i + 1
      pdir%index(i) = f
      pdir%value(i) = dot_product(grad, mesh%face_centroid(:,f))
      pneu%index(i) = f
      pneu%value(i) = dot_product(grad, mesh%unit_normal(:,f))
      zero_normal%index(i) = f
      zero_normal%value(i) = 0.0_r8
      vdir%index(i) = f
      vdir%value(:,i) = velocity
    end do

    call ops%derivative_cf(state%p_cc, derivative, dirichlet_bc=pdir)
    dirichlet_ok = .true.
    do i = 1, nboundary
      f = pdir%index(i)
      if (check_exact) dirichlet_ok = dirichlet_ok .and. &
          abs(derivative(f) - dot_product(grad, mesh%unit_normal(:,f))) < 1.0e-12_r8
    end do
    call require(dirichlet_ok, name // ': pressure Dirichlet derivative is incorrect')

    call ops%derivative_cf(state%p_cc, derivative, normal_flux_bc=pneu)
    neumann_ok = .true.
    do i = 1, nboundary
      f = pneu%index(i)
      neumann_ok = neumann_ok .and. abs(derivative(f) - pneu%value(i)) < 1.0e-12_r8
    end do
    call require(neumann_ok, name // ': pressure Neumann derivative is incorrect')

    call ops%gradient_cc(state%p_cc, gradient_c, dirichlet_bc=pdir)
    gradient_ok = .true.
    do c = 1, mesh%ncell_onP
      gradient_ok = gradient_ok .and. maxval(abs(gradient_c(:,c) - grad)) < 1.0e-12_r8
    end do
    call require(gradient_ok, name // ': pressure Dirichlet gradient is incorrect')

    call ops%gradient_cc(state%p_cc, gradient_c, normal_flux_bc=pneu)
    gradient_ok = .true.
    do c = 1, mesh%ncell_onP
      gradient_ok = gradient_ok .and. maxval(abs(gradient_c(:,c) - grad)) < 1.0e-12_r8
    end do
    call require(gradient_ok, name // ': pressure Neumann gradient is incorrect')

    call ops%interpolate_cf(state%vel_cc, velocity_f, dirichlet_bc=vdir)
    velocity_ok = .true.
    do i = 1, nboundary
      f = vdir%index(i)
      velocity_ok = velocity_ok .and. abs(velocity_f(f) - &
          dot_product(mesh%unit_normal(:,f), velocity)) < 1.0e-12_r8
    end do
    call require(velocity_ok, name // ': velocity Dirichlet interpolation is incorrect')

    call ops%interpolate_cf(state%vel_cc, velocity_f, zero_normal_bc=zero_normal)
    velocity_ok = .true.
    do i = 1, nboundary
      velocity_ok = velocity_ok .and. abs(velocity_f(zero_normal%index(i))) < 1.0e-12_r8
    end do
    call require(velocity_ok, name // ': zero-normal velocity interpolation is incorrect')
  end subroutine


  subroutine report_derivative_accuracy(mesh, ops, state)
    type(unstr_2d_mesh), target, intent(in) :: mesh
    type(flow_2d_operators), intent(in) :: ops
    type(flow_2d_state), intent(inout) :: state

    real(r8), allocatable :: derivative(:)
    real(r8), parameter :: exact_grad(2) = [1.0_r8, 2.0_r8]
    real(r8) :: exact, error
    integer :: c, f

    do c = 1, mesh%ncell
      state%p_cc(c) = dot_product(exact_grad, mesh%cell_centroid(:,c))
    end do
    allocate(derivative(mesh%nface))
    call ops%derivative_cf(state%p_cc, derivative)
    error = 0.0_r8
    do f = 1, mesh%nface_onP
      if (mesh%fcell(2,f) == 0) cycle
      exact = dot_product(exact_grad, mesh%unit_normal(:,f))
      error = max(error, abs(derivative(f) - exact))
    end do
    if (is_IOP) print '("external quad mesh: max linear face derivative error = ",es12.5)', error
  end subroutine


  subroutine rotate_mesh(mesh, angle)
    type(unstr_2d_mesh), intent(inout) :: mesh
    real(r8), intent(in) :: angle

    real(r8) :: theta, rotation(2,2)

    theta = angle*acos(-1.0_r8)/180.0_r8
    rotation = reshape([cos(theta), sin(theta), -sin(theta), cos(theta)], [2,2])
    mesh%x = matmul(rotation, mesh%x)
    call mesh%compute_geometry()
  end subroutine


  subroutine test_divergence(mesh, ops, state, name)
    type(unstr_2d_mesh), target, intent(in) :: mesh
    type(flow_2d_operators), intent(in) :: ops
    type(flow_2d_state), intent(inout) :: state
    character(*), intent(in) :: name

    real(r8), allocatable :: div(:)
    real(r8), parameter :: vel(2) = [1.5_r8, -0.75_r8]
    integer :: f

    do f = 1, mesh%nface
      state%vel_fn(f) = dot_product(vel, mesh%unit_normal(:,f))
    end do
    allocate(div(mesh%ncell_onP))
    call ops%divergence(state%vel_fn, div)
    call require(maxval(abs(div)) < 1.0e-12_r8, name // ': constant velocity divergence is nonzero')
  end subroutine


  subroutine test_derivative(mesh, ops, state, name, check_exact)
    type(unstr_2d_mesh), target, intent(in) :: mesh
    type(flow_2d_operators), intent(in) :: ops
    type(flow_2d_state), intent(inout) :: state
    character(*), intent(in) :: name
    logical, intent(in) :: check_exact

    real(r8), allocatable :: derivative(:)
    integer :: c, f, nface
    logical :: derivative_ok, exact_ok
    real(r8) :: exact_grad(2), expected, dx

    do c = 1, mesh%ncell
      state%p_cc(c) = mesh%cell_centroid(1,c) + 2.0_r8*mesh%cell_centroid(2,c)
    end do
    allocate(derivative(mesh%nface))
    call ops%derivative_cf(state%p_cc, derivative)
    exact_grad = [1.0_r8, 2.0_r8]
    nface = 0
    derivative_ok = .true.
    exact_ok = .true.
    do f = 1, mesh%nface_onP
      if (mesh%fcell(2,f) > 0) then
        nface = nface + 1
        dx = 2.0_r8*min( &
            dot_product(mesh%face_centroid(:,f) - mesh%cell_centroid(:,mesh%fcell(1,f)), mesh%unit_normal(:,f)), &
            -dot_product(mesh%face_centroid(:,f) - mesh%cell_centroid(:,mesh%fcell(2,f)), mesh%unit_normal(:,f)))
        expected = dot_product(exact_grad, mesh%cell_centroid(:,mesh%fcell(2,f)) - &
            mesh%cell_centroid(:,mesh%fcell(1,f)))/dx
        derivative_ok = derivative_ok .and. abs(derivative(f) - expected) < 1.0e-12_r8
        if (check_exact) exact_ok = exact_ok .and. &
            abs(derivative(f) - dot_product(exact_grad, mesh%unit_normal(:,f))) < 1.0e-12_r8
      end if
    end do
    call require(derivative_ok, name // ': linear pressure derivative is incorrect')
    if (check_exact) call require(exact_ok, name // ': linear pressure derivative is not rotation invariant')
    call require(nface > 0, name // ': no interior faces were tested')
  end subroutine


  subroutine test_interpolation(mesh, ops, state, name)
    type(unstr_2d_mesh), target, intent(in) :: mesh
    type(flow_2d_operators), intent(in) :: ops
    type(flow_2d_state), intent(inout) :: state
    character(*), intent(in) :: name

    real(r8), allocatable :: scalar_f(:)
    real(r8), parameter :: velocity(2) = [1.5_r8, -0.75_r8]
    real(r8) :: dcc(2), factor, expected
    integer :: c, f, c1, c2
    logical :: scalar_ok, velocity_ok

    do c = 1, mesh%ncell
      state%p_cc(c) = mesh%cell_centroid(1,c) + 2.0_r8*mesh%cell_centroid(2,c)
      state%vel_cc(:,c) = velocity
    end do
    allocate(scalar_f(mesh%nface))
    call ops%interpolate_cf(state%p_cc, scalar_f)
    call ops%interpolate_cf(state%vel_cc, state%vel_fn)
    scalar_ok = .true.
    velocity_ok = .true.
    do f = 1, mesh%nface_onP
      c1 = mesh%fcell(1,f)
      c2 = mesh%fcell(2,f)
      factor = 1.0_r8
      if (c2 > 0) then
        dcc = mesh%cell_centroid(:,c2) - mesh%cell_centroid(:,c1)
        factor = 1.0_r8 - dot_product(mesh%face_centroid(:,f) - mesh%cell_centroid(:,c1), dcc)/dot_product(dcc, dcc)
      end if
      expected = factor*state%p_cc(c1)
      if (c2 > 0) expected = expected + (1.0_r8-factor)*state%p_cc(c2)
      scalar_ok = scalar_ok .and. abs(scalar_f(f) - expected) < 1.0e-12_r8
      velocity_ok = velocity_ok .and. abs(state%vel_fn(f) - dot_product(velocity, mesh%unit_normal(:,f))) < 1.0e-12_r8
    end do
    call require(scalar_ok, name // ': scalar interpolation is incorrect')
    call require(velocity_ok, name // ': velocity interpolation is incorrect')
  end subroutine


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (global_any(.not.condition)) then
      if (is_IOP) print '("ERROR: ",a)', message
      status = 1
    end if
  end subroutine

end program test_flow_2d_operators
