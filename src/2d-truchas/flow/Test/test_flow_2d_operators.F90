program test_flow_2d_operators

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication
  use truchas_env, only: prefix, overwrite_output
  use truchas_logging_services
  use unstr_2d_mesh_type
  use unstr_2d_mesh_factory
  use flow_2d_state_type
  use flow_2d_operators_type
  implicit none

  integer :: status

  call init_parallel_communication
  prefix = 'run'
  overwrite_output = .true.
  call TLS_initialize
  call TLS_set_verbosity(TLS_VERB_NORMAL)

  status = 0
  call test_mesh(0.0_r8, 'quad mesh')
  call test_mesh(0.0_r8, 'rotated quad mesh', 30.0_r8)
  call test_mesh(0.5_r8, 'mixed triangle/quad mesh')

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
