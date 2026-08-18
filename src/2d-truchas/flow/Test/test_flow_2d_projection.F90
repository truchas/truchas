program test_flow_2d_projection

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication
  use truchas_env, only: prefix, overwrite_output
  use truchas_logging_services
  use unstr_2d_mesh_type
  use unstr_2d_mesh_factory
  use flow_2d_operators_type
  use flow_2d_projection_type
  use pcsr_matrix_type
  implicit none

  integer :: status

  call init_parallel_communication
  prefix = 'run'
  overwrite_output = .true.
  call TLS_initialize
  call TLS_set_verbosity(TLS_VERB_NORMAL)

  status = 0
  call test_projection

  call halt_parallel_communication
  stop status

contains

  subroutine test_projection
    type(unstr_2d_mesh), pointer :: mesh
    type(flow_2d_operators), target :: operators
    type(flow_2d_projection), target :: projection
    type(pcsr_matrix), pointer :: matrix
    real(r8), allocatable :: inv_density_f(:), p(:), one(:), result(:)
    real(r8), parameter :: gradient(2) = [1.0_r8, 2.0_r8]
    integer :: c, f, pin_face
    logical :: constant_nullspace, linear_harmonic, pinned

    mesh => new_unstr_2d_mesh([0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8], [8, 8], 0.0_r8, 0.0_r8)
    call rotate_mesh(mesh, 30.0_r8)
    call operators%init(mesh)
    call projection%init(mesh, operators)
    allocate(inv_density_f(mesh%nface), p(mesh%ncell), one(mesh%ncell), result(mesh%ncell_onP))
    inv_density_f = 1.0_r8
    one = 1.0_r8
    do c = 1, mesh%ncell
      p(c) = dot_product(gradient, mesh%cell_centroid(:,c))
    end do

    call projection%assemble(inv_density_f)
    matrix => projection%matrix()
    call matrix%matvec(one, result)
    constant_nullspace = maxval(abs(result)) < 1.0e-12_r8
    call require(constant_nullspace, 'unreferenced pressure does not have a constant nullspace')

    call matrix%matvec(p, result)
    linear_harmonic = .true.
    do c = 1, mesh%ncell_onP
      if (all(mesh%cnhbr(mesh%cstart(c):mesh%cstart(c+1)-1) > 0)) then
        linear_harmonic = linear_harmonic .and. abs(result(c)) < 1.0e-12_r8
      end if
    end do
    call require(linear_harmonic, 'projection operator is not rotation invariant on an orthogonal mesh')

    pin_face = 0
    if (this_PE == 1) then
      do f = 1, mesh%nface_onP
        if (mesh%fcell(2,f) == 0) then
          pin_face = f
          exit
        end if
      end do
    end if
    call projection%assemble(inv_density_f, pin_face)
    call matrix%matvec(one, result)
    pinned = maxval(abs(result)) > 1.0e-12_r8
    call require(global_any(pinned), 'pressure pin did not remove the constant nullspace')
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


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (global_any(.not.condition)) then
      if (is_IOP) print '("ERROR: ",a)', message
      status = 1
    end if
  end subroutine

end program test_flow_2d_projection
