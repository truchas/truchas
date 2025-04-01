#include "f90_assert.fpp"

module thes_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use complex_lin_op_class
  use simpl_mesh_type
  use complex_pcsr_matrix_type
  use pcsr_matrix_type
  use pcsr_precon_class
  use cs_minres_solver_type
  use parameter_list_type
  implicit none
  private

  type, extends(complex_lin_op), public :: thes_solver
    type(simpl_mesh), pointer :: mesh => null() ! unowned reference
    type(complex_pcsr_matrix) :: A
    complex(r8), allocatable :: rhs(:)
    type(pcsr_matrix), pointer :: M => null() ! pointer to avoid dangling pointer
    class(pcsr_precon), allocatable :: pc
    type(cs_minres_solver) :: minres
  contains
    procedure :: init
    procedure :: solve
    ! deferred procedures from complex_lin_op class
    procedure :: matvec
    procedure :: precon
  end type

contains

  subroutine matvec(this, x, y)
    class(thes_solver), intent(inout) :: this
    complex(r8) :: x(:), y(:)
    call this%mesh%node_imap%gather_offp(x)
    call this%A%matvec(x, y)
  end subroutine

  subroutine precon(this, x, y)
    class(thes_solver), intent(inout) :: this
    complex(r8) :: x(:), y(:)
    call this%mesh%node_imap%gather_offp(x)
    y = x
    call this%pc%apply(y%re)
    call this%pc%apply(y%im)
    call this%mesh%node_imap%gather_offp(y) ! necessary?
  end subroutine

  subroutine init(this, mesh, eps, bc, params, stat, errmsg)

    use thes_bc_type

    class(thes_solver), intent(out) :: this
    type(simpl_mesh), intent(in), target :: mesh
    complex(r8), intent(in) :: eps(:)
    type(thes_bc), intent(in) :: bc
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: j

    ASSERT(size(eps) == mesh%ncell)

    this%mesh => mesh

    !! Raw system matrix ignoring boundary conditions
    block
      use mimetic_discretization, only: W1_matrix_WE, cell_grad
      use upper_packed_matrix_procs, only: upm_cong_prod
      type(pcsr_graph), pointer :: g
      real(r8) :: m1(21), gtm1g(10)
      allocate(g)
      call g%init(this%mesh%node_imap)
      call g%add_clique(this%mesh%cnode)
      call g%add_complete
      call this%A%init(g, take_graph=.true.)
      do j = 1, this%mesh%ncell
        m1 = W1_matrix_WE(this%mesh, j)
        gtm1g = upm_cong_prod(6, 4, m1, cell_grad)
        call this%A%add_to(this%mesh%cnode(:,j), eps(j)*gtm1g)
      end do
    end block

    !! RHS vector
    allocate(this%rhs(this%mesh%nnode))
    this%rhs = 0.0_r8

    !TODO: RHS contribution from Dirichlet conditions
    if (allocated(bc%dirichlet)) then
      block
        !complex(r8) :: r(this%mesh%nnode)
        complex(r8) :: r(mesh%nnode)
        associate (index => bc%dirichlet%index, value => bc%dirichlet%value)
          do j = 1, size(index)
            this%rhs(index(j)) = value(j)
          end do
          call mesh%node_imap%gather_offp(this%rhs)
          call this%A%matvec(this%rhs, r)
          do j = 1, size(index)
            r(index(j)) = 0.0_r8
          end do
        end associate
        this%rhs(:this%mesh%nnode_onp) = this%rhs(:this%mesh%nnode_onp) - r(:this%mesh%nnode_onp)
        call mesh%node_imap%gather_offp(this%rhs)
      end block
    end if

    !TODO: modify system matrix for Dirichlet conditions
    if (allocated(bc%dirichlet)) then
      associate (index => bc%dirichlet%index)
        do j = 1, size(index)
          call this%A%project_out(index(j))
          call this%A%set(index(j), index(j), cmplx(1,0,kind=r8))
        end do
      end associate
    end if

    !! Initialize the preconditioner; based on real part of the system matrix
    block
      use pcsr_precon_factory, only: alloc_pcsr_precon
      type(parameter_list), pointer :: plist
      allocate(this%M)
      call this%M%init(this%A%graph, take_graph=.false.)
      this%M%values = this%A%values%re
      plist => params%sublist('preconditioner', stat, errmsg)
      if (stat /= 0) return
      call alloc_pcsr_precon(this%pc, this%M, plist, stat, errmsg)
      if (stat /= 0) return
      call this%pc%compute
    end block

    call this%minres%init(this%mesh%nnode_onp, params)

  end subroutine init

  subroutine solve(this, phi, stat, errmsg)

    use string_utilities, only: i_to_c

    class(thes_solver), intent(inout) :: this
    complex(r8), intent(out) :: phi(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    call this%minres%solve(this, this%rhs, phi)
    select case (this%minres%flag)
    case (0,1,3) ! conventional success cases
      stat = 0
      return
    case (-3) ! preconditioner is not positive definite
      stat = 1
      errmsg = 'CS-MINRES solve failed: preconditioner not positive definite'
      return
    case (2,4) ! "successful" in some sense, but with an incompatible system
      stat = 1
    case (-1,5) ! found an eigenvector instead
      stat = 1
    case (6:8) ! Gave up iteration due to limits
      stat = 1
    case (9) ! system is singular
      stat = 1
    case default
      stat = 1
    end select
    errmsg = 'CS-MINRES solve failed: stat=' // i_to_c(this%minres%flag)
  end subroutine

end module thes_solver_type
