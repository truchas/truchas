#include "f90_assert.fpp"

module fdme_minres_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use complex_lin_op_class
  use fdme_model_type
  use fdme_precon_class
  use vector_class
  use cs_minres_solver_type
  implicit none
  private

  type, extends(complex_lin_op) :: fdme_lin_op
    type(fdme_model), pointer :: model => null() ! unowned reference
    class(fdme_precon), pointer :: my_precon => null() ! unowned reference
    real(r8), allocatable :: dinv(:)
  contains
    procedure :: matvec
    procedure :: precon
  end type

  type, public :: fdme_minres_solver
    type(fdme_model), pointer :: model => null() ! unowned reference
    type(cs_minres_solver) :: minres
    type(fdme_lin_op) :: lin_op
  contains
    procedure :: init
    procedure :: solve
  end type

contains

  subroutine matvec(this, x, y)
    class(fdme_lin_op), intent(inout) :: this
    complex(r8) :: x(:), y(:)
    call this%model%mesh%edge_imap%gather_offp(x)
    call this%model%A%matvec(x, y)
  end subroutine

  subroutine precon(this, x, y)
    class(fdme_lin_op), intent(inout) :: this
    complex(r8) :: x(:), y(:)
    integer :: nedge_onP
    if (.not.allocated(this%dinv)) then
      block
        integer :: j
        associate (A => this%model%A)
          call A%kdiag_init
          allocate(this%dinv(A%nrow))
          do j = 1, A%nrow
            this%dinv(j) = 1.0_r8 / A%values(A%kdiag(j))%re
          end do
        end associate
      end block
    end if
    !y = x ! no preconditioning
    !y = this%dinv*x ! diagonal preconditioning
    !block ! doesn't work with hiptmair (doesn't satisfy requirements)
    !  real(r8) :: xarray(2,size(x))
    !  xarray(1,:) = x%re; xarray(2,:) = x%im
    !  call this%my_precon%apply(xarray)
    !  y%re = xarray(1,:); y%im = xarray(2,:)
    !end block
    nedge_onP = this%model%mesh%nedge_onP
    call this%model%mesh%edge_imap%gather_offp(x)
    y = 0.0_r8
    call gs_relaxation(this%model%A, x(:nedge_onP), y, 'fb')
    call this%model%mesh%edge_imap%scatter_offp_sum(y)
    call this%model%mesh%edge_imap%gather_offp(y)
  end subroutine

  subroutine gs_relaxation(A, f, u, pattern)

    use complex_pcsr_matrix_type

    type(complex_pcsr_matrix), intent(inout) :: A
    complex(r8), intent(in) :: f(:)
    complex(r8), intent(inout) :: u(:)
    character(*), intent(in) :: pattern

    integer :: i, i1, i2, di, j, k, n
    complex(r8) :: s

    ASSERT(A%nrow == A%ncol)
    ASSERT(size(u) >= A%ncol)

    n = min(A%nrow, size(f))

    if (.not.allocated(A%kdiag)) call A%kdiag_init

    do j = 1, len(pattern)
      call loop_range(pattern(j:j), n, i1, i2, di)
      do i = i1, i2, di
        s = f(i)
        do k = A%graph%xadj(i), A%graph%xadj(i+1)-1
          s = s - A%values(k)%re * u(A%graph%adjncy(k))
        end do
        u(i) = u(i) + s / A%values(A%kdiag(i))%re
      end do
    end do

  end subroutine gs_relaxation

  subroutine loop_range(direction, len, i1, i2, di)
    character(1), intent(in) :: direction
    integer, intent(in) :: len
    integer, intent(out) :: i1, i2, di
    select case (direction)
    case ('f', 'F') ! forward sweep
      i1 = 1
      i2 = len
      di = 1
    case ('b', 'B') ! backward sweep
      i1 = len
      i2 = 1
      di = -1
    case default
      INSIST(.false.)
    end select

  end subroutine

  subroutine init(this, model, precon, params, stat, errmsg)
    use parameter_list_type
    class(fdme_minres_solver), intent(out) :: this
    type(fdme_model), pointer :: model !TODO: don't make a pointer
    class(fdme_precon), intent(in), target :: precon
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    this%model => model
    call this%minres%init(model%mesh%nedge_onP, params)
    this%lin_op%model => model
    this%lin_op%my_precon => precon
    stat = 0
  end subroutine

  subroutine solve(this, efield, stat)
    class(fdme_minres_solver), intent(inout) :: this
    complex(r8), intent(inout) :: efield(:)
    integer, intent(out) :: stat
    call this%minres%solve(this%lin_op, this%model%rhs, efield)
    stat = 0 !FIXME: need to extract from minres
  end subroutine

end module fdme_minres_solver_type
