#include "f90_assert.fpp"

module fdme_precon_gs_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use fdme_precon_class
  use fdme_model_type
  implicit none
  private

  type, extends(fdme_precon), public :: fdme_precon_gs
    real(r8), allocatable :: dinv(:)
    integer :: num_iter
  contains
    procedure :: init
    procedure :: setup
    procedure :: apply
  end type

contains

  subroutine init(this, model, params, stat, errmsg)
    use parameter_list_type
    class(fdme_precon_gs), intent(out) :: this
    type(fdme_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    character(:), allocatable :: context
    this%model => model
    stat = 0
    call params%get('num-cycles', this%num_iter, stat, errmsg, default=1)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    if (this%num_iter <= 0) then
      stat = 1
      errmsg = context // '"num-cycles" must be > 0'
      return
    end if
  end subroutine

  subroutine setup(this)
    class(fdme_precon_gs), intent(inout) :: this
    integer :: j
    associate (A => this%model%A)
      call A%kdiag_init
      allocate(this%dinv(A%nrow))
      do j = 1, A%nrow
        this%dinv(j) = 1.0_r8 / A%values(A%kdiag(j))%re
      end do
    end associate
  end subroutine

  subroutine apply(this, x, y)
    class(fdme_precon_gs), intent(inout) :: this
    complex(r8), intent(in)  :: x(:)
    complex(r8), intent(out) :: y(:)
    integer :: nedge_onP
    nedge_onP = this%model%mesh%nedge_onP
    y = 0
    call gs_relaxation(this%model%A, x(:nedge_onP), y, 'fb', this%num_iter)
  end subroutine


  subroutine gs_relaxation(A, f, u, pattern, num_iter)

    use complex_pcsr_matrix_type

    type(complex_pcsr_matrix), intent(inout) :: A
    complex(r8), intent(in) :: f(:)
    complex(r8), intent(inout) :: u(:)
    character(*), intent(in) :: pattern
    integer, intent(in) :: num_iter

    integer :: i, i1, i2, di, j, k, n, iter
    complex(r8) :: s

    ASSERT(A%nrow == A%ncol)
    ASSERT(size(u) >= A%ncol)

    n = min(A%nrow, size(f))

    if (.not.allocated(A%kdiag)) call A%kdiag_init

    do iter = 1, num_iter
      do j = 1, len(pattern)
        call loop_range(pattern(j:j), n, i1, i2, di)
        do i = i1, i2, di
          s = f(i)
          do k = A%graph%xadj(i), A%graph%xadj(i+1)-1
            s = s - A%values(k)%re * u(A%graph%adjncy(k))
          end do
          u(i) = u(i) + s / A%values(A%kdiag(i))%re
        end do
        call A%graph%row_imap%gather_offp(u)
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

end module fdme_precon_gs_type
