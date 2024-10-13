#include "f90_assert.fpp"

module fdme_minres_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use complex_lin_op_class
  use fdme_vector_type
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
    call this%model%matvec2(x, y)
  end subroutine

  subroutine precon(this, x, y)
    class(fdme_lin_op), intent(inout) :: this
    complex(r8) :: x(:), y(:)
    if (.not.allocated(this%dinv)) then
      block
        integer :: j
        associate (A => this%model%AA)
          call A%kdiag_init
          allocate(this%dinv(A%nrow))
          do j = 1, A%nrow
            this%dinv(j) = 1.0_r8 / A%values(A%kdiag(j))%re
          end do
        end associate
      end block
    end if
    !y = x ! no preconditioning
    y = this%dinv*x ! diagonal preconditioning
    !block ! doesn't work with hiptmair (doesn't satisfy requirements)
    !  real(r8) :: xarray(2,size(x))
    !  xarray(1,:) = x%re; xarray(2,:) = x%im
    !  call this%my_precon%apply(xarray)
    !  y%re = xarray(1,:); y%im = xarray(2,:)
    !end block
  end subroutine

  subroutine init(this, vec, model, precon, params, stat, errmsg)
    use parameter_list_type
    class(fdme_minres_solver), intent(out) :: this
    type(fdme_vector), intent(in) :: vec
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
    type(fdme_vector), intent(inout) :: efield
    integer, intent(out) :: stat
    complex(r8), allocatable :: x(:)!, b(:)
    integer :: n
    allocate(x(size(efield%array,dim=2)))!, b(size(this%model%rhs%array,dim=2)))
    call this%minres%solve(this%lin_op, this%model%crhs, x)
    efield%array(1,:) = x%re; efield%array(2,:) = x%im
  end subroutine

end module fdme_minres_solver_type
