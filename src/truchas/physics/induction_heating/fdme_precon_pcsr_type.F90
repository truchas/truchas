#include "f90_assert.fpp"

module fdme_precon_pcsr_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use fdme_precon_class
  use fdme_model_type
  use pcsr_matrix_type
  use pcsr_precon_class
  implicit none
  private

  type, extends(fdme_precon), public :: fdme_precon_pcsr
    type(pcsr_matrix), pointer :: matrix => null()  ! pointer to avoid dangling pointer
    class(pcsr_precon), allocatable :: precon
    real(r8) :: beta
  contains
    procedure :: init
    procedure :: setup
    procedure :: apply
    final :: fdme_precon_pcsr_delete
  end type

contains

  subroutine fdme_precon_pcsr_delete(this)
    type(fdme_precon_pcsr), intent(inout) :: this
    if (associated(this%matrix)) deallocate(this%matrix)
  end subroutine

  subroutine init(this, model, params, stat, errmsg)

    use parameter_list_type
    use pcsr_precon_factory, only: alloc_pcsr_precon

    class(fdme_precon_pcsr), intent(out) :: this
    type(fdme_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    this%model => model

    allocate(this%matrix)
    call this%matrix%init(this%model%A%graph, take_graph=.false.)

    call params%get('beta', this%beta, stat, errmsg, default=0.0_r8)
    if (stat /= 0) return

    call alloc_pcsr_precon(this%precon, this%matrix, params, stat, errmsg)
    if (stat /= 0) return

  end subroutine

  subroutine setup(this)
    class(fdme_precon_pcsr), intent(inout) :: this
    this%matrix%values(:) = this%model%A%values%re + this%beta*this%model%M%values%re
    call this%precon%compute
  end subroutine

  subroutine apply(this, x, y)
    class(fdme_precon_pcsr), intent(inout) :: this
    complex(r8), intent(in)  :: x(:)
    complex(r8), intent(out) :: y(:)
    y = x
    call this%precon%apply(y%re)
    call this%precon%apply(y%im)
  end subroutine

end module fdme_precon_pcsr_type
