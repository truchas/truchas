#include "f90_assert.fpp"

module fdme_precon_boomer_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use new_fdme_precon_class
  use fdme_model_type
  use pcsr_matrix_type
  use pcsr_precon_boomer_type
  implicit none
  private

  type, extends(fdme_precon), public :: fdme_precon_boomer
    type(pcsr_matrix), pointer :: matrix => null()  ! pointer to avoid dangling pointer
    type(pcsr_precon_boomer) :: boomer
  contains
    procedure :: init
    procedure :: setup
    procedure :: apply
    final :: fdme_precon_boomer_delete
  end type

contains

  subroutine fdme_precon_boomer_delete(this)
    type(fdme_precon_boomer), intent(inout) :: this
    if (associated(this%matrix)) deallocate(this%matrix)
  end subroutine


  subroutine init(this, model, params, stat, errmsg)

    use parameter_list_type

    class(fdme_precon_boomer), intent(out) :: this
    type(fdme_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: plist

    this%model => model

    allocate(this%matrix)
    call this%matrix%init(this%model%A%graph, take_graph=.false.)

    plist => params%sublist('boomer') !FIXME
    call plist%set('num-cycles', 1) !HACK
    call this%boomer%init(this%matrix, plist, stat, errmsg)
    if (stat /= 0) return
    
  end subroutine

  subroutine setup(this)
    class(fdme_precon_boomer), intent(inout) :: this
    this%matrix%values(:) = this%model%A%values%re
    call this%boomer%compute
  end subroutine

  subroutine apply(this, x, y)
    class(fdme_precon_boomer), intent(inout) :: this
    complex(r8), intent(in)  :: x(:)
    complex(r8), intent(out) :: y(:)
    y = x
    call this%boomer%apply(y%re)
    call this%boomer%apply(y%im)
  end subroutine

end module fdme_precon_boomer_type
