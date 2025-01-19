module fdme_precon_class

  use fdme_model_type
  implicit none
  private

  type, abstract, public :: fdme_precon
    type(fdme_model), pointer :: model => null() ! unowned reference
  contains
    procedure(init),  deferred :: init
    procedure(setup), deferred :: setup
    procedure(apply), deferred :: apply
  end type

  abstract interface
    subroutine init(this, model, params, stat, errmsg)
      use parameter_list_type
      import
      class(fdme_precon), intent(out) :: this
      type(fdme_model), intent(in), target :: model
      type(parameter_list), intent(inout) :: params
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
    end subroutine
    subroutine setup(this)
      import
      class(fdme_precon), intent(inout) :: this
    end subroutine
    subroutine apply(this, x, y)
      use,intrinsic :: iso_fortran_env, only: r8 => real64
      import
      class(fdme_precon), intent(inout) :: this
      complex(r8), intent(in)  :: x(:)
      complex(r8), intent(out) :: y(:)
    end subroutine
  end interface

end module fdme_precon_class
