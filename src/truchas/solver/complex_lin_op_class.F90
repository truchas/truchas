module complex_lin_op_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: complex_lin_op
  contains
    procedure(matvec), deferred :: matvec
    procedure(matvec), deferred :: precon
  end type

  abstract interface
    subroutine matvec(this, x, y)
      import complex_lin_op, r8
      class(complex_lin_op), intent(inout) :: this
      complex(r8) :: x(:), y(:)
    end subroutine
  end interface

end module complex_lin_op_class
