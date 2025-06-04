module complex_lin_op_class

  use zvector_class
  implicit none
  private

  type, abstract, public :: complex_lin_op
  contains
    procedure(matvec), deferred :: matvec
    procedure(matvec), deferred :: precon
  end type

  abstract interface
    subroutine matvec(this, x, y)
      import complex_lin_op, zvector
      class(complex_lin_op), intent(inout) :: this
      class(zvector) :: x, y
    end subroutine
  end interface

end module complex_lin_op_class
