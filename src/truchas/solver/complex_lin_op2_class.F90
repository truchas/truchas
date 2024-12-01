module complex_lin_op2_class

  use zvector_class
  implicit none
  private

  type, abstract, public :: complex_lin_op2
  contains
    procedure(matvec), deferred :: matvec
    procedure(matvec), deferred :: precon
  end type

  abstract interface
    subroutine matvec(this, x, y)
      import complex_lin_op2, zvector
      class(complex_lin_op2), intent(inout) :: this
      class(zvector) :: x, y
    end subroutine
  end interface

end module complex_lin_op2_class
