!!
!! DL_SCALAR_FUNC_TYPE
!!
!! A concrete implementation of the abstract base class SCALAR_FUNC.  This
!! implementation defines a function dynamically loaded from a shared library.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Adapted for F2008, April 2014
!!

#include "f90_assert.fpp"

module dl_scalar_func_type

  use kinds, only: r8
  use scalar_func_class
  use fortran_dynamic_loader
  implicit none
  private

  type, extends(scalar_func), public :: dl_scalar_func
    ! private
    type(shlib) :: so = shlib()                   ! shared library handle
    procedure(f), pointer, nopass :: f => null()  ! the function
    real(r8), allocatable :: p(:)                 ! function parameters
  contains
    procedure :: eval
    final :: dl_scalar_func_delete
  end type dl_scalar_func

  abstract interface
    function f (x, p) result (fx)
      import r8
      real(r8), intent(in) :: x(*), p(*)
      real(r8) :: fx
    end function
  end interface

  !! Defined constructor
  interface dl_scalar_func
    procedure dl_scalar_func_value
  end interface

contains

  !! Final subroutine for DL_SCALAR_FUNC objects
  subroutine dl_scalar_func_delete (this)
    type(dl_scalar_func), intent(inout) :: this
    !call this%so%close () NNC: this does not play nice with assignment -- temporary fix
    this%f => null()
  end subroutine dl_scalar_func_delete

  !! Constructor for DL_SCALAR_FUNC objects
  function dl_scalar_func_value (lib, sym, p) result (dlf)

    use,intrinsic :: iso_c_binding, only: c_funptr, c_f_procpointer

    character(*), intent(in) :: lib, sym
    real(r8), intent(in), optional :: p(:)
    type(dl_scalar_func) :: dlf

    type(c_funptr) :: funptr

    !! Open the shared library.
    if (scan(lib, '/') == 0) then
      call dlf%so%open ('./' // lib, RTLD_NOW)
    else
      call dlf%so%open (lib, RTLD_NOW)
    end if

    !! Bind to the library function.
    call dlf%so%func (sym, funptr)
    call c_f_procpointer (funptr, dlf%f)

    !! Store the function parameters, if any.
    if (present(p)) then
      dlf%p = p
    else
      allocate(dlf%p(0))
    end if

  end function dl_scalar_func_value

  function eval (this, x) result (fx)
    class(dl_scalar_func), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: fx
    fx = this%f (x, this%p)
  end function eval

end module dl_scalar_func_type
