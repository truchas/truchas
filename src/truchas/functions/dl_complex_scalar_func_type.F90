!!
!! DL_COMPLEX_SCALAR_FUNC_TYPE
!!
!! A concrete implementation of the abstract base class COMPLEX_SCALAR_FUNC.
!! This implementation defines a function dynamically loaded from a shared
!! object library.
!!
!! Neil N. Carlson <neil.n.carlson@gmail.com>
!! December 2024
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module dl_complex_scalar_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use complex_scalar_func_class
  use fortran_dynamic_loader
  implicit none
  private

  public :: alloc_dl_complex_scalar_func

  type, extends(complex_scalar_func), public :: dl_complex_scalar_func
    ! private
    type(shlib) :: so = shlib()                   ! shared library handle
    procedure(f), pointer, nopass :: f => null()  ! the function
    real(r8), allocatable :: p(:)                 ! function parameters
  contains
    procedure :: eval
    final :: dl_complex_scalar_func_delete
  end type

  abstract interface
    subroutine f(x, p, r)
      import r8
      real(r8), intent(in) :: x(*), p(*)
      complex(r8) :: r
    end subroutine
  end interface

  !! Defined constructor
  interface dl_complex_scalar_func
    procedure dl_complex_scalar_func_value
  end interface

contains

  subroutine alloc_dl_complex_scalar_func(f, lib, sym, p)
    class(complex_scalar_func), allocatable, intent(out) :: f
    character(*), intent(in) :: lib, sym
    real(r8), intent(in), optional :: p(:)
    allocate(f, source=dl_complex_scalar_func(lib, sym, p))
  end subroutine

  !! Final subroutine for DL_SCALAR_FUNC objects
  subroutine dl_complex_scalar_func_delete(this)
    type(dl_complex_scalar_func), intent(inout) :: this
    !call this%so%close() NNC: this does not play nice with assignment -- temporary fix
    this%f => null()
  end subroutine

  !! Constructor for DL_SCALAR_FUNC objects
  function dl_complex_scalar_func_value(lib, sym, p) result(dlf)

    use,intrinsic :: iso_c_binding, only: c_funptr, c_f_procpointer

    character(*), intent(in) :: lib, sym
    real(r8), intent(in), optional :: p(:)
    type(dl_complex_scalar_func) :: dlf

    type(c_funptr) :: funptr

    !! Open the shared library.
    if (scan(lib, '/') == 0) then
      call dlf%so%open('./' // lib, RTLD_NOW)
    else
      call dlf%so%open(lib, RTLD_NOW)
    end if

    !! Bind to the library function.
    call dlf%so%func(sym, funptr)
    call c_f_procpointer(funptr, dlf%f)

    !! Store the function parameters, if any.
    if (present(p)) then
      dlf%p = p
    else
      allocate(dlf%p(0))
    end if

  end function dl_complex_scalar_func_value

  function eval(this, x) result(fx)
    class(dl_complex_scalar_func), intent(in) :: this
    real(r8), intent(in) :: x(:)
    complex(r8) :: fx
    call this%f(x, this%p, fx)
  end function

end module dl_complex_scalar_func_type
