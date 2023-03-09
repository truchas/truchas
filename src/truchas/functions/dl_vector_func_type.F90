!!
!! DL_VECTOR_FUNC_TYPE
!!
!! A concrete implementation of the abstract base class VECTOR_FUNC.  This
!! implementation defines a function dynamically loaded from a shared library.
!!
!! Zach Jibben <zjibben@lanl.gov>
!! March 2023
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module dl_vector_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use vector_func_class
  use fortran_dynamic_loader
  implicit none
  private

  type, extends(vector_func), public :: dl_vector_func
    ! private
    type(shlib) :: so = shlib()                   ! shared library handle
    procedure(f), pointer, nopass :: f => null()  ! the function
    real(r8), allocatable :: p(:)                 ! function parameters
  contains
    procedure :: eval
    procedure :: eval_comp
    final :: dl_vector_func_delete
  end type dl_vector_func

  abstract interface
    subroutine f(x, p, r)
      import r8
      real(r8), intent(in) :: x(*), p(*)
      real(r8), intent(out) :: r(*)
    end subroutine
  end interface

  !! Defined constructor
  interface dl_vector_func
    procedure dl_vector_func_value
  end interface

contains

  !! Final subroutine for DL_VECTOR_FUNC objects
  subroutine dl_vector_func_delete(this)
    type(dl_vector_func), intent(inout) :: this
    !call this%so%close() NNC: this does not play nice with assignment -- temporary fix
    this%f => null()
  end subroutine dl_vector_func_delete

  !! Constructor for DL_VECTOR_FUNC objects
  function dl_vector_func_value(lib, sym, dim, p) result(dlf)

    use,intrinsic :: iso_c_binding, only: c_funptr, c_f_procpointer

    character(*), intent(in) :: lib, sym
    integer, intent(in) :: dim
    real(r8), intent(in), optional :: p(:)
    type(dl_vector_func) :: dlf

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

    dlf%dim = dim

  end function dl_vector_func_value

  function eval(this, x) result(fx)
    class(dl_vector_func), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: fx(this%dim)
    call this%f(x, this%p, fx)
  end function eval

  function eval_comp(this, i, x) result(fxi)
    class(dl_vector_func), intent(in) :: this
    integer,  intent(in) :: i
    real(r8), intent(in) :: x(:)
    real(r8) :: fxi, fx(this%dim)
    fx = this%eval(x)
    fxi = fx(i)
  end function eval_comp

end module dl_vector_func_type
