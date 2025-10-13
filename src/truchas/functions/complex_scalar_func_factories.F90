!!
!! COMPLEX_SCALAR_FUNC_FACTORIES
!!
!! Procedures for instantiating new COMPLEX_SCALAR_FUNC class objects. All
!! a new CLASS(COMPLEX_SCALAR_FUNC) object, but the dynamic type is determined
!! by the particular function and arguments.
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

module complex_scalar_func_factories

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use complex_scalar_func_class
  use const_complex_scalar_func_type
  use dl_complex_scalar_func_type
  implicit none
  private

  public :: complex_scalar_func, alloc_complex_scalar_func
  
  !! Export type-specific allocation procedures from the specific type modules
  public :: alloc_const_complex_scalar_func, alloc_dl_complex_scalar_func

  interface alloc_complex_scalar_func
    procedure alloc_complex_scalar_func, get_complex_scalar_func
  end interface

contains

  !TODO: this needs error handling
  subroutine alloc_complex_scalar_func(f, params)

    use parameter_list_type

    class(complex_scalar_func), allocatable, intent(out) :: f
    type(parameter_list) :: params

    complex(r8) :: const
    character(:), allocatable :: ftype

    call params%get('type', ftype)
    select case (ftype)
    case ('constant')
      call params%get('value', const)
      call alloc_const_complex_scalar_func(f, const)
    case default
      INSIST(.false.) !TODO: need proper error handling
    end select
    ASSERT(allocated(f))

  end subroutine alloc_complex_scalar_func

  !! This subroutine gets the complex scalar function specified by the value of
  !! the parameter PARAM in the parameter list PLIST. The parameter value is
  !! either a complex array, a character string that is the name of a complex
  !! scalar function in the complex scalar function table, or a parameter list
  !! that defines the function.

  subroutine get_complex_scalar_func(plist, param, f, stat, errmsg)

    use parameter_list_type
    use func_table, only: lookup_func  !TODO: pass underlying object as argument

    type(parameter_list), intent(inout) :: plist
    character(*), intent(in) :: param
    class(complex_scalar_func), allocatable, intent(out) :: f
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    complex(r8) :: const
    character(:), allocatable :: fname
    type(parameter_list), pointer :: func_params

    if (plist%is_sublist(param)) then
      func_params => plist%sublist(param)
      call alloc_complex_scalar_func(f, func_params)  !TODO: should return stat, errmsg
    else if (plist%is_scalar(param)) then
#ifdef GNU_PR93762
      block
        character(:), allocatable :: dummy
        call plist%get(param, fname, stat, errmsg=dummy)
      end block
#else
      call plist%get(param, fname, stat)
#endif
      if (stat == 0) then ! name of a function
        call lookup_func(fname, f)
        if (.not.allocated(f)) then
          stat = 1
          errmsg = 'unknown complex scalar function name: ' // fname
          return
        end if
      else ! it might be a constant value
        call plist%get(param, const, stat, errmsg)
        if (stat /= 0) return
        call alloc_const_complex_scalar_func(f, const)
      end if
    else if (plist%is_vector(param)) then ! it might be complex constant as a 2-vector
      block
        class(*), allocatable :: tmp(:)
        call plist%get_any(param, tmp, stat, errmsg)
        if (stat /= 0) return
        stat = 1
        if (size(tmp) == 2) then
          select type (tmp)
          type is (real(r8))
            const = cmplx(tmp(1),tmp(2),kind=r8)
            stat = 0
          end select
        end if
        if (stat /= 0) then
          errmsg = 'not a complex(real64) parameter'
          return
        end if
        call alloc_const_complex_scalar_func(f, const)
      end block
    else
      stat = 1
      errmsg = 'invalid parameter value'
      return
    end if

  end subroutine get_complex_scalar_func

end module complex_scalar_func_factories
