!!
!! SCALAR_FUNC_FACTORIES
!!
!! Procedures for instantiating new SCALAR_FUNC class objects.  All return
!! return a new CLASS(SCALAR_FUNC) object, but the dynamic type is determined
!! by the particular function and arguments.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Adapted for F2008, April 2014
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! Two flavors of procedures are provided.  In one, the functions return a
!! pointer to a newly allocated CLASS(SCALAR_FUNC) object, and in the second,
!! subroutines allocate an allocatable CLASS(SCALAR_FUNC) argument.  It is
!! typical for a CLASS(SCALAR_FUNC) object to be passed from one procedure
!! to another before ultimately being stored as a data structure component.
!! Pointers have been used to do this to avoid copying of the objects, however
!! the new MOVE_ALLOC instrinsic allows allocatable variables to be used for
!! the same purpose, bringing the advantages of allocatable over pointer.
!!

#include "f90_assert.fpp"

module scalar_func_factories

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use scalar_func_class
  implicit none
  private

  public :: scalar_func ! re-export

  !! These return CLASS(SCALAR_FUNC) pointers
  public :: new_const_scalar_func
  public :: new_fptr_scalar_func
  public :: new_mpoly_scalar_func
  public :: new_poly_scalar_func
  public :: new_tabular_scalar_func, new_tabular_ad_scalar_func
#ifdef ENABLE_DYNAMIC_LOADING
  public :: new_dl_scalar_func
#endif

  !! These subroutines allocate an allocatable CLASS(SCALAR_FUNC) argument
  public :: alloc_const_scalar_func
  public :: alloc_fptr_scalar_func
  public :: alloc_mpoly_scalar_func
  public :: alloc_poly_scalar_func
  public :: alloc_sstep_scalar_func
  public :: alloc_tabular_scalar_func, alloc_tabular_ad_scalar_func
#ifdef ENABLE_DYNAMIC_LOADING
  public :: alloc_dl_scalar_func
#endif

  !! These higher-level procedures take a parameter list as input.
  public :: new_scalar_func
  public :: alloc_scalar_func

  interface alloc_scalar_func
    procedure alloc_scalar_func, get_scalar_func
  end interface

contains

  subroutine alloc_const_scalar_func(f, const)
    use const_scalar_func_type
    class(scalar_func), allocatable, intent(out) :: f
    real(r8), intent(in) :: const
    allocate(f, source=const_scalar_func(const))
  end subroutine alloc_const_scalar_func

#ifdef ENABLE_DYNAMIC_LOADING
  subroutine alloc_dl_scalar_func(f, lib, sym, p)
    use dl_scalar_func_type
    class(scalar_func), allocatable, intent(out) :: f
    character(*), intent(in) :: lib, sym
    real(r8), intent(in), optional :: p(:)
    allocate(f, source=dl_scalar_func(lib, sym, p))
  end subroutine alloc_dl_scalar_func
#endif

  subroutine alloc_fptr_scalar_func(f, fptr, p)
    use fptr_scalar_func_type
    class(scalar_func), allocatable, intent(out) :: f
    procedure(fptr_func) :: fptr
    real(r8), intent(in), optional :: p(:)
    allocate(f, source=fptr_scalar_func(fptr, p))
  end subroutine alloc_fptr_scalar_func

  subroutine alloc_mpoly_scalar_func(f, c, e, x0)
    use mpoly_scalar_func_type
    class(scalar_func), allocatable, intent(out) :: f
    real(r8), intent(in) :: c(:)
    integer,  intent(in) :: e(:,:)
    real(r8), intent(in), optional :: x0(:)
    allocate(f, source=mpoly_scalar_func(c, e, x0))
  end subroutine alloc_mpoly_scalar_func

  subroutine alloc_poly_scalar_func(f, c, e, x0)
    use poly_scalar_func_type
    class(scalar_func), allocatable, intent(out) :: f
    real(r8), intent(in) :: c(:)
    integer,  intent(in) :: e(:)
    real(r8), intent(in), optional :: x0
    allocate(f, source=poly_scalar_func(c, e, x0))
  end subroutine alloc_poly_scalar_func

  subroutine alloc_tabular_scalar_func(f, x, y, dim, smooth, extrap)
    use tabular_scalar_func_type
    class(scalar_func), allocatable, intent(out) :: f
    real(r8), intent(in) :: x(:), y(:)
    integer, intent(in), optional :: dim
    logical, intent(in), optional :: smooth
    character(*), intent(in), optional :: extrap
    allocate(f, source=tabular_scalar_func(x, y, dim, smooth, extrap))
  end subroutine alloc_tabular_scalar_func

  subroutine alloc_tabular_ad_scalar_func(f, df, x0, y0)
    use tabular_scalar_func_type
    class(scalar_func), allocatable, intent(out) :: f
    type(tabular_scalar_func), intent(in) :: df
    real(r8), intent(in) :: x0, y0
    allocate(f, source=tabular_ad_scalar_func(df, x0, y0))
  end subroutine alloc_tabular_ad_scalar_func

  ! A temporary stand-in (?) until we move to parameter list driven instantiation
  subroutine alloc_sstep_scalar_func(f, x0, y0, x1, y1)
    use fptr_scalar_func_type
    class(scalar_func), allocatable, intent(out) :: f
    real(r8), intent(in) :: x0, y0, x1, y1
    call alloc_fptr_scalar_func(f, smooth_ramp, [x0, y0, x1, y1])
  end subroutine alloc_sstep_scalar_func

  function new_const_scalar_func(const) result(f)
    use const_scalar_func_type
    real(r8), intent(in) :: const
    class(scalar_func), pointer :: f
    allocate(f, source=const_scalar_func(const))
  end function new_const_scalar_func

#ifdef ENABLE_DYNAMIC_LOADING
  function new_dl_scalar_func(lib, sym, p) result(f)
    use dl_scalar_func_type
    character(*), intent(in) :: lib, sym
    real(r8), intent(in), optional :: p(:)
    class(scalar_func), pointer :: f
    allocate(f, source=dl_scalar_func(lib, sym, p))
  end function new_dl_scalar_func
#endif

  function new_fptr_scalar_func(fptr, p) result(f)
    use fptr_scalar_func_type
    procedure(fptr_func) :: fptr
    real(r8), intent(in), optional :: p(:)
    class(scalar_func), pointer :: f
    allocate(f, source=fptr_scalar_func(fptr, p))
  end function new_fptr_scalar_func

  function new_mpoly_scalar_func(c, e, x0) result(f)
    use mpoly_scalar_func_type
    real(r8), intent(in) :: c(:)
    integer,  intent(in) :: e(:,:)
    real(r8), intent(in), optional :: x0(:)
    class(scalar_func), pointer :: f
    allocate(f, source=mpoly_scalar_func(c, e, x0))
  end function new_mpoly_scalar_func

  function new_poly_scalar_func(c, e, x0) result(f)
    use poly_scalar_func_type
    real(r8), intent(in) :: c(:)
    integer,  intent(in) :: e(:)
    real(r8), intent(in), optional :: x0
    class(scalar_func), pointer :: f
    allocate(f, source=poly_scalar_func(c, e, x0))
  end function new_poly_scalar_func

  function new_tabular_scalar_func(x, y, dim, smooth, extrap) result(f)
    use tabular_scalar_func_type
    real(r8), intent(in) :: x(:), y(:)
    integer, intent(in), optional :: dim
    logical, intent(in), optional :: smooth
    character(*), intent(in), optional :: extrap
    class(scalar_func), pointer :: f
    allocate(f, source=tabular_scalar_func(x, y, dim, smooth, extrap))
  end function new_tabular_scalar_func

  function new_tabular_ad_scalar_func(df, x0, y0) result(f)
    use tabular_scalar_func_type
    type(tabular_scalar_func), intent(in) :: df
    real(r8), intent(in) :: x0, y0
    class(scalar_func), pointer :: f
    allocate(f, source=tabular_ad_scalar_func(df, x0, y0))
  end function new_tabular_ad_scalar_func

  !! A higher-level constructors that takes a parameter list
  !! that describes the function to be instantiated.

  function new_scalar_func(params) result(f)

    use parameter_list_type
    type(parameter_list), intent(inout) :: params
    class(scalar_func), pointer :: f

    real(r8) :: x0, t0, v0, t1, v1
    integer, allocatable :: expon(:), expon2(:,:)
    real(r8), allocatable :: coef(:), x(:), y(:), x0_def(:), x02(:)
    character(:), allocatable :: ftype

    call params%get('type', ftype)
    select case (ftype)
    case ('constant')
      call params%get('value', v0)
      f => new_const_scalar_func(v0)
    case ('polynomial')
      call params%get ('poly-coef', coef)
      if (params%is_vector('poly-powers')) then
        call params%get('poly-powers', expon)
        call params%get('poly-center', x0, default=0.0_r8)
        ASSERT(size(coef) == size(expon))
        f => new_poly_scalar_func(coef, expon, x0)
      else if (params%is_matrix('poly-powers')) then
        call params%get('poly-powers', expon2)
        ASSERT(size(coef) == size(expon2,2))
        allocate(x0_def(size(expon2,1)))
        x0_def = 0.0_r8
        call params%get ('poly-center', x02, default=x0_def)
        f => new_mpoly_scalar_func(coef, expon2, x02)
      else
        INSIST(.false.)
      end if
    case ('tabular')
      call params%get('tabular x', x)
      call params%get('tabular y', y)
      ASSERT(size(x) == size(y))
      f => new_tabular_scalar_func(x, y)
    case ('smooth ramp')
      call params%get('begin time', t0)
      call params%get('begin value', v0)
      call params%get('end time', t1)
      call params%get('end value', v1)
      f => new_fptr_scalar_func(smooth_ramp,[t0,v0,t1,v1])
    case default
      INSIST(.false.)
    end select
    ASSERT(associated(f))

  end function new_scalar_func


  subroutine alloc_scalar_func(f, params)

    use parameter_list_type
    class(scalar_func), allocatable, intent(out) :: f
    type(parameter_list), intent(inout) :: params

    real(r8) :: x0, t0, v0, t1, v1
    integer, allocatable :: expon(:), expon2(:,:)
    real(r8), allocatable :: coef(:), x(:), y(:), x0_def(:), x02(:)
    character(:), allocatable :: ftype

    call params%get('type', ftype)
    select case (ftype)
    case ('constant')
      call params%get('value', v0)
      call alloc_const_scalar_func(f, v0)
    case ('polynomial')
      call params%get('poly-coef', coef)
      if (params%is_vector('poly-powers')) then
        call params%get('poly-powers', expon)
        call params%get('poly-center', x0, default=0.0_r8)
        ASSERT(size(coef) == size(expon))
        call alloc_poly_scalar_func (f, coef, expon, x0)
      else if (params%is_matrix('poly-powers')) then
        call params%get('poly-powers', expon2)
        ASSERT(size(coef) == size(expon2,2))
        allocate(x0_def(size(expon2,1)))
        x0_def = 0.0_r8
        call params%get ('poly-center', x02, default=x0_def)
        call alloc_mpoly_scalar_func(f, coef, expon2, x02)
      else
        INSIST(.false.)
      end if
    case ('tabular')
      call params%get('tabular x', x)
      call params%get('tabular y', y)
      ASSERT(size(x) == size(y))
      call alloc_tabular_scalar_func(f, x, y)
    case ('smooth ramp')
      call params%get('begin time', t0)
      call params%get('begin value', v0)
      call params%get('end time', t1)
      call params%get('end value', v1)
      call alloc_fptr_scalar_func(f, smooth_ramp,[t0,v0,t1,v1])
    case default
      INSIST(.false.)
    end select
    ASSERT(allocated(f))

  end subroutine alloc_scalar_func

  real(r8) function smooth_ramp(x, p)
    real(r8), intent(in) :: x(*), p(*)
    real(r8) :: s
    associate (t => x(1), t0 => p(1), v0 => p(2), t1 => p(3), v1 => p(4))
      if (t <= t0) then
        smooth_ramp = v0
      else if (t >= t1) then
        smooth_ramp = v1
      else
        s = (t - t0) / (t1 - t0)
        smooth_ramp = v0 + (v1-v0) * (s**2) * (3 - 2*s)
      end if
    end associate
  end function smooth_ramp

  !! This subroutine gets the scalar function specified by the value of the
  !! parameter PARAM in the parameter list PLIST. The parameter value is either
  !! a real scalar, a character string that is the name of a function in the
  !! function table, or a parameter list that defines the function.

  subroutine get_scalar_func(plist, param, f, stat, errmsg)

    use parameter_list_type
    use func_table, only: lookup_func  !TODO: pass underlying object as argument

    type(parameter_list), intent(inout) :: plist
    character(*), intent(in) :: param
    class(scalar_func), allocatable, intent(out) :: f
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8) :: const
    character(:), allocatable :: fname
    type(parameter_list), pointer :: func_params

    if (plist%is_sublist(param)) then
      func_params => plist%sublist(param)
      call alloc_scalar_func(f, func_params)  !TODO: should return stat, errmsg
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
          errmsg = 'unknown function name: ' // fname
          return
        end if
      else  ! it must be a constant value
        call plist%get(param, const, stat, errmsg)
        if (stat /= 0) return
        call alloc_const_scalar_func(f, const)
      end if
    else
      stat = 1
      errmsg = 'invalid parameter value'
      return
    end if

  end subroutine get_scalar_func

end module scalar_func_factories
