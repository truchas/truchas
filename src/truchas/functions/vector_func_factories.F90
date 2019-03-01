!!
!! VECTOR_FUNC_FACTORIES
!!
!! Procedures for instantiating new VECTOR_FUNC class objects.  All return
!! return a new CLASS(VECTOR_FUNC) object, but the dynamic type is determined
!! by the particular function and arguments.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! April 2014
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

module vector_func_factories

  use kinds, only: r8
  use vector_func_class
  implicit none
  private

  public :: vector_func ! re-export

  !! These functions return CLASS(VECTOR_FUNC) pointers
  public :: new_const_vector_func
  public :: new_tabular_vector_func

  !! These subroutines allocate an allocatable CLASS(VECTOR_FUNC) argument
  public :: alloc_const_vector_func
  public :: alloc_tabular_vector_func

contains

  subroutine alloc_const_vector_func (f, const)
    use const_vector_func_type
    class(vector_func), allocatable, intent(out) :: f
    real(r8), intent(in) :: const(:)
    allocate(f, source=const_vector_func(const))
  end subroutine alloc_const_vector_func

  subroutine alloc_tabular_vector_func (f, x, y)
    use tabular_vector_func_type
    class(vector_func), allocatable, intent(out) :: f
    real(r8), intent(in) :: x(:), y(:,:)
    allocate(f, source=tabular_vector_func(x, y))
  end subroutine alloc_tabular_vector_func

  function new_const_vector_func (const) result (f)
    use const_vector_func_type
    real(r8), intent(in) :: const(:)
    class(vector_func), pointer :: f
    allocate(f, source=const_vector_func(const))
  end function new_const_vector_func

  function new_tabular_vector_func (x, y) result (f)
    use tabular_vector_func_type
    real(r8), intent(in) :: x(:), y(:,:)
    class(vector_func), pointer :: f
    allocate(f, source=tabular_vector_func(x,y))
  end function new_tabular_vector_func


  subroutine alloc_vector_func (f, params)

    use parameter_list_type

    class(vector_func), allocatable, intent(out) :: f
    type(parameter_list) :: params

    real(r8), allocatable :: v0(:), x(:), y(:,:)
    character(:), allocatable :: ftype

    call params%get ('type', ftype)
    select case (ftype)
    case ('constant')
      call params%get ('value', v0)
      call alloc_const_vector_func (f, v0)
    case ('tabular')
      call params%get ('tabular-x', x)
      call params%get ('tabular-y', y)
      INSIST(size(x) == size(y,dim=2))  !TODO: need proper error handling
      call alloc_tabular_vector_func (f, x, y)
    case default
      INSIST(.false.) !TODO: need proper error handling
    end select
    ASSERT(allocated(f))

  end subroutine alloc_vector_func

end module vector_func_factories
