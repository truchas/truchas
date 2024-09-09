!!
!! COMPLEX_VECTOR_FUNC_CONTAINERS
!!
!! This module defines several ad hoc containers for working with collections
!! of COMPLEX_VECTOR_FUNC objects whose dynamic types may differ.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module complex_vector_func_containers

  use complex_vector_func_class
  implicit none
  private

  public :: complex_vector_func ! re-export base class
  public :: complex_vector_func_box, complex_vector_func_list ! types
  public :: complex_vector_func_list_to_box_array ! procedure

  type :: complex_vector_func_box
    class(complex_vector_func), allocatable :: f
  contains
    procedure :: eval
    procedure :: eval_comp
  end type

  type :: complex_vector_func_list
    private
    integer :: n = 0
    type(list_func), pointer :: first => null()
  contains
    procedure :: append
    final :: complex_vector_func_list_delete
  end type

  type :: list_func
    class(complex_vector_func), allocatable :: f
    type(list_func), pointer :: next => null()
  end type list_func

contains

  function eval(this, x) result(fx)
    use,intrinsic :: iso_fortran_env, only: r8 => real64
    class(complex_vector_func_box), intent(in) :: this
    real(r8), intent(in) :: x(:)
    complex(r8) :: fx(this%f%dim)
    fx = this%f%eval(x)
  end function

  function eval_comp(this, i, x) result(fx)
    use,intrinsic :: iso_fortran_env, only: r8 => real64
    class(complex_vector_func_box), intent(in) :: this
    integer, intent(in) :: i
    real(r8), intent(in) :: x(:)
    complex(r8) :: fx(this%f%dim)
    fx = this%f%eval(x)
  end function

  !! Final subroutine for COMPLEX_VECTOR_FUNC_LIST objects.
  subroutine complex_vector_func_list_delete (this)
    type(complex_vector_func_list), intent(inout) :: this
    type(list_func), pointer :: rest
    do while (associated(this%first))
      rest => this%first%next
      deallocate(this%first)
      this%first => rest
    end do
  end subroutine

  subroutine append(this, f)
    class(complex_vector_func_list) :: this
    class(complex_vector_func), allocatable :: f
    type(list_func), pointer :: last
    ASSERT(allocated(f))
    if (associated(this%first)) then
      last => this%first
      do while (associated(last%next))
        last => last%next
      end do
      allocate(last%next)
      last => last%next
    else
      allocate(this%first)
      last => this%first
    end if
    call move_alloc (f, last%f)
    this%n = this%n + 1
  end subroutine

  subroutine complex_vector_func_list_to_box_array(list, array)
    type(complex_vector_func_list), intent(inout) :: list
    type(complex_vector_func_box), allocatable, intent(out) :: array(:)
    integer :: j
    type(list_func), pointer :: rest
    allocate(array(list%n))
    do j = 1, list%n
      ASSERT(associated(list%first))
      call move_alloc (list%first%f, array(j)%f)
      rest => list%first%next
      deallocate(list%first)
      list%first => rest
    end do
    list%n = 0
    ASSERT(.not.associated(list%first))
  end subroutine

end module complex_vector_func_containers
