!!
!! VECTOR_FUNC_CONTAINERS
!!
!! This module defines several containers for working with collections of
!! VECTOR_FUNC objects whose dynamic types may differ.  The functionality
!! provided is limited to just that required by the needs of the immediate
!! application code; this does not aim to be anything general.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! February 2014
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  A polymorphic VECTOR_FUNC array must have a single dynamic type, and thus
!!  cannot be used to hold a collection of VECTOR_FUNC objects whose dynamic
!!  types may differ.  This module provides two derived types for addressing
!!  this problem.
!!
!!  The type VECTOR_FUNC_BOX wraps a public allocatable VECTOR_FUNC component:
!!
!!    TYPE VECTOR_FUNC_BOX
!!      CLASS(VECTOR_FUNC), ALLOCATABLE :: F
!!    END TYPE
!!
!!  Use an array of this type to hold a collection of VECTOR_FUNC objects of
!!  differing dynamic type.
!!
!!  The type VECTOR_FUNC_LIST provides a dynamic linked-list of VECTOR_FUNC
!!  objects.  Its components are private and it has a single type bound
!!  procedure APPEND(F) which appends the passed VECTOR_FUNC object F to the
!!  list.  The passed F must be allocatable, and the procedure moves the
!!  allocation from F to an internal list object; F is returned unallocated.
!!
!!  The subroutine VECTOR_FUNC_LIST_TO_BOX_ARRAY (LIST, ARRAY) converts the
!!  VECTOR_FUNC_LIST object LIST into a rank-1 VECTOR_FUNC_BOX array object.
!!  ARRAY must be allocatable; it is allocated to the correct size by the
!!  subroutine.  LIST is returned empty, the allocation for its VECTOR_FUNC
!!  objects having been moved into the returned ARRAY.
!!
!!  The type VECTOR_FUNC_VLIST provides a dynamic linked-list of vectors of
!!  VECTOR_FUNC objects, or more precisely rank-1 VECTOR_FUNC_BOX arrays.
!!  Its components are private and it has a single type bound procedure
!!  APPEND(V) which appends the passed rank-1 VECTOR_FUNC_BOX array V to the
!!  list.  The passed V must be allocatable, and the procedure moves the
!!  allocation from V to in internal list object; V is returned unallocated.
!!  All vectors V must have the same size; the first vector appended
!!  establishes that size.
!!
!!  The subroutine VECTOR_FUNC_VLIST_TO_BOX_ARRAY (VLIST, ARRAY) converts the
!!  VECTOR_FUNC_VLIST object VLIST into a rank-2 VECTOR_FUNC_BOX array object.
!!  ARRAY must be allocatable; is is allocated to the correct size by the
!!  subroutine.  VLIST is returned empty, the allocation for its VECTOR_FUNC
!!  objects having been moved into the returned ARRAY.  The first dimension of
!!  ARRAY corresponds to the elements of the vectors.
!!

#include "f90_assert.fpp"

module vector_func_containers

  use vector_func_class
  implicit none
  private

  public :: vector_func ! re-export (necessary?)
  public :: vector_func_box, vector_func_list, vector_func_vlist ! types
  public :: vector_func_list_to_box_array, vector_func_vlist_to_box_array ! procedures

  type :: vector_func_box
    class(vector_func), allocatable :: f
  contains
    procedure :: eval => vector_func_box_eval
    procedure :: eval_comp => vector_func_box_eval_comp
  end type vector_func_box

  type :: vector_func_list
    private
    integer :: n = 0
    type(list_func), pointer :: first => null()
  contains
    procedure :: append => vector_func_list_append
    final :: vector_func_list_delete
  end type vector_func_list

  type :: list_func
    class(vector_func), allocatable :: f
    type(list_func), pointer :: next => null()
  end type list_func

  type :: vector_func_vlist
    private
    integer :: m = 0, n = 0
    type(vlist_func), pointer :: first => null()
  contains
    procedure :: append => vector_func_vlist_append
    final :: vector_func_vlist_delete
  end type vector_func_vlist

  type :: vlist_func
    type(vector_func_box), allocatable :: v(:)
    type(vlist_func), pointer :: next => null()
  end type vlist_func

contains

  function vector_func_box_eval (this, x) result (fx)
    use,intrinsic :: iso_fortran_env, only: r8 => real64
    class(vector_func_box), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: fx(this%f%dim)
    fx = this%f%eval(x)
  end function vector_func_box_eval

  function vector_func_box_eval_comp (this, i, x) result (fx)
    use,intrinsic :: iso_fortran_env, only: r8 => real64
    class(vector_func_box), intent(in) :: this
    integer, intent(in) :: i
    real(r8), intent(in) :: x(:)
    real(r8) :: fx(this%f%dim)
    fx = this%f%eval(x)
  end function vector_func_box_eval_comp


  !! Final subroutine for VECTOR_FUNC_LIST objects.
  subroutine vector_func_list_delete (this)
    type(vector_func_list), intent(inout) :: this
    type(list_func), pointer :: rest
    do while (associated(this%first))
      rest => this%first%next
      deallocate(this%first)
      this%first => rest
    end do
  end subroutine

  subroutine vector_func_list_append (this, f)
    class(vector_func_list) :: this
    class(vector_func), allocatable :: f
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
  end subroutine vector_func_list_append

  subroutine vector_func_list_to_box_array (list, array)
    type(vector_func_list) :: list
    type(vector_func_box), allocatable, intent(out) :: array(:)
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
  end subroutine vector_func_list_to_box_array

  !! Final subroutine for VECTOR_FUNC_VLIST objects.
  subroutine vector_func_vlist_delete (this)
    type(vector_func_vlist), intent(inout) :: this
    type(vlist_func), pointer :: rest
    do while (associated(this%first))
      rest => this%first%next
      deallocate(this%first)
      this%first => rest
    end do
  end subroutine

  subroutine vector_func_vlist_append (this, v)
    class(vector_func_vlist) :: this
    type(vector_func_box), allocatable :: v(:)
    type(vlist_func), pointer :: last
    ASSERT(allocated(v))
    if (associated(this%first)) then
      INSIST(size(v) == this%m)
      last => this%first
      do while (associated(last%next))
        last => last%next
      end do
      allocate(last%next)
      last => last%next
    else
      this%m = size(v)
      allocate(this%first)
      last => this%first
    end if
    call move_alloc (v, last%v)
    this%n = this%n + 1
  end subroutine vector_func_vlist_append

  subroutine vector_func_vlist_to_box_array (vlist, array)
    type(vector_func_vlist) :: vlist
    type(vector_func_box), allocatable, intent(out) :: array(:,:)
    integer :: i, j
    type(vlist_func), pointer :: rest
    allocate(array(vlist%m,vlist%n))
    do j = 1, vlist%n
      ASSERT(associated(vlist%first))
      do i = 1, vlist%m
        call move_alloc (vlist%first%v(i)%f, array(i,j)%f)
      end do
      rest => vlist%first%next
      deallocate(vlist%first)
      vlist%first => rest
    end do
    ASSERT(.not.associated(vlist%first))
    vlist%m = 0
    vlist%n = 0
  end subroutine vector_func_vlist_to_box_array

end module vector_func_containers
