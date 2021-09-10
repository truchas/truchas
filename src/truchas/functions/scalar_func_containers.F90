!!
!! SCALAR_FUNC_CONTAINERS
!!
!! This module defines several containers for working with collections of
!! SCALAR_FUNC objects whose dynamic types may differ.  The functionality
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
!!  A polymorphic SCALAR_FUNC array must have a single dynamic type, and thus
!!  cannot be used to hold a collection of SCALAR_FUNC objects whose dynamic
!!  types may differ.  This module provides two derived types for addressing
!!  this problem.
!!
!!  The type SCALAR_FUNC_BOX wraps a public allocatable SCALAR_FUNC component:
!!
!!    TYPE SCALAR_FUNC_BOX
!!      CLASS(SCALAR_FUNC), ALLOCATABLE :: F
!!    END TYPE
!!
!!  Use an array of this type to hold a collection of SCALAR_FUNC objects of
!!  differing dynamic type.
!!
!!  The type SCALAR_FUNC_LIST provides a dynamic linked-list of SCALAR_FUNC
!!  objects.  Its components are private and it has a single type bound
!!  procedure APPEND(F) which appends the passed SCALAR_FUNC object F to the
!!  list.  The passed F must be allocatable, and the procedure moves the
!!  allocation from F to an internal list object; F is returned unallocated.
!!
!!  The subroutine SCALAR_FUNC_LIST_TO_BOX_ARRAY(LIST, ARRAY) converts the
!!  SCALAR_FUNC_LIST object LIST into a rank-1 SCALAR_FUNC_BOX array object.
!!  ARRAY must be allocatable; it is allocated to the correct size by the
!!  subroutine.  LIST is returned empty, the allocation for its SCALAR_FUNC
!!  objects having been moved into the returned ARRAY.
!!
!!  The type SCALAR_FUNC_VLIST provides a dynamic linked-list of vectors of
!!  SCALAR_FUNC objects, or more precisely rank-1 SCALAR_FUNC_BOX arrays.
!!  Its components are private and it has a single type bound procedure
!!  APPEND(V) which appends the passed rank-1 SCALAR_FUNC_BOX array V to the
!!  list.  The passed V must be allocatable, and the procedure moves the
!!  allocation from V to in internal list object; V is returned unallocated.
!!  All vectors V must have the same size; the first vector appended
!!  establishes that size.
!!
!!  The subroutine SCALAR_FUNC_VLIST_TO_BOX_ARRAY(VLIST, ARRAY) converts the
!!  SCALAR_FUNC_VLIST object VLIST into a rank-2 SCALAR_FUNC_BOX array object.
!!  ARRAY must be allocatable; is is allocated to the correct size by the
!!  subroutine.  VLIST is returned empty, the allocation for its SCALAR_FUNC
!!  objects having been moved into the returned ARRAY.  The first dimension of
!!  ARRAY corresponds to the elements of the vectors.
!!

#include "f90_assert.fpp"

module scalar_func_containers

  use scalar_func_class
  implicit none
  private

  public :: scalar_func ! re-export (necessary?)
  public :: scalar_func_box, scalar_func_ptr, scalar_func_list, scalar_func_vlist ! types
  public :: scalar_func_list_to_box_array, scalar_func_vlist_to_box_array ! procedures

  type :: scalar_func_box
    class(scalar_func), allocatable :: f
  contains
    procedure :: eval => scalar_func_box_eval
  end type scalar_func_box

  type :: scalar_func_ptr
    class(scalar_func), pointer :: f => null() ! unowned reference
  contains
    procedure :: eval => scalar_func_ptr_eval
  end type

  type :: scalar_func_list
    private
    integer :: n = 0
    type(list_func), pointer :: first => null()
  contains
    procedure :: append => scalar_func_list_append
    final :: scalar_func_list_delete
  end type

  type :: list_func
    class(scalar_func), allocatable :: f
    type(list_func), pointer :: next => null()
  end type

  type :: scalar_func_vlist
    private
    integer :: m = 0, n = 0
    type(vlist_func), pointer :: first => null()
  contains
    procedure :: append => scalar_func_vlist_append
    final :: scalar_func_vlist_delete
  end type

  type :: vlist_func
    type(scalar_func_box), allocatable :: v(:)
    type(vlist_func), pointer :: next => null()
  end type

contains

  function scalar_func_box_eval(this, x) result(fx)
    use kinds, only: r8
    class(scalar_func_box), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: fx
    fx = this%f%eval(x)
  end function scalar_func_box_eval

  function scalar_func_ptr_eval(this, x) result(fx)
    use kinds, only: r8
    class(scalar_func_ptr), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: fx
    fx = this%f%eval(x)
  end function scalar_func_ptr_eval

  !! Final subroutine for SCALAR_FUNC_LIST objects.
  subroutine scalar_func_list_delete(this)
    type(scalar_func_list), intent(inout) :: this
    type(list_func), pointer :: rest
    do while (associated(this%first))
      rest => this%first%next
      deallocate(this%first)
      this%first => rest
    end do
  end subroutine

  subroutine scalar_func_list_append(this, f)
    class(scalar_func_list) :: this
    class(scalar_func), allocatable :: f
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
    call move_alloc(f, last%f)
    this%n = this%n + 1
  end subroutine scalar_func_list_append

  subroutine scalar_func_list_to_box_array(list, array)
    type(scalar_func_list) :: list
    type(scalar_func_box), allocatable, intent(out) :: array(:)
    integer :: j
    type(list_func), pointer :: rest
    allocate(array(list%n))
    do j = 1, list%n
      ASSERT(associated(list%first))
      call move_alloc(list%first%f, array(j)%f)
      rest => list%first%next
      deallocate(list%first)
      list%first => rest
    end do
    list%n = 0
    ASSERT(.not.associated(list%first))
  end subroutine scalar_func_list_to_box_array

  !! Final subroutine for SCALAR_FUNC_VLIST objects.
  subroutine scalar_func_vlist_delete(this)
    type(scalar_func_vlist), intent(inout) :: this
    type(vlist_func), pointer :: rest
    do while (associated(this%first))
      rest => this%first%next
      deallocate(this%first)
      this%first => rest
    end do
  end subroutine

  subroutine scalar_func_vlist_append(this, v)
    class(scalar_func_vlist) :: this
    type(scalar_func_box), allocatable :: v(:)
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
    call move_alloc(v, last%v)
    this%n = this%n + 1
  end subroutine scalar_func_vlist_append

  subroutine scalar_func_vlist_to_box_array(vlist, array)
    type(scalar_func_vlist) :: vlist
    type(scalar_func_box), allocatable, intent(out) :: array(:,:)
    integer :: i, j
    type(vlist_func), pointer :: rest
    allocate(array(vlist%m,vlist%n))
    do j = 1, vlist%n
      ASSERT(associated(vlist%first))
      do i = 1, vlist%m
        call move_alloc(vlist%first%v(i)%f, array(i,j)%f)
      end do
      rest => vlist%first%next
      deallocate(vlist%first)
      vlist%first => rest
    end do
    ASSERT(.not.associated(vlist%first))
    vlist%m = 0
    vlist%n = 0
  end subroutine scalar_func_vlist_to_box_array

end module scalar_func_containers
