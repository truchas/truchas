!!
!! VSA_MIN_HEAP
!!
!! This module implements a minimum heap data structure for prioritizing
!! the faces. It is used in the VSA patching algorithm.
!!
!! David Neill-Asanza <dhna@lanl.gov>
!! 2 May 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  The module defines the VSA_MIN_HEAP derived type. It has the
!!  following type bound procedures.
!!
!!  INIT (CAPACITY) initializes the min-heap. CAPACITY is the maximum number of
!!    entries that can be stored in the heap. For the VSA algorithm, CAPACITY
!!    must be at least NEDGE_PER_FACE*NFACE, since a face may be put on the
!!    heap by each of its neighbors.
!!
!!  PUT (FACEID, PATCH, WEIGHT) adds an entry to the heap. FACEID is the ID
!!    of face being added, PATCH is the ID of the patch the face is being
!!    tested against, and WEIGHT is the cost of adding the face to that patch.
!!
!!  POP () returns the element of least weight in O(1) time. The returned
!!    element is also removed from the heap. This function assumes that the heap
!!    is not empty.
!!
!!  PEEK () returns the element of least weight in O(1) time. The heap is not
!!    modified. This function assumes that the heap is not empty.
!!
!!  EMPTY () returns true if the heap is empty, false otherwise.
!!


#include "f90_assert.fpp"

module vsa_min_heap

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, public :: heap_entry
    integer  :: faceid
    integer  :: patch   ! patch this face is being tested against
    real(r8) :: weight  ! cost of adding this face to patch
  end type

  type, public :: min_heap
    private
    integer :: n  ! number of elements in heap
    integer :: capacity
    type(heap_entry), allocatable :: heap(:)

    contains
      procedure, public  :: init
      procedure, public  :: put
      procedure, public  :: pop
      procedure, public  :: peek
      procedure, public  :: empty
      procedure, private :: sift_up
      procedure, private :: sift_down
  end type

contains


  subroutine init(this, capacity)

    class(min_heap), intent(out) :: this
    integer, intent(in) :: capacity

    allocate(this%heap(capacity))
    this%n = 0
    this%capacity = capacity

  end subroutine init


  !! Inserts an element into the heap
  subroutine put(this, faceid, patch, weight)

    class(min_heap), intent(inout) :: this
    integer, intent(in) :: faceid, patch
    real(r8), intent(in) :: weight

    integer :: i

    ASSERT(this%n < this%capacity)

    this%n = this%n + 1
    i = this%n

    this%heap(i)%faceid = faceid
    this%heap(i)%patch  = patch
    this%heap(i)%weight = weight

    call this%sift_up

  end subroutine put


  !! Removes the minimum element from the heap
  function pop(this) result(ret)

    class(min_heap), intent(inout) :: this
    type(heap_entry) :: ret

    ASSERT(this%n > 0)

    ret = this%heap(1)
    this%n = this%n - 1

    !! Avoid sifting if heap is empty
    if (this%n > 0) then
      !! Replace root with last element and sift it down
      this%heap(1) = this%heap(this%n + 1)
      call this%sift_down(1)
    end if

  end function


  !! Returns the minimum element from the heap without removal
  function peek(this) result(ret)

    class(min_heap), intent(in) :: this
    type(heap_entry) :: ret

    ASSERT(this%n > 0)

    ret = this%heap(1)

  end function


  !! Checks whether the heap is empty
  logical function empty(this)

    class(min_heap), intent(in) :: this

    empty = this%n <= 0

  end function


  !! Sifts new value up to maintain heap invariant
  subroutine sift_up(this)

  class(min_heap), intent(inout) :: this

    type(heap_entry) :: temp
    integer :: i, p

    i = this%n  ! current index
    p = i/2     ! index of the parent

    !! Swap with parent until element is smaller
    do while (i > 1)
      !! Gets around evaluation of all operands of .AND. with NAG compiler
      if (this%heap(i)%weight >= this%heap(p)%weight) exit
      temp = this%heap(i)
      this%heap(i) = this%heap(p)
      this%heap(p) = temp

      i = p
      p = i/2
    end do

  end subroutine sift_up


  !! Sifts down the value at given index to maintain the heap invariant
  recursive subroutine sift_down(this, i)

    class(min_heap), intent(inout) :: this
    integer, intent(in) :: i

    type(heap_entry) :: temp
    integer :: l, r, min_idx

    l = 2*i      ! left child
    r = 2*i + 1  ! right child

    min_idx = i
    if (l < this%n) then
      !! Gets around evaluation of all operands of .AND. with NAG compiler
      if (this%heap(l)%weight < this%heap(min_idx)%weight) min_idx = l
    end if
    if (r < this%n) then
      !! Gets around evaluation of all operands of .AND. with NAG compiler
      if (this%heap(r)%weight < this%heap(min_idx)%weight) min_idx = r
    end if

    !! Swap root with minimal index and fix sub-tree
    if (min_idx /= i) then
      temp = this%heap(i)
      this%heap(i) = this%heap(min_idx)
      this%heap(min_idx) = temp

      call this%sift_down(min_idx)
    end if

  end subroutine sift_down


end module vsa_min_heap
