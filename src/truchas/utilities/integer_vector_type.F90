!!
!! INTEGER_VECTOR
!!
!! This module defines an integer vector class in the spirit of C++'s standard vector.
!! It has the ability to perform reallocations if more elements than it's current
!! capacity is requested. Note, be aware of this. Constantly forcing reallocations
!! can be costly for performance.
!!
!! Robert Chiodi  <robertchiodi@lanl.gov>
!! June 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module integer_vector_type
  use,intrinsic :: iso_fortran_env, only: r8 => real64 
  implicit none
  private

  type, public :: integer_vector
     private
     integer, pointer, private :: data_m(:) => NULL()
     integer, private :: size_m = 0
   contains
     procedure :: init => integer_vector_init
     procedure :: at => integer_vector_at
     procedure :: set => integer_vector_set
     procedure :: size => integer_vector_size
     procedure :: resize => integer_vector_resize
     procedure :: capacity => integer_vector_capacity
     procedure :: reserve => integer_vector_reserve
     procedure :: shrink_to_fit => integer_vector_shrink_to_fit
     procedure :: push_back => integer_vector_push_back
     procedure :: pop_back => integer_vector_pop_back
     procedure, private :: enlarge_allocation => integer_vector_type_enlarge_allocation
     final ::  integer_vector_delete
  end type integer_vector

  public :: assignment(=)
  interface assignment(=)
     module procedure integer_vector_copy_assignment
  end interface
  
contains

  impure elemental subroutine integer_vector_init(this)

    class(integer_vector), intent(inout) :: this

    this%data_m => NULL()
    this%size_m = 0

  end subroutine integer_vector_init

  function integer_vector_at(this, a_index) result(a_int)

    class(integer_vector), intent(in) :: this
    integer, intent(in) :: a_index
    integer :: a_int

    ASSERT(a_index <= this%size())
    a_int = this%data_m(a_index)
    return

  end function integer_vector_at

  subroutine integer_vector_set(this, a_index, a_value)

    class(integer_vector), intent(inout) :: this
    integer, intent(in) :: a_index
    integer, intent(in) :: a_value

    ASSERT(a_index <= this%size())
    this%data_m(a_index) = a_value

  end subroutine integer_vector_set

  subroutine integer_vector_push_back(this, a_value)

    class(integer_vector), intent(inout) :: this
    integer, intent(in) :: a_value

    this%size_m = this%size_m + 1
    if(this%size() > this%capacity()) then
      call this%enlarge_allocation(this%size())
    end if
    this%data_m(this%size_m) = a_value
    
  end subroutine integer_vector_push_back

  subroutine integer_vector_pop_back(this)
    class(integer_vector), intent(inout) :: this

    ASSERT(this%size() > 0)
    this%size_m = this%size_m - 1
    
  end subroutine integer_vector_pop_back

  function integer_vector_size(this) result(a_size)
  
    class(integer_vector), intent(in) :: this
    integer :: a_size
    
    a_size = this%size_m
    return
  end function integer_vector_size

  subroutine integer_vector_resize(this, a_new_size)
    
    class(integer_vector), intent(inout) :: this
    integer, intent(in) :: a_new_size

    ASSERT(a_new_size >= 0)    
    if(a_new_size > this%capacity()) then
      call this%enlarge_allocation(a_new_size)
   end if
   this%size_m = a_new_size

  end subroutine integer_vector_resize

  function integer_vector_capacity(this) result(a_capacity)  
    class(integer_vector), intent(in) :: this
    integer :: a_capacity
    if(associated(this%data_m)) then      
      a_capacity = size(this%data_m,1)
    else
      a_capacity = 0
    end if
    return
  end function integer_vector_capacity

  subroutine integer_vector_reserve(this, a_new_capacity)    
    class(integer_vector), intent(inout) :: this
    integer,intent(in) :: a_new_capacity

    integer, pointer :: tmp_ptr(:)

    ASSERT(a_new_capacity >= 0)

    if(a_new_capacity > this%capacity()) then
      ASSERT(a_new_capacity >= this%size())
      if(associated(this%data_m)) then
        tmp_ptr => this%data_m
        nullify(this%data_m)
        allocate(this%data_m(a_new_capacity))
        this%data_m(1:this%size()) = tmp_ptr(1:this%size())
        deallocate(tmp_ptr)
      else
        allocate(this%data_m(a_new_capacity))
      end if
      ASSERT(associated(this%data_m))
    end if
   
  end subroutine integer_vector_reserve

  subroutine integer_vector_shrink_to_fit(this)    
    class(integer_vector), intent(inout) :: this

    integer, pointer :: tmp_ptr(:)

    if(.not.(associated(this%data_m))) then
      ASSERT(this%size() == 0)
      return
    end if
    tmp_ptr => this%data_m
    nullify(this%data_m)
    allocate(this%data_m(this%size()))
    this%data_m = tmp_ptr(1:this%size())
    deallocate(tmp_ptr)
    
  end subroutine integer_vector_shrink_to_fit

  subroutine integer_vector_type_enlarge_allocation(this, a_min_size)
    class(integer_vector), intent(inout) :: this
    integer, intent(in) :: a_min_size

    integer :: heuristic_growth

    heuristic_growth = ceiling(2.0*real(this%capacity(), r8))

    if(heuristic_growth >= a_min_size) then
      call this%reserve(heuristic_growth)
    else
      call this%reserve(a_min_size)
    end if
   
    
  end subroutine integer_vector_type_enlarge_allocation

  subroutine integer_vector_delete(this)
    type(integer_vector), intent(inout) :: this

    if(associated(this%data_m)) then
       deallocate(this%data_m)
       nullify(this%data_m)
    end if
    ASSERT(.not.associated(this%data_m))
    this%size_m = 0

  end subroutine integer_vector_delete

  subroutine integer_vector_copy_assignment(this, a_other)
    type(integer_vector), intent(inout) :: this
    type(integer_vector), intent(in) :: a_other

    if(a_other%size() == 0) then
      this%size_m = 0
      return
    end if

    ASSERT(associated(a_other%data_m))
    call this%resize(a_other%size())
    ASSERT(associated(this%data_m))   
    this%data_m(1:a_other%size()) = a_other%data_m(1:a_other%size())    
  end subroutine integer_vector_copy_assignment
  
end module integer_vector_type
