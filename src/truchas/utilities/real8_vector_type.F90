!!
!! REAL8_VECTOR
!!
!! This module defines a real8 vector class in the spirit of C++'s standard vector.
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

module real8_vector_type
  use,intrinsic :: iso_fortran_env, only: r8 => real64  
  implicit none
  private

  type, public :: real8_vector
     private
     real(r8), allocatable, private :: data_m(:)
     integer, private :: size_m = 0
   contains
     procedure :: at => real8_vector_at
     procedure :: set => real8_vector_set
     procedure :: size => real8_vector_size
     procedure :: resize => real8_vector_resize
     procedure :: capacity => real8_vector_capacity
     procedure :: reserve => real8_vector_reserve
     procedure :: shrink_to_fit => real8_vector_shrink_to_fit
     procedure :: push_back => real8_vector_push_back
     procedure :: pop_back => real8_vector_pop_back
     procedure, private :: enlarge_allocation => real8_vector_type_enlarge_allocation
  end type real8_vector

  public :: assignment(=)
  interface assignment(=)
     module procedure real8_vector_copy_assignment
  end interface
  
contains

  function real8_vector_at(this, a_index) result(a_value)

    class(real8_vector), intent(in) :: this
    integer, intent(in) :: a_index
    real(r8) :: a_value

    ASSERT(a_index <= this%size())
    a_value = this%data_m(a_index)
    return

  end function real8_vector_at

  subroutine real8_vector_set(this, a_index, a_value)

    class(real8_vector), intent(inout) :: this
    integer, intent(in) :: a_index
    real(r8), intent(in) :: a_value

    ASSERT(a_index <= this%size())
    this%data_m(a_index) = a_value

  end subroutine real8_vector_set

  subroutine real8_vector_push_back(this, a_value)

    class(real8_vector), intent(inout) :: this
    real(r8), intent(in) :: a_value

    this%size_m = this%size_m + 1
    if(this%size() > this%capacity()) then
      call this%enlarge_allocation(this%size())
    end if
    this%data_m(this%size_m) = a_value
    
  end subroutine real8_vector_push_back

  subroutine real8_vector_pop_back(this)
    class(real8_vector), intent(inout) :: this

    ASSERT(this%size() > 0)
    this%size_m = this%size_m - 1
    
  end subroutine real8_vector_pop_back

  function real8_vector_size(this) result(a_size)
  
    class(real8_vector), intent(in) :: this
    integer :: a_size
    
    a_size = this%size_m
    return
  end function real8_vector_size

  subroutine real8_vector_resize(this, a_new_size)
    
    class(real8_vector), intent(inout) :: this
    integer, intent(in) :: a_new_size

    ASSERT(a_new_size >= 0)    
    if(a_new_size > this%capacity()) then
      call this%enlarge_allocation(a_new_size)
   end if
   this%size_m = a_new_size

  end subroutine real8_vector_resize

  function real8_vector_capacity(this) result(a_capacity)  
    class(real8_vector), intent(in) :: this
    integer :: a_capacity
    if(allocated(this%data_m)) then      
      a_capacity = size(this%data_m,1)
    else
      a_capacity = 0
    end if
    return
  end function real8_vector_capacity

  subroutine real8_vector_reserve(this, a_new_capacity)    
    class(real8_vector), intent(inout) :: this
    integer,intent(in) :: a_new_capacity

    real(r8), allocatable :: tmp(:)

    ASSERT(a_new_capacity >= 0)

    if(a_new_capacity > this%capacity()) then
       ASSERT(a_new_capacity >= this%size())
       if(allocated(this%data_m)) then
          allocate(tmp(a_new_capacity))
          if(this%size() > 0) then
             tmp(1:this%size()) = this%data_m(this%size())
          end if
          call move_alloc(tmp, this%data_m)
       else
          if(allocated(this%data_m)) then
             deallocate(this%data_m)
          end if
          allocate(this%data_m(a_new_capacity))
       end if
    end if
   
  end subroutine real8_vector_reserve

  subroutine real8_vector_shrink_to_fit(this)    
    class(real8_vector), intent(inout) :: this

    real(r8), allocatable :: tmp(:)

    if(.not.(allocated(this%data_m))) then
      ASSERT(this%size() == 0)
      return
    end if
    allocate(tmp(this%size()))    
    tmp = this%data_m(1:this%size())
    call move_alloc(tmp, this%data_m)
    
  end subroutine real8_vector_shrink_to_fit

  subroutine real8_vector_type_enlarge_allocation(this, a_min_size)
    class(real8_vector), intent(inout) :: this
    integer, intent(in) :: a_min_size

    integer :: heuristic_growth

    heuristic_growth = ceiling(2.0*real(this%capacity(), r8))

    if(heuristic_growth >= a_min_size) then
      call this%reserve(heuristic_growth)
    else
      call this%reserve(a_min_size)
    end if
   
    
  end subroutine real8_vector_type_enlarge_allocation

  subroutine real8_vector_copy_assignment(this, a_other)
    type(real8_vector), intent(inout) :: this
    type(real8_vector), intent(in) :: a_other

    if(a_other%size() == 0) then
      this%size_m = 0
      return
    end if

    ASSERT(allocated(a_other%data_m))
    call this%resize(a_other%size())
    ASSERT(allocated(this%data_m))   
    this%data_m(1:a_other%size()) = a_other%data_m(1:a_other%size())    
  end subroutine real8_vector_copy_assignment
  
end module real8_vector_type
