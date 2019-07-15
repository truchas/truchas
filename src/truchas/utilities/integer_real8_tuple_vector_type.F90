!!
!! INTEGER_REAL8_TUPLE_VECTOR
!!
!! This module combines the integer_vector and real8_vector classes. The size of each vector
!! will remain the same during use, allowing the integer at an index and the real8 at an
!! index to be treated as a tuple.
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

module integer_real8_tuple_vector_type
  use,intrinsic :: iso_fortran_env, only: r8 => real64  
  use integer_vector_type
  use real8_vector_type
  implicit none
  private

  type, public :: integer_real8_tuple_vector
     private
     type(integer_vector) :: vec_int_m
     type(real8_vector) :: vec_real8_m
   contains
     procedure :: at_int => integer_real8_tuple_vector_at_int
     procedure :: at_r8 => integer_real8_tuple_vector_at_r8
     procedure :: set => integer_real8_tuple_vector_set
     procedure :: size => integer_real8_tuple_vector_size
     procedure :: resize => integer_real8_tuple_vector_resize
     procedure :: capacity => integer_real8_tuple_vector_capacity
     procedure :: reserve => integer_real8_tuple_vector_reserve
     procedure :: shrink_to_fit => integer_real8_tuple_vector_shrink_to_fit
     procedure :: push_back => integer_real8_tuple_vector_push_back
     procedure :: pop_back => integer_real8_tuple_vector_pop_back
  end type integer_real8_tuple_vector

  public :: assignment(=)
  interface assignment(=)
     module procedure integer_real8_tuple_vector_copy_assignment
  end interface
  
contains

  function integer_real8_tuple_vector_at_int(this, a_index) result(a_value)

    class(integer_real8_tuple_vector), intent(in) :: this
    integer, intent(in) :: a_index
    integer :: a_value

    ASSERT(a_index <= this%size())
    a_value = this%vec_int_m%at(a_index)
    return

  end function integer_real8_tuple_vector_at_int

  function integer_real8_tuple_vector_at_r8(this, a_index) result(a_value)

    class(integer_real8_tuple_vector), intent(in) :: this
    integer, intent(in) :: a_index
    real(r8) :: a_value

    ASSERT(a_index <= this%size())
    a_value = this%vec_real8_m%at(a_index)
    return

  end function integer_real8_tuple_vector_at_r8
  
  subroutine integer_real8_tuple_vector_set(this, a_index, a_int_value, a_r8_value)

    class(integer_real8_tuple_vector), intent(inout) :: this
    integer, intent(in) :: a_index
    integer, intent(in) :: a_int_value
    real(r8), intent(in) :: a_r8_value

    ASSERT(a_index <= this%size())
    call this%vec_int_m%set(a_index, a_int_value)
    call this%vec_real8_m%set(a_index, a_r8_value)    

  end subroutine integer_real8_tuple_vector_set

  subroutine integer_real8_tuple_vector_push_back(this, a_int_value, a_r8_value)

    class(integer_real8_tuple_vector), intent(inout) :: this
    integer, intent(in) :: a_int_value
    real(r8), intent(in) :: a_r8_value    

    call this%vec_int_m%push_back(a_int_value)
    call this%vec_real8_m%push_back(a_r8_value)    
    
  end subroutine integer_real8_tuple_vector_push_back

  subroutine integer_real8_tuple_vector_pop_back(this)
    class(integer_real8_tuple_vector), intent(inout) :: this

    call this%vec_int_m%pop_back()
    call this%vec_real8_m%pop_back()    
    
  end subroutine integer_real8_tuple_vector_pop_back

  function integer_real8_tuple_vector_size(this) result(a_size)
  
    class(integer_real8_tuple_vector), intent(in) :: this
    integer :: a_size

    ASSERT(this%vec_int_m%size() == this%vec_real8_m%size())
    a_size = this%vec_int_m%size()
    return
  end function integer_real8_tuple_vector_size

  subroutine integer_real8_tuple_vector_resize(this, a_new_size)
    
    class(integer_real8_tuple_vector), intent(inout) :: this
    integer, intent(in) :: a_new_size

    call this%vec_int_m%resize(a_new_size)
    call this%vec_real8_m%resize(a_new_size)    

  end subroutine integer_real8_tuple_vector_resize

  function integer_real8_tuple_vector_capacity(this) result(a_capacity)  
    class(integer_real8_tuple_vector), intent(in) :: this
    integer :: a_capacity

    ASSERT(this%vec_int_m%capacity() == this%vec_real8_m%capacity())
    a_capacity = this%vec_int_m%capacity()
    return
  end function integer_real8_tuple_vector_capacity

  subroutine integer_real8_tuple_vector_reserve(this, a_new_capacity)    
    class(integer_real8_tuple_vector), intent(inout) :: this
    integer,intent(in) :: a_new_capacity

    call this%vec_int_m%reserve(a_new_capacity)
    call this%vec_real8_m%reserve(a_new_capacity)    
   
  end subroutine integer_real8_tuple_vector_reserve

  subroutine integer_real8_tuple_vector_shrink_to_fit(this)    
    class(integer_real8_tuple_vector), intent(inout) :: this

    call this%vec_int_m%shrink_to_fit()
    call this%vec_real8_m%shrink_to_fit()    
    
  end subroutine integer_real8_tuple_vector_shrink_to_fit

  subroutine integer_real8_tuple_vector_copy_assignment(this, a_other)
    type(integer_real8_tuple_vector), intent(inout) :: this
    type(integer_real8_tuple_vector), intent(in) :: a_other
    
    this%vec_int_m = a_other%vec_int_m
    this%vec_real8_m = a_other%vec_real8_m
    
  end subroutine integer_real8_tuple_vector_copy_assignment
  
end module integer_real8_tuple_vector_type
