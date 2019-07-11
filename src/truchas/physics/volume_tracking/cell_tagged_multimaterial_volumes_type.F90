!!
!! Simple class to store the different material fluxes and where they originated from.
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

module cell_tagged_mm_volumes_type
  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use integer_vector_type
  use integer_real8_tuple_vector_type

  implicit none
  private

  type, public :: cell_tagged_mm_volumes
     type(integer_vector), private :: cell_id_m
     type(integer_real8_tuple_vector), pointer, private :: phase_moments_m(:) => null()
   contains
     procedure :: get_number_of_cells => cell_tagged_mm_volumes_get_number_of_cells
     procedure :: get_cell_id => cell_tagged_mm_volumes_get_cell_id
     procedure :: add_cell_fluxes => cell_tagged_mm_volumes_add_cell_fluxes
     procedure :: get_cell_fluxes => cell_tagged_mm_volumes_get_cell_fluxes
     procedure :: reserve => cell_tagged_mm_volumes_reserve
     procedure :: clear => cell_tagged_mm_volumes_clear
     procedure :: negate_in_place => cell_tagged_mm_volumes_negate_in_place
     procedure, private :: phase_moments_reserve => cell_tagged_mm_volumes_phase_moments_reserve
     final :: cell_tagged_mm_volumes_delete
  end type cell_tagged_mm_volumes

  public :: assignment(=)
  interface assignment(=)
     module procedure cell_tagged_mm_volumes_copy_assignment
  end interface

contains

  function cell_tagged_mm_volumes_get_number_of_cells(this) result(a_size)
    class(cell_tagged_mm_volumes), intent(in) :: this
    integer :: a_size
    
    a_size = this%cell_id_m%size()
    return
  end function cell_tagged_mm_volumes_get_number_of_cells

  function cell_tagged_mm_volumes_get_cell_id(this, a_index) result(a_id)
    class(cell_tagged_mm_volumes), intent(inout) :: this
    integer, intent(in) :: a_index
    integer :: a_id

    a_id = this%cell_id_m%at(a_index)
    return
  end function cell_tagged_mm_volumes_get_cell_id

  subroutine cell_tagged_mm_volumes_add_cell_fluxes(this, a_cell_id, a_cell_fluxes)
    class(cell_tagged_mm_volumes), intent(inout) :: this
    integer, intent(in) :: a_cell_id
    type(integer_real8_tuple_vector), intent(in) :: a_cell_fluxes

    ASSERT(associated(this%phase_moments_m))
    call this%cell_id_m%push_back(a_cell_id)
    if(this%cell_id_m%capacity() > size(this%phase_moments_m,1)) then
      ! Reallocation was forced during cell_id push_back,
      ! match this in phase_moments
      call this%phase_moments_reserve()
    end if
    this%phase_moments_m(this%cell_id_m%size()) = a_cell_fluxes

  end subroutine cell_tagged_mm_volumes_add_cell_fluxes
  
  function cell_tagged_mm_volumes_get_cell_fluxes(this, a_index) result(a_phase_volumes_ptr)
    class(cell_tagged_mm_volumes), intent(in) :: this
    integer, intent(in) :: a_index
    type(integer_real8_tuple_vector), pointer :: a_phase_volumes_ptr

    ASSERT(a_index <= this%cell_id_m%size())
    a_phase_volumes_ptr => this%phase_moments_m(a_index)
    return
    
  end function cell_tagged_mm_volumes_get_cell_fluxes

  subroutine cell_tagged_mm_volumes_reserve(this, a_new_size)
    class(cell_tagged_mm_volumes), intent(inout) :: this
    integer, intent(in) :: a_new_size

    call this%cell_id_m%reserve(a_new_size)
    call this%phase_moments_reserve()

  end subroutine cell_tagged_mm_volumes_reserve

  subroutine cell_tagged_mm_volumes_phase_moments_reserve(this)
    class(cell_tagged_mm_volumes), intent(inout) :: this

    type(integer_real8_tuple_vector), pointer :: tmp_ptr(:)

    if(associated(this%phase_moments_m)) then
      if(size(this%phase_moments_m,1) < this%cell_id_m%capacity()) then
        tmp_ptr => this%phase_moments_m
        nullify(this%phase_moments_m)
        allocate(this%phase_moments_m(this%cell_id_m%capacity()))
        this%phase_moments_m(1:this%cell_id_m%size()) = tmp_ptr(1:this%cell_id_m%size())
        deallocate(tmp_ptr)
      end if
    else
      allocate(this%phase_moments_m(this%cell_id_m%capacity()))
    end if
    call this%phase_moments_m%init()       
    
  end subroutine cell_tagged_mm_volumes_phase_moments_reserve

  subroutine cell_tagged_mm_volumes_clear(this)
    class(cell_tagged_mm_volumes), intent(inout) :: this

    call this%cell_id_m%resize(0)

  end subroutine cell_tagged_mm_volumes_clear

  subroutine cell_tagged_mm_volumes_negate_in_place(this)
    class(cell_tagged_mm_volumes), intent(inout) :: this

    integer :: n, i
    type(integer_real8_tuple_vector), pointer :: elem

    do n = 1, this%cell_id_m%size()
      elem => this%phase_moments_m(n)
      do i = 1, elem%size()
        call elem%set(i, elem%at_int(i), -elem%at_r8(i))
      end do
    end do

  end subroutine cell_tagged_mm_volumes_negate_in_place
  
  subroutine cell_tagged_mm_volumes_delete(this)
    type(cell_tagged_mm_volumes), intent(inout) :: this

    if(associated(this%phase_moments_m)) then
      deallocate(this%phase_moments_m)
    end if
    nullify(this%phase_moments_m)
    this%phase_moments_m => NULL()
    
  end subroutine cell_tagged_mm_volumes_delete

  subroutine cell_tagged_mm_volumes_copy_assignment(this, a_other)
    type(cell_tagged_mm_volumes), intent(inout) :: this
    type(cell_tagged_mm_volumes), intent(in) :: a_other

    if(a_other%get_number_of_cells() == 0) then
      call this%clear()
      return
    end if
    
    this%cell_id_m = a_other%cell_id_m
    call this%phase_moments_reserve() ! Just matches capacity of this%cell_id_m
    ASSERT(associated(this%phase_moments_m))
    ASSERT(associated(a_other%phase_moments_m))
    this%phase_moments_m(1:a_other%get_number_of_cells()) = a_other%phase_moments_m(1:a_other%get_number_of_cells())
    
  end subroutine cell_tagged_mm_volumes_copy_assignment

end module cell_tagged_mm_volumes_type

