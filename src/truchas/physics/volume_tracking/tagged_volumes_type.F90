!!
!! Simple class to store the different material fluxes and where they originated from.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module tagged_volumes_type
  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, public :: tagged_volumes
     integer, private :: capacity
     integer, private :: size
     integer, allocatable, private :: index_id(:)
     real(r8), allocatable, private :: flux_volumes(:,:)
   contains
     procedure, public :: init
     final :: delete     
     procedure, public :: clear
     procedure, public :: add_flux     
     procedure, public :: number_of_fluxes
     procedure, public :: get_flux_id
     procedure, public :: get_flux
   
  end type tagged_volumes

contains
  subroutine init(this, a_number_of_phases, a_capacity)

    class(tagged_volumes), intent(inout) :: this
    integer, intent(in) :: a_number_of_phases
    integer, intent(in) :: a_capacity

    this%capacity = a_capacity
    allocate(this%index_id(this%capacity))
    allocate(this%flux_volumes(a_number_of_phases, a_capacity))
    this%size = 0
  end subroutine init

  subroutine delete(this)

    type(tagged_volumes), intent(inout) :: this

    deallocate(this%index_id, this%flux_volumes)

  end subroutine delete

  subroutine clear(this)

    class(tagged_volumes), intent(inout) :: this

    this%size = 0
    this%index_id = -1
    this%flux_volumes = 0.0_r8

  end subroutine clear

  subroutine add_flux(this, a_tag_id, a_nmat_flux)

    class(tagged_volumes), intent(inout) :: this
    integer, intent(in) :: a_tag_id
    real(r8), intent(in) :: a_nmat_flux(:)

    ASSERT(this%size < this%capacity)
    this%size = this%size + 1
    this%index_id(this%size) = a_tag_id
    this%flux_volumes(:,this%size) = a_nmat_flux

  end subroutine add_flux

  function number_of_fluxes(this) result(a_number_of_fluxes)
    
    class(tagged_volumes), intent(in) :: this
    integer :: a_number_of_fluxes

    a_number_of_fluxes = this%size
    return
  end function number_of_fluxes

  function get_flux_id(this, a_index) result(a_id)

    class(tagged_volumes), intent(in) :: this
    integer, intent(in) :: a_index
    integer :: a_id

    ASSERT(a_index < this%size)
    a_id = this%index_id(a_index)
    return
  end function get_flux_id
  
  function get_flux(this, a_index) result(a_flux)

    class(tagged_volumes), intent(in) :: this
    integer, intent(in) :: a_index
    real(r8) :: a_flux(1:size(this%flux_volumes,1))

    ASSERT(a_index < this%size)    
    a_flux = this%flux_volumes(:,a_index)
    return
  end function get_flux
  
  
end module tagged_volumes_type

