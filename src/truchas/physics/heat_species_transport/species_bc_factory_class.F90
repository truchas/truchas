!!
!! SPECIES_BC_FACTORY_CLASS
!!
!! An abstract factory class that defines the interface for creating abstract
!! boundary condition objects for the species advection-diffusion model.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! March 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module species_bc_factory_class

  use bndry_func1_class
  use bndry_func2_class
  implicit none
  private

  type, abstract, public :: species_bc_factory
  contains
    procedure(alloc_bf1), deferred :: alloc_dir_bc
    procedure(alloc_bf1), deferred :: alloc_flux_bc
    procedure(alloc_bf2), deferred :: alloc_mtc_bc
  end type

  abstract interface
    subroutine alloc_bf1(this, comp, bc, stat, errmsg)
      import
      class(species_bc_factory), intent(inout) :: this    !TODO: intent(in)?
      integer, intent(in) :: comp
      class(bndry_func1), allocatable, intent(out) :: bc
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
    end subroutine
    subroutine alloc_bf2(this, comp, bc, stat, errmsg)
      import
      class(species_bc_factory), intent(inout) :: this    !TODO: intent(in)?
      integer, intent(in) :: comp
      class(bndry_func2), allocatable, intent(out) :: bc
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
    end subroutine
  end interface

end module species_bc_factory_class
