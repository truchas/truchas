!!
!! THERMAL_BC_FACTORY_CLASS
!!
!! An abstract factory class that defines the interface for creating abstract
!! boundary condition objects for the heat conduction model.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! November 2018
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module thermal_bc_factory_class

  use bndry_vfunc_class
  use bndry_func1_class
  use bndry_func2_class
  use intfc_func2_class
  implicit none
  private

  type, abstract, public :: thermal_bc_factory
  contains
    procedure(alloc_bf1), deferred :: alloc_dir_bc
    procedure(alloc_bf1), deferred :: alloc_flux_bc
    procedure(alloc_vbf), deferred :: alloc_vflux_bc
    procedure(alloc_bf2), deferred :: alloc_htc_bc
    procedure(alloc_bf2), deferred :: alloc_rad_bc
    procedure(alloc_if2), deferred :: alloc_htc_ic
    procedure(alloc_if2), deferred :: alloc_rad_ic
  end type

  abstract interface
    subroutine alloc_vbf(this, bc, stat, errmsg)
      import
      class(thermal_bc_factory), intent(inout) :: this    !TODO: intent(in)?
      class(bndry_func2), allocatable, intent(out) :: bc
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
    end subroutine
    subroutine alloc_bf1(this, bc, stat, errmsg)
      import
      class(thermal_bc_factory), intent(inout) :: this    !TODO: intent(in)?
      class(bndry_func1), allocatable, intent(out) :: bc
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
    end subroutine
    subroutine alloc_bf2(this, bc, stat, errmsg)
      import
      class(thermal_bc_factory), intent(inout) :: this    !TODO: intent(in)?
      class(bndry_func2), allocatable, intent(out) :: bc
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
    end subroutine
    subroutine alloc_if2(this, ic, stat, errmsg)
      import
      class(thermal_bc_factory), intent(inout) :: this    !TODO: intent(in)?
      class(intfc_func2), allocatable, intent(out) :: ic
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
    end subroutine
  end interface

end module thermal_bc_factory_class
