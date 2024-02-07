!!
!! EM_BC_type
!!
!! This module provides a type encapsulating EM BC data.
!!
!! Zach Jibben <zjibben@lanl.gov>
!! April 2023
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module em_bc_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use,intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use truchas_logging_services
  use simpl_mesh_type
  use scalar_func_class
  use vector_func_class
  use bndry_func1_class
  implicit none
  private

  type, public :: em_bc
    private
    class(bndry_func1), allocatable :: ebc, hbc
    logical, allocatable, public :: is_ebc_edge(:)
    real(r8), allocatable, public :: efield(:)
    real(r8), allocatable, public :: hsource(:)
  contains
    procedure :: init
    procedure :: compute
  end type em_bc

contains

  subroutine init(this, mesh, ebc, hbc)

    use parameter_list_type

    class(em_bc), intent(out) :: this
    type(simpl_mesh), intent(in) :: mesh
    class(bndry_func1), allocatable, intent(inout) :: ebc, hbc

    integer :: j

    call move_alloc(ebc, this%ebc)
    call move_alloc(hbc, this%hbc)

    allocate(this%is_ebc_edge(mesh%nedge), source=.false.)
    if (allocated(this%ebc)) then !TODO: get rid of this mask array
      do j = 1, size(this%ebc%index)
        this%is_ebc_edge(this%ebc%index(j)) = .true.
      end do
    end if

    allocate(this%efield(mesh%nedge), this%hsource(mesh%nedge))
    this%efield = 0
    this%hsource = 0

  end subroutine init


  subroutine compute(this, t)
    class(em_bc), intent(inout) :: this
    real(r8), intent(in) :: t
    integer :: j
    ASSERT(allocated(this%is_ebc_edge))
    if (allocated(this%ebc)) then
      call this%ebc%compute(t)
      do j = 1, size(this%ebc%index)
        this%efield(this%ebc%index(j)) = this%ebc%value(j)
      end do
    end if
    if (allocated(this%hbc)) then
      call this%hbc%compute(t)
      this%hsource = 0.0_r8
      do j = 1, size(this%hbc%index)
        this%hsource(this%hbc%index(j)) = this%hbc%value(j)
      end do
      !where (this%is_ebc_edge) this%hsource = 0 ! already ensured by init of hbc
    end if
    ASSERT(all(ieee_is_finite(this%efield)))
    ASSERT(all(ieee_is_finite(this%hsource)))
  end subroutine compute

end module em_bc_type
