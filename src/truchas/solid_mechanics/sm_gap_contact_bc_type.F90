!!
!! Zach Jibben <zjibben@lanl.gov>
!! October 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module sm_gap_contact_bc_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use bndry_func1_class
  use bndry_face_func_type
  use integration_geometry_type
  implicit none
  private

  type, extends(bndry_func1), public :: sm_gap_contact_bc
    private
    real(r8), public, allocatable :: rotation_matrix(:,:,:)
    logical, public :: enabled

    type(bndry_face_func), allocatable :: bff
    real(r8), allocatable :: area_ip(:)
    integer, allocatable :: fini(:), xfini(:)
  contains
    procedure :: init
    procedure :: compute
  end type sm_gap_contact_bc

contains

  subroutine init(this, mesh, ig, bff)

    use parallel_communication, only: global_any

    class(sm_gap_contact_bc), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(integration_geometry), intent(in) :: ig
    type(bndry_face_func), intent(inout), allocatable :: bff

    !! TODO-WARN
    allocate(this%index(0))

    this%enabled = global_any(size(this%index) > 0)
    
  end subroutine init


  !! TODO-WARN
  subroutine compute(this, t)

    class(sm_gap_contact_bc), intent(inout) :: this
    real(r8), intent(in) :: t

    !! TODO-WARN

  end subroutine compute

end module sm_gap_contact_bc_type
