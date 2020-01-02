!!
!! VSA_PATCH_TYPE
!!
!! This module provides a derived type that encapsulates the patch objects used
!! by the Variational Shape Approximation (VSA) algorithm.
!!
!! David Neill-Asanza <dhna@lanl.gov>
!! 8 May 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  INIT (THIS, FACEID, AREA, CENTER, NORMAL) initializes THIS with a single
!!    face.  FACEID is the global index of the face. AREA, CENTER, and NORMAL
!!    are the corresponding 3D geometrical properties of the face, stored
!!    respectively as a real scalar, and two length 3 rank-1 real arrays.
!!
!!  RESET_PATCH (THIS, FACEID, AREA, WEIGHT) resets THIS to a single face. The
!!    patch proxy is NOT modified.  FACEID is the global index of the face.
!!    AREA is the area of the face.  WEIGHT is a real scalar representing the
!!    weight of the face relative to the patch.
!!
!!  REPLACE (THIS, OTHER) replaces THIS with OTHER, a TYPE(VSA_PATCH) object.
!!    The dynamic memory of OTHER is safely deallocated.
!!
!!  ADD_FACE(THIS, FACEID, AREA, WEIGHT) adds a face to THIS.  FACEID is the
!!    global index of the face.  AREA is a real scalar representing the area of
!!    the face.  WEIGHT is a real scalar representing the weight of the face
!!    relative to the patch.
!!
!!  REMOVE_FACE (THIS, FACEID, AREA) removes a face from THIS.  FACEID is the
!!    global index of the face.  AREA is a real scalar representing the area of
!!    the face.
!!
!!  GET_WEIGHT (THIS, AREA, CENTER, NORMAL) computes the weight of a face
!!    relative to THIS.  AREA, CENTER, and NORMAL are the corresponding 3D
!!    geometrical properties of the face, stored respectively as a real scalar,
!!    and two length 3 rank-1 real arrays.
!!


#include "f90_assert.fpp"

module vsa_patch_type

  use kinds, only: r8
  implicit none

  !! Coefficient for Voronoi distance bias
  real(r8), parameter :: DIST_COEFF = 5_r8

  !! Minimum size of face data arrays
  integer, parameter :: MIN_PATCH_CAPACITY = 8

  type, public :: vsa_patch
    real(r8) :: center(3), normal(3)
    real(r8) :: area
    real(r8) :: total_weight
    real(r8) :: normal_mag
    integer  :: nface
    integer, allocatable  :: face(:)
    real(r8), allocatable :: weight(:)

    contains
      procedure, public :: init => init_patch
      procedure, public :: reset => reset_patch
      procedure, public :: replace
      procedure, public :: add_face
      procedure, public :: remove_face
      procedure, public :: get_weight
  end type

contains


  !! Initialize patch data
  subroutine init_patch(this, faceid, area, center, normal)

    class(vsa_patch), intent(inout) :: this
    real(r8), intent(in) :: area, center(3), normal(3)
    integer, intent(in) :: faceid

    if (.not. allocated(this%face))   allocate(this%face(MIN_PATCH_CAPACITY))
    if (.not. allocated(this%weight)) allocate(this%weight(MIN_PATCH_CAPACITY))

    this%nface        = 1
    this%area         = area
    this%center       = center
    this%normal       = normal
    this%normal_mag   = 1_r8    ! Faces have unit normals
    this%face(1)      = faceid
    this%weight(1)    = 0.0_r8
    this%total_weight = 0.0_r8

  end subroutine init_patch


  !! Reset patch to a single face
  subroutine reset_patch(this, faceid, area, weight)

    class(vsa_patch), intent(inout) :: this
    integer, intent(in) :: faceid
    real(r8), intent(in) :: weight, area

    this%nface     = 1
    this%area      = area
    this%face(1)   = faceid
    this%weight(1) = weight
    this%total_weight = weight

  end subroutine reset_patch


  !! Add a face to the patch
  subroutine add_face(this, faceid, area, weight)

    class(vsa_patch), intent(inout) :: this
    integer, intent(in) :: faceid
    real(r8), intent(in) :: area, weight

    integer, allocatable :: ftemp(:)
    real(r8), allocatable :: wtemp(:)
    integer :: n

    this%nface = this%nface + 1
    this%area = this%area + area
    this%total_weight = this%total_weight + weight

    !! If capacity is exceeded, grow data arrays
    if (this%nface > size(this%face)) then
      n = size(this%face)
      allocate(ftemp(2*n), wtemp(2*n))

      !! Move face array
      ftemp(1:n) = this%face
      call move_alloc(ftemp, this%face)

      !! Move weight array
      wtemp(1:n) = this%weight
      call move_alloc(wtemp, this%weight)
    end if

    this%face(this%nface) = faceid
    this%weight(this%nface) = weight

  end subroutine add_face


  !! Removes the face at the given index
  subroutine remove_face(this, faceidx, area)

    class(vsa_patch), intent(inout) :: this
    integer, intent(in) :: faceidx
    real(r8), intent(in) :: area

    ASSERT(faceidx > 0)
    ASSERT(faceidx <= this%nface)

    this%area = this%area - area
    this%total_weight = this%total_weight - this%weight(faceidx)

    this%face(faceidx) = this%face(this%nface)
    this%weight(faceidx) = this%weight(this%nface)
    this%nface = this%nface - 1

  end subroutine remove_face


  !! Calculate weight of adding a face to the patch
  function get_weight(this, center, normal, radius) result(ret)

    class(vsa_patch), intent(in) :: this
    real(r8), intent(in) :: center(3), normal(3), radius
    real(r8) :: ret

    !! L^{2,1} error
    ret = sum((normal - this%normal)**2)

    !! L^2 error
    ret = ret + sum((center - this%center)**2) / radius**2

  end function get_weight


  !! Replaces this patch with other, taking care to deallocate any memory
  subroutine replace(this, other)

    class(vsa_patch), intent(inout) :: this, other

    call move_alloc(other%face, this%face)
    call move_alloc(other%weight, this%weight)
    this%nface        = other%nface
    this%area         = other%area
    this%center       = other%center
    this%normal       = other%normal
    this%normal_mag   = other%normal_mag
    this%total_weight = other%total_weight

  end subroutine replace


end module vsa_patch_type
