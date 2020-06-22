!!
!! RE_VF_DIFF_TYPE
!!
!! This module provides a derived type for describing the difference of two
!! distributed view factor matrices of a radiation enclosure and methods that
!! operate on instances of this type.
!!
!! David Neill-Asanza <dhna@lanl.gov>
!! 24 Feb 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  The VF_DIFF type measures the difference between two view factor data sets
!!  of the same enclosure.  Specifically, VF_DIFF computes various matrix norms
!!  of the difference matrix D = B - A, where A and B are given view factor
!!  matrices.  The matrix D is always a face-based matrix, and VF_DIFF takes
!!  care of expanding A or B if either is patch-based.
!!
!!  The VF_DIFF type has the following type bound subroutines.
!!
!!  INIT(DATASET1, DATASET2) initializes this VF_DIFF object with the view
!!    factor data read from the radiation enclosure datasets DATASET1 and
!!    DATASET2.  Each dataset must contain view factor data created by the
!!    genre program with the same ENCLOSURE namelist.  The datasets DATASET1
!!    and DATASET2 correspond to matrices A and B, respectively. This is a
!!    collective procedure.
!!
!!  COMPUTE_NORMS(MAX_NORM, ONE_NORM, TWO_NORM, MAX_VAL, RLOC, CLOC) computes
!!    various matrix norms of the VF difference matrix D = B - A.  The matrix D
!!    is not formed explicitly.  This is a collective procedure.
!!
!!    The following matrix norms are stored in the return variables:
!!
!!    MAX_NORM: the max norm of D as an operator.  If e is the irradiance
!!      discrepancy resulting from D and a given radiosity q, then
!!      |e|_max <= |D|_max |q|_max, where |D|_max is maximum absolute row
!!      sum of D.
!!
!!    ONE_NORM: the L1 norm of D as an operator.  If e is the irradiance
!!      discrepancy resulting from D and a given radiosity q, then
!!      |e|_1 <= |D|_1 |q|_1, where |D|_1 is the maximum absolute column
!!      sum of D.
!!
!!    TWO_NORM: the L2 norm of D as an operator.  If e is the irradiance
!!      discrepancy resulting from D and a given radiosity q, then
!!      |e|_2 <= |D|_2 |q|_2, where |D|_2 is the square root of the sum over j
!!      of the dot product of the jth row with the jth column of D.
!!
!!    MAX_VAL, RLOC, CLOC: the element-wise max norm of D, and the location of
!!      the maximum element.  In other words, MAX_VAL is the element of D of
!!      maximum absolute value and its row RLOC and column CLOC location.
!!


#include "f90_assert.fpp"

module re_vf_diff_type

  use kinds, only: r8
  use scl
  use re_patch_type
  use re_dist_vf_type
  implicit none
  private

  type, public :: vf_diff
    integer :: nface_tot = 0    ! total number of patches (number of matrix columns)
    type(dist_vf) :: A, B       ! view factor matrices, D=B-A
    type(re_patch) :: epA, epB  ! enclosure patches for matrices A and B
  contains
    procedure, public :: init
    procedure, public :: compute_norms
  end type vf_diff

contains

  !! Initializes VF_DIFF by reading two enclosure radiation data sets.
  subroutine init(this, dataset1, dataset2)

    class(vf_diff), target, intent(out) :: this
    character(:), allocatable, intent(in) :: dataset1, dataset2

    call init_VF_fields(this%A, this%epA, dataset1)
    call init_VF_fields(this%B, this%epB, dataset2)
    ASSERT(this%epA%nface == this%epB%nface)

    this%nface_tot = this%epA%nface

  end subroutine init


  !! Initializes all the fields of VF_DIFF associated with a particular VF matrix
  subroutine init_VF_fields(dvf, ep, path)

    class(dist_vf), target, intent(out) :: dvf
    class(re_patch), intent(out) :: ep
    character(:), allocatable, intent(in) :: path

    integer :: i

    !! Read view factor data
    call read_dist_vf(dvf, path)

    !! Read patch data
    call ep%read(path)
    call ep%bcast()

    !! Face-to-patch map is not initialized for face-based matrices
    if (.not. ep%has_patches) ep%f2p_map = [(i, i=1,ep%nface)]

    !! All other data is only needed on rank 1
    if (scl_rank() == 1) then

      !! Get patch areas
      if (.not. dvf%has_area) then
        !! Patch areas are always written, so dvf MUST be face-based.
        INSIST(.not. ep%has_patches)
        !! Compute face areas
        block
          use re_utilities, only: compute_face_area
          use rad_encl_file_type
          type(rad_encl_file) :: file
          integer,  allocatable :: fnode(:), fsize(:), xface(:)
          real(r8), allocatable :: x(:,:)
          call file%get_surface(fsize, fnode, x)
          allocate(xface(ep%nface+1))
          xface(1) = 1
          do i = 1, ep%nface
            xface(i+1) = xface(i) + fsize(i)
          end do
          call compute_face_area(xface, fnode, x, dvf%area)
          deallocate(fsize, fnode, xface, x)
        end block
      end if

      ASSERT(all(dvf%area > 0.0_r8))

      !! Get face weights.
      if (.not. dvf%has_weight) then
        dvf%has_weight = .true.
        !! Face weights are read by a patch-based dvf, so it MUST be face-based.
        INSIST(.not. ep%has_patches)
        ASSERT(.not. allocated(dvf%w))
        allocate(dvf%w(ep%nface), source=1.0_r8)
      end if

      ASSERT(all(dvf%w > 0.0_r8))
      ASSERT(all(dvf%w <= 1.0_r8))
    end if

  end subroutine init_VF_fields


  !! Converts VF row n into VF column n using reciprocity
  subroutine row_to_col(val, n, area)

    real, intent(inout) :: val(:)
    integer, intent(in) :: n
    real(r8), intent(in) :: area(:)

    ASSERT(size(val) == size(area))
    ASSERT(n <= size(val))

    val = val * area(n) / area

  end subroutine


  !! Computes the max, one and two norm of D=B-A without forming D explicitly.
  subroutine compute_norms(this, max_norm, one_norm, two_norm, max_val, rloc, cloc)

    class(vf_diff), intent(in) :: this
    real, intent(out) :: max_norm, one_norm, two_norm, max_val
    integer, intent(out) :: rloc, cloc  ! row and column of maximum value

    integer :: i, n, my_rank, pA, pB
    real, allocatable :: pvalA(:), pvalB(:)  ! unpacked column/row of the patch VF matrices
    real, allocatable :: frowA(:), fcolA(:)  ! unpacked column/row of the face VF matrix A
    real, allocatable :: frowB(:), fcolB(:)  ! unpacked column/row of the face VF matrix B
    real, allocatable :: rowD(:), colD(:)    ! unpacked column/row of VF matrix D

    my_rank = scl_rank()

    n = merge(this%nface_tot, 0, my_rank==1)
    allocate(frowA(n), fcolA(n), frowB(n), fcolB(n), rowD(n), colD(n))

    max_norm = -huge(0.0)
    one_norm = -huge(0.0)
    two_norm = 0.0
    max_val = -huge(0.0)
    rloc = 0
    cloc = 0

    !! Compute D = B - A row by row
    do i = 1, this%nface_tot
      !! Get unpacked patch-based rows on rank 1
      pA = this%epA%f2p_map(i)
      pB = this%epB%f2p_map(i)
      pvalA = unpack_dvf_row(this%A, pA)
      pvalB = unpack_dvf_row(this%B, pB)

      if (my_rank == 1) then
        !! Compute face-based rows
        call this%epA%patch_to_face_array(pvalA, frowA)
        call this%epB%patch_to_face_array(pvalB, frowB)
        ASSERT(size(frowA) == size(frowB))
        frowA = frowA * this%A%w
        frowB = frowB * this%B%w
        rowD = abs(frowB - frowA)

        !! Compute face-based columns using reciprocity
        call row_to_col(pvalA, pA, this%A%area)
        call row_to_col(pvalB, pB, this%B%area)
        call this%epA%patch_to_face_array(pvalA, fcolA)
        call this%epB%patch_to_face_array(pvalB, fcolB)
        ASSERT(size(fcolA) == size(fcolB))
        fcolA = fcolA * this%A%w(i)
        fcolB = fcolB * this%B%w(i)
        colD = abs(fcolB - fcolA)

        !! Compute norms
        max_norm = max(max_norm, sum(rowD))
        one_norm = max(one_norm, sum(colD))
        !TODO: this is the Frobenius norm. Change docs?
        two_norm = two_norm + dot_product(rowD, rowD)

        !! Find maximum entry
        n = maxloc(rowD, dim=1)
        if (max_val < rowD(n)) then
          max_val = rowD(n)
          rloc = i
          cloc = n
        end if
      end if
    end do

    two_norm = sqrt(two_norm)

    call scl_bcast(max_norm)
    call scl_bcast(one_norm)
    call scl_bcast(two_norm)
    call scl_bcast(max_val)
    call scl_bcast(rloc)
    call scl_bcast(cloc)

  end subroutine compute_norms

end module re_vf_diff_type
