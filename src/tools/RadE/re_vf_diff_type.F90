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
    integer :: nface_tot = 0        ! total number of patches (number of matrix columns)
    type(dist_vf) :: A, B           ! view factor matrices, D=B-A
    type(re_patch) :: epA, epB      ! enclosure patches for matrices A and B
    real(r8), allocatable :: wA(:)  ! ratio of face area to patch area for matrix A
    real(r8), allocatable :: wB(:)  ! ratio of face area to patch area for matrix B
    real(r8), pointer :: areaA(:)   ! area of patches for matrix A
    real(r8), pointer :: areaB(:)   ! area of patches for matrix B
    real(r8), allocatable :: farea(:)  ! area of enclosure faces
  contains
    procedure, public :: init => vf_diff_init
    procedure, public :: compute_norms
  end type vf_diff

contains

  !! Initializes VF_DIFF by reading two enclosure radiation data sets.
  subroutine vf_diff_init(this, dataset1, dataset2)

    use cell_geometry, only: face_normal, vector_length
    use rad_encl_file_type

    class(vf_diff), target, intent(out) :: this
    character(:), allocatable, intent(in) :: dataset1, dataset2

    type(rad_encl_file) :: file
    integer,  allocatable :: fnode(:), fsize(:)
    real(r8), allocatable :: x(:,:)
    integer :: i, nface, npatch
    integer :: xface(2) ! range for face nodes

    !! Get face areas
    if (scl_rank() == 1) then
      call file%open_ro(dataset1)
      call file%get_patch_dims(nface, npatch)
      allocate(this%farea(nface))
      this%nface_tot = nface

      if ((.not. file%has_patches()) .and. file%has_area()) then
        !! Face areas are stored in file, read them.
        call file%get_area(this%farea)
      else
        !! Compute face areas
        call file%get_surface(fsize, fnode, x)
        nface = size(fsize)
        xface = 0
        do i = 1, nface
          xface(1) = xface(2) + 1
          xface(2) = xface(2) + fsize(i)
          associate(face_nodes => fnode(xface(1):xface(2)))
            this%farea(i) = vector_length(face_normal(x(:,face_nodes)))
          end associate
        end do
      end if
    else
      allocate(this%farea(0))
    end if
    call scl_bcast(this%nface_tot)

    call init_VF_fields(this%A, this%epA, this%wA, this%areaA, this%farea, dataset1)
    ASSERT(this%nface_tot == this%epA%nface)
    call init_VF_fields(this%B, this%epB, this%wB, this%areaB, this%farea, dataset2)
    ASSERT(this%nface_tot == this%epB%nface)

  end subroutine vf_diff_init


  !! Initializes all the fields of VF_DIFF associated with a particular VF matrix
  subroutine init_VF_fields(dvf, ep, w, parea, farea, path)

    class(dist_vf), target, intent(out) :: dvf
    class(re_patch), intent(out) :: ep
    real(r8), allocatable, intent(out) :: w(:)
    real(r8), pointer, intent(out) :: parea(:)
    real(r8), target, intent(in) :: farea(:)
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
      ASSERT(all(farea > 0.0_r8))

      !! Get patch areas
      if (.not. ep%has_patches) then
        parea => farea
      else
        !! Patch areas MUST have been written.
        INSIST(dvf%has_area)
        parea => dvf%area
      end if

      ASSERT(all(parea > 0.0_r8))

      !! Compute face weights
      allocate(w(ep%nface))
      if (.not. ep%has_patches) then
        w = 1.0_r8
      else
        do i = 1, ep%nface
          !! We use the patch areas computed by Chaparral. These should almost exactly match
          !! the sum of the face areas, but we don't check that explicitly.
          w(i) = farea(i) / parea(ep%f2p_map(i))
        end do
      end if

      ASSERT(all(w > 0.0_r8))
      ASSERT(all(w <= 1.0_r8 + epsilon(0.0_r8)))
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
        frowA = frowA * this%wA
        frowB = frowB * this%wB
        rowD = abs(frowB - frowA)

        !! Compute face-based columns using reciprocity
        call row_to_col(pvalA, pA, this%areaA)
        call row_to_col(pvalB, pB, this%areaB)
        call this%epA%patch_to_face_array(pvalA, fcolA)
        call this%epB%patch_to_face_array(pvalB, fcolB)
        ASSERT(size(fcolA) == size(fcolB))
        fcolA = fcolA * this%wA(i)
        fcolB = fcolB * this%wB(i)
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
