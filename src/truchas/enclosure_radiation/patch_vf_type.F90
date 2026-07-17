!!
!! PATCH_VF_TYPE
!!
!! A concrete implementation for the abstract base class STATIC_VF.
!! This implementation operates on patch-based view factor data, with
!! clusters of enclosure mesh faces comprising patches.
!!
!! David Neill-Asanza <dhna@lanl.gov>, July 2019
!! Neil N. Carlson <nnc@lanl.gov>, refactoring June 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module patch_vf_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use static_vf_class
  use vf_data_type
  use parallel_communication
  implicit none
  private

  type, extends(static_vf), public :: patch_vf
    private
    type(vf_data) :: vf
    real(r8), allocatable :: w(:)
    integer,  allocatable :: f2p_map(:)
    ! In parent classes:
    ! integer, public :: nface, nface_tot
    ! integer, allocatable :: face_map(:)
    ! real, allocatable, public :: amb_vf(:)
    ! logical :: has_ambient
  contains
    procedure :: init
    procedure :: matvec
  end type

contains

  subroutine init(this, file)

    use rad_encl_file_type
    use permutations

    class(patch_vf), intent(out) :: this
    type(rad_encl_file), intent(in) :: file  ! only referenced on IOP

    integer :: n, offset, bsize(nPE)
    integer, allocatable :: f2p_map(:), face_map(:), color_p(:), color(:)
    real(r8), allocatable :: w(:)

    call this%vf%init(file)

    if (is_IOP) call file%get_patch_dims(this%nface_tot, n)
    call broadcast(this%nface_tot)

    !! Get global face to patch map
    allocate(f2p_map(merge(this%nface_tot,0,is_IOP)))
    if (is_IOP) call file%get_f2p_map(f2p_map)

    !! Coloring of the enclosure patches from the block-row matrix partitioning
    allocate(color_p(merge(this%vf%npatch_tot,0,is_IOP)))
    call gather(this%vf%npatch, bsize)
    if (is_IOP) then
      offset = 0
      do n = 1, nPE
        color_p(offset+1:offset+bsize(n)) = n
        offset = offset + bsize(n)
      end do
    end if

    !! Subordinate coloring of the enclosure faces
    allocate(color(merge(this%nface_tot,0,is_IOP)))
    if (is_IOP) color = color_p(f2p_map)

    !! Subordinate face renumbering (new-to-old face map)
    allocate(face_map(merge(this%nface_tot,0,is_IOP)))
    if (is_IOP) call blocked_coloring_map(color, bsize, face_map)
    call scatter(bsize, this%nface)
    allocate(this%face_map(this%nface))
    call scatter(face_map, this%face_map)

    !! Face-to-patch map (reorder, distribute)
    if (is_IOP) call reorder(f2p_map, face_map)
    allocate(this%f2p_map(this%nface))
    call scatter(f2p_map, this%f2p_map)
    deallocate(f2p_map)

    !! Convert global patch indices to local
    call gather(this%vf%npatch, bsize)
    call broadcast(bsize)
    offset = sum(bsize(1:this_PE-1))
    this%f2p_map = this%f2p_map - offset
    ASSERT(all(1 <= this%f2p_map))
    ASSERT(all(this%f2p_map <= this%vf%npatch))

    !! Face weights (read, reorder, distribute)
    allocate(w(merge(this%nface_tot,0,is_IOP)))
    if (is_IOP) then
      call file%get_face_weight(w)
      call reorder(w, face_map)
    end if
    allocate(this%w(this%nface))
    call scatter(w, this%w)
    deallocate(w)

    !! Expand the patch-based ambient view factors into the face-based version
    this%has_ambient = allocated(this%vf%amb_vf)
    allocate(this%amb_vf(this%nface))
    if (this%has_ambient) then
      this%amb_vf(:) = this%vf%amb_vf(this%f2p_map)
    else
      this%amb_vf = 0
    end if

  end subroutine init

  !! Computes the matrix-vector product Y = PHI*X.
  subroutine matvec(this, x, y)

    class(patch_vf), intent(in) :: this
    real(r8), intent(in)  :: x(:)
    real(r8), intent(out) :: y(:)

    integer :: j
    real(r8) :: xp(this%vf%npatch), yp(this%vf%npatch)

    ASSERT(size(x) == this%nface)
    ASSERT(size(y) == this%nface)

    !! Restrict to patch-based space
    xp = 0
    do j = 1, this%nface
      xp(this%f2p_map(j)) = xp(this%f2p_map(j)) + this%w(j)*x(j)
    end do

    call this%vf%matvec(xp, yp)

    !! Interpolate to face-based space
    do j = 1, this%nface
      y(j) = yp(this%f2p_map(j))
    end do

  end subroutine matvec

  !! Given a coloring vector COLOR, this auxiliary routine computes a mapping
  !! (or renumbering) that makes the coloring a blocked coloring, and the sizes
  !! of the resulting blocks.  That is, the array COLOR(MAP(:)) is such that
  !! the first BSIZE(1) elements are color 1, the next BSIZE(2) elements are
  !! color 2, and so on.  The mapping preserves the relative order of elements
  !! having the same color.

  subroutine blocked_coloring_map(color, bsize, map)

    use permutations

    integer, intent(in)  :: color(:)  ! coloring
    integer, intent(out) :: bsize(:)  ! block size
    integer, intent(out) :: map(:)    ! mapping permutation

    integer :: j, n, next(size(bsize))

    ASSERT(size(color) == size(map))
    ASSERT(minval(color) >= 1 .and. maxval(color) <= size(bsize))

    !! Compute the block size of each color.
    bsize = 0
    do j = 1, size(color)
      bsize(color(j)) = bsize(color(j)) + 1
    end do

    !! NEXT(j) is the next free index for block j.
    next(1) = 1
    do n = 2, size(bsize)
      next(n) = next(n-1) + bsize(n-1)
    end do

    !! Generate the mapping (new-to-old)
    do j = 1, size(color)
      map(next(color(j))) = j
      next(color(j)) = next(color(j)) + 1
    end do
    ASSERT(is_perm(map))

  end subroutine blocked_coloring_map

end module patch_vf_type
