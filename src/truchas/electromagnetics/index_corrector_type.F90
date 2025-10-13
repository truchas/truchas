!!
!! INDEX_CORRECTOR_TYPE
!!
!! This module provides a helper type for mapping indices for a parallel block-
!! structured matrix with blocks of differing sizes (e.g., an edge*edge block,
!! followd by a edge*node block). See fdme_model_type.F90 for usage.
!!
!! Zach Jibben <zjibben@lanl.gov>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module index_corrector_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use index_map_type
  implicit none
  private

  type, public :: index_corrector
    private
    type(index_map), pointer, public :: imap => null() ! new imap for combined indices
    logical :: imap_owned = .false.

    integer :: total_onp_size
    integer, allocatable :: onp_size(:), offp_size(:), onp_offset(:), offp_offset(:), offp_map(:)
  contains
    final :: delete_index_corrector
    procedure :: init
    procedure :: eval
  end type index_corrector

contains

  elemental subroutine delete_index_corrector(this)
    type(index_corrector), intent(inout) :: this
    if (this%imap_owned .and. associated(this%imap)) deallocate(this%imap)
  end subroutine delete_index_corrector


  subroutine init(this, block_imaps)

    class(index_corrector), intent(out) :: this
    type(index_map), intent(in) :: block_imaps(:)

    integer :: bsize, i, b, nblocks
    integer, allocatable :: gid(:), offp_index_original(:)

    call compute_mixed_offp_index(block_imaps, offp_index_original)
    bsize = sum(block_imaps(:)%onp_size)
    allocate(this%imap)
    call this%imap%init(bsize, offp_index_original)
    this%imap_owned = .true.

    nblocks = size(block_imaps)
    allocate(this%onp_size(nblocks), this%offp_size(nblocks), &
        this%onp_offset(nblocks), this%offp_offset(nblocks))

    do b = 1, nblocks
      this%onp_size(b) = block_imaps(b)%onp_size
      this%offp_size(b) = block_imaps(b)%offp_size
    end do
    this%total_onp_size = sum(this%onp_size)

    do b = 1, nblocks
      this%onp_offset(b) = sum(this%onp_size(:b-1))
      this%offp_offset(b) = sum(this%offp_size(:b-1)) - this%onp_size(b)
    end do

    ! The new imap internally sorts a copy of offp_index after handing it off,
    ! so we need to map to the new sorting.
    allocate(this%offp_map(this%imap%offp_size), gid(this%imap%local_size))
    do i = 1, this%imap%onp_size
      gid(i) = this%imap%global_index(i)
    end do
    call this%imap%gather_offp(gid)
    call compute_map(offp_index_original, this%imap%offp_index, this%offp_map)

  end subroutine init


  function eval(this, i, b) result(corrected_index)
    class(index_corrector), intent(in) :: this
    integer, intent(in) :: i, b
    integer :: corrected_index
    ASSERT(b <= size(this%onp_size))
    INSIST(i > 0)
    INSIST(i <= this%onp_size(b) + this%offp_size(b))
    if (i <= this%onp_size(b)) then
      corrected_index = i + this%onp_offset(b)
    else
      corrected_index = this%total_onp_size + this%offp_map(i + this%offp_offset(b))
    end if
  end function eval


  subroutine compute_mixed_offp_index(block_imaps, offp_index)

    use mpi
    use parallel_communication, only: npe, comm, this_pe

    type(index_map), intent(in) :: block_imaps(:)
    integer, allocatable, intent(out) :: offp_index(:)

    integer :: ierr, offp_size, i1, i2, b
    integer :: offset(npe), first_gid(npe)
    integer :: Nonp(size(block_imaps)), Noffp(size(block_imaps))
    integer :: first_gid_x(npe,size(block_imaps)), n_onp(npe,size(block_imaps))

    if (allocated(offp_index)) deallocate(offp_index)

    do b = 1, size(block_imaps)
      Nonp(b) = block_imaps(b)%onp_size
      Noffp(b) = block_imaps(b)%offp_size
      call MPI_Allgather(block_imaps(b)%first_gid, 1, MPI_INTEGER, &
          first_gid_x(:,b), 1, MPI_INTEGER, comm, ierr)
      ASSERT(ierr == 0)
      call MPI_Allgather(block_imaps(b)%onp_size, 1, MPI_INTEGER, n_onp(:,b), 1, MPI_INTEGER, &
          comm, ierr)
      ASSERT(ierr == 0)
    end do

    offp_size = sum(Noffp)
    allocate(offp_index(offp_size))

    if (offp_size > 0) then
      first_gid = sum(first_gid_x - 1, dim=2) + 1

      i2 = 0
      offset = 0
      do b = 1, size(block_imaps)
        i1 = i2 + 1
        i2 = i1 + Noffp(b) - 1
        associate (nb_offp => offp_index(i1:i2))
          nb_offp = block_imaps(b)%offp_index
          call correct_offp(nb_offp, first_gid_x(:,b), first_gid, offset)
          offset = offset + n_onp(:,b)
        end associate
      end do
    end if

  contains

    subroutine correct_offp(n_offp, first_gid_x, first_gid, offset)

      use parallel_communication, only: this_pe

      integer, intent(inout) :: n_offp(:)
      integer, intent(in) :: first_gid_x(:), first_gid(:), offset(:)

      integer :: i, other_pe

      do i = 1, size(n_offp)
        ! TODO: bisectional search likely faster than findloc
        other_pe = findloc(first_gid_x <= n_offp(i), .true., back=.true., dim=1)
        ASSERT(other_pe > 0 .and. other_pe <= size(first_gid) .and. other_pe /= this_pe)
        n_offp(i) = n_offp(i) - first_gid_x(other_pe) + first_gid(other_pe) + offset(other_pe)
      end do

    end subroutine correct_offp

  end subroutine compute_mixed_offp_index


  ! Compute a map between two given arrays, such that a(:) = b(map(:)).
  ! It is assumed that a and b consist of an identical list of unique
  ! integers, but in different orders.
  !
  ! TODO: This heap sort implementation is O(NlogN), but very cheap on tests
  ! so far. An implementation based on integer_map would be O(logN).
  subroutine compute_map(a, b, map)

    use sort_utilities

    integer, intent(in) :: a(:), b(:)
    integer, intent(out) :: map(:)

    integer :: i
    integer :: xa(size(a)), xb(size(b))

    ASSERT(size(a) == size(b))
    ASSERT(size(a) == size(map))

    call heap_sort(a, xa)
    call heap_sort(b, xb)

    do i = 1, size(a)
      map(xa(i)) = xb(i)
    end do

    ASSERT(all(a == b(map)))

  end subroutine compute_map

end module index_corrector_type
