!! Implementation of INDEX_MAP GATHER Subroutines
!!
!! Copyright 2022 Neil N. Carlson <neil.n.carlson@gmail.com>
!! Use subject to the MIT license: https://opensource.org/licenses/MIT
!!

#include "f90_assert.fpp"

submodule(index_map_type) gather_impl
implicit none
contains

!!!! RANK-1 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine gath_i4_1(this, src, dest)
    class(index_map), intent(in) :: this
    integer(i4), intent(in) :: src(:)
    integer(i4), intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(dest) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(src) >= this%onp_size)
    call MPI_Gatherv(src, this%onp_size, MPI_INTEGER4, &
        dest, this%counts, this%displs, MPI_INTEGER4, this%root, this%comm, ierr)
  end subroutine

  module subroutine gath_i8_1(this, src, dest)
    class(index_map), intent(in) :: this
    integer(i8), intent(in) :: src(:)
    integer(i8), intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(dest) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(src) >= this%onp_size)
    call MPI_Gatherv(src, this%onp_size, MPI_INTEGER8, &
        dest, this%counts, this%displs, MPI_INTEGER8, this%root, this%comm, ierr)
  end subroutine

  module subroutine gath_r4_1(this, src, dest)
    class(index_map), intent(in) :: this
    real(r4), intent(in) :: src(:)
    real(r4), intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(dest) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(src) >= this%onp_size)
    call MPI_Gatherv(src, this%onp_size, MPI_REAL4, &
        dest, this%counts, this%displs, MPI_REAL4, this%root, this%comm, ierr)
  end subroutine

  module subroutine gath_r8_1(this, src, dest)
    class(index_map), intent(in) :: this
    real(r8), intent(in) :: src(:)
    real(r8), intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(dest) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(src) >= this%onp_size)
    call MPI_Gatherv(src, this%onp_size, MPI_REAL8, &
        dest, this%counts, this%displs, MPI_REAL8, this%root, this%comm, ierr)
  end subroutine

  module subroutine gath_dl_1(this, src, dest)
    class(index_map), intent(in) :: this
    logical, intent(in) :: src(:)
    logical, intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(dest) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(src) >= this%onp_size)
    call MPI_Gatherv(src, this%onp_size, MPI_LOGICAL, &
        dest, this%counts, this%displs, MPI_LOGICAL, this%root, this%comm, ierr)
  end subroutine

!!!! RANK-2 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! This auxiliary subroutine performs the gather operation for elements that
  !! are blocks of values of the specified MPI data type. The block size need
  !! only be specified on the root process. Note that the array dummy arguments
  !! are assumed-size, so if the actual arguments are not contiguous, copy-in/
  !! copy-out of contiguous temporaries will occur.

  subroutine gath_aux2(this, block_size, mpi_type, src, dest)

    type(index_map), intent(in) :: this
    integer, intent(inout) :: block_size
    integer, intent(in) :: mpi_type
    type(*), intent(in) :: src(*)
    type(*), intent(inout) :: dest(*)

    interface ! explicit interface needed to pass assumed-type arguments
      subroutine MPI_Gatherv(sendbuf, sendcount, sendtype, &
          recvbuf, recvcounts, displs, recvtype, root, comm, ierr)
        type(*), intent(in) :: sendbuf(*)
        type(*), intent(inout) :: recvbuf(*)
        integer, intent(in) :: sendcount, sendtype, recvcounts(*), displs(*), recvtype, root, comm
        integer, intent(out) :: ierr
      end subroutine
    end interface

    integer :: ierr
    integer :: block_type

    call MPI_Bcast(block_size, 1, MPI_INTEGER, this%root, this%comm, ierr)
    call MPI_Type_contiguous(block_size, mpi_type, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Gatherv(src, this%onp_size, block_type, &
        dest, this%counts, this%displs, block_type, this%root, this%comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine gath_aux2

  module subroutine gath_i4_2(this, src, dest)
    class(index_map), intent(in) :: this
    integer(i4), intent(in) :: src(:,:)
    integer(i4), intent(inout) :: dest(:,:)
    integer :: block_size
    ASSERT(size(src,2) >= this%onp_size)
    ASSERT(size(dest,2) >= merge(this%global_size,0,this%is_root))
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(dest(:,1))
    call gath_aux2(this, block_size, MPI_INTEGER4, src, dest)
  end subroutine

  module subroutine gath_i8_2(this, src, dest)
    class(index_map), intent(in) :: this
    integer(i8), intent(in) :: src(:,:)
    integer(i8), intent(inout) :: dest(:,:)
    integer :: block_size
    ASSERT(size(src,2) >= this%onp_size)
    ASSERT(size(dest,2) >= merge(this%global_size,0,this%is_root))
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(dest(:,1))
    call gath_aux2(this, block_size, MPI_INTEGER8, src, dest)
  end subroutine

  module subroutine gath_r4_2(this, src, dest)
    class(index_map), intent(in) :: this
    real(r4), intent(in) :: src(:,:)
    real(r4), intent(inout) :: dest(:,:)
    integer :: block_size
    ASSERT(size(src,2) >= this%onp_size)
    ASSERT(size(dest,2) >= merge(this%global_size,0,this%is_root))
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(dest(:,1))
    call gath_aux2(this, block_size, MPI_REAL4, src, dest)
  end subroutine

  module subroutine gath_r8_2(this, src, dest)
    class(index_map), intent(in) :: this
    real(r8), intent(in) :: src(:,:)
    real(r8), intent(inout) :: dest(:,:)
    integer :: block_size
    ASSERT(size(src,2) >= this%onp_size)
    ASSERT(size(dest,2) >= merge(this%global_size,0,this%is_root))
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(dest(:,1))
    call gath_aux2(this, block_size, MPI_REAL8, src, dest)
  end subroutine

  module subroutine gath_dl_2(this, src, dest)
    class(index_map), intent(in) :: this
    logical, intent(in) :: src(:,:)
    logical, intent(inout) :: dest(:,:)
    integer :: block_size
    ASSERT(size(src,2) >= this%onp_size)
    ASSERT(size(dest,2) >= merge(this%global_size,0,this%is_root))
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(dest(:,1))
    call gath_aux2(this, block_size, MPI_LOGICAL, src, dest)
  end subroutine

!!!! RANK-3 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine gath_i4_3(this, src, dest)
    class(index_map), intent(in) :: this
    integer(i4), intent(in) :: src(:,:,:)
    integer(i4), intent(inout) :: dest(:,:,:)
    integer :: block_size
    ASSERT(size(src,3) >= this%onp_size)
    ASSERT(size(dest,3) >= merge(this%global_size,0,this%is_root))
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(dest(:,:,1))
    call gath_aux2(this, block_size, MPI_INTEGER4, src, dest)
  end subroutine

  module subroutine gath_i8_3(this, src, dest)
    class(index_map), intent(in) :: this
    integer(i8), intent(in) :: src(:,:,:)
    integer(i8), intent(inout) :: dest(:,:,:)
    integer :: block_size
    ASSERT(size(src,3) >= this%onp_size)
    ASSERT(size(dest,3) >= merge(this%global_size,0,this%is_root))
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(dest(:,:,1))
    call gath_aux2(this, block_size, MPI_INTEGER8, src, dest)
  end subroutine

  module subroutine gath_r4_3(this, src, dest)
    class(index_map), intent(in) :: this
    real(r4), intent(in) :: src(:,:,:)
    real(r4), intent(inout) :: dest(:,:,:)
    integer :: block_size
    ASSERT(size(src,3) >= this%onp_size)
    ASSERT(size(dest,3) >= merge(this%global_size,0,this%is_root))
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(dest(:,:,1))
    call gath_aux2(this, block_size, MPI_REAL4, src, dest)
  end subroutine

  module subroutine gath_r8_3(this, src, dest)
    class(index_map), intent(in) :: this
    real(r8), intent(in) :: src(:,:,:)
    real(r8), intent(inout) :: dest(:,:,:)
    integer :: block_size
    ASSERT(size(src,3) >= this%onp_size)
    ASSERT(size(dest,3) >= merge(this%global_size,0,this%is_root))
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(dest(:,:,1))
    call gath_aux2(this, block_size, MPI_REAL8, src, dest)
  end subroutine

  module subroutine gath_dl_3(this, src, dest)
    class(index_map), intent(in) :: this
    logical, intent(in) :: src(:,:,:)
    logical, intent(inout) :: dest(:,:,:)
    integer :: block_size
    ASSERT(size(src,3) >= this%onp_size)
    ASSERT(size(dest,3) >= merge(this%global_size,0,this%is_root))
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(dest(:,:,1))
    call gath_aux2(this, block_size, MPI_LOGICAL, src, dest)
  end subroutine

end submodule gather_impl
