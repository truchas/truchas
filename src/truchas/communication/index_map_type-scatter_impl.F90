!! Implementation of INDEX_MAP SCATTER Subroutines
!!
!! Copyright 2022 Neil N. Carlson <neil.n.carlson@gmail.com>
!! Use subject to the MIT license: https://opensource.org/licenses/MIT
!!

#include "f90_assert.fpp"

submodule(index_map_type) scatter_impl
implicit none
contains

!!!! RANK-1 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine scat_i4_1(this, src, dest)
    class(index_map), intent(in) :: this
    integer(i4), intent(in) :: src(:)
    integer(i4), intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(src) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(dest) >= this%onp_size)
    call MPI_Scatterv(src, this%counts, this%displs, MPI_INTEGER4, &
        dest, this%onp_size, MPI_INTEGER4, this%root, this%comm, ierr)
  end subroutine

  module subroutine scat_i8_1(this, src, dest)
    class(index_map), intent(in) :: this
    integer(i8), intent(in) :: src(:)
    integer(i8), intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(src) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(dest) >= this%onp_size)
    call MPI_Scatterv(src, this%counts, this%displs, MPI_INTEGER8, &
        dest, this%onp_size, MPI_INTEGER8, this%root, this%comm, ierr)
  end subroutine

  module subroutine scat_r4_1(this, src, dest)
    class(index_map), intent(in) :: this
    real(r4), intent(in) :: src(:)
    real(r4), intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(src) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(dest) >= this%onp_size)
    call MPI_Scatterv(src, this%counts, this%displs, MPI_REAL4, &
        dest, this%onp_size, MPI_REAL4, this%root, this%comm, ierr)
  end subroutine

  module subroutine scat_r8_1(this, src, dest)
    class(index_map), intent(in) :: this
    real(r8), intent(in) :: src(:)
    real(r8), intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(src) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(dest) >= this%onp_size)
    call MPI_Scatterv(src, this%counts, this%displs, MPI_REAL8, &
        dest, this%onp_size, MPI_REAL8, this%root, this%comm, ierr)
  end subroutine

  module subroutine scat_c4_1(this, src, dest)
    class(index_map), intent(in) :: this
    complex(r4), intent(in) :: src(:)
    complex(r4), intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(src) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(dest) >= this%onp_size)
    call MPI_Scatterv(src, this%counts, this%displs, MPI_COMPLEX8, &
        dest, this%onp_size, MPI_COMPLEX8, this%root, this%comm, ierr)
  end subroutine

  module subroutine scat_c8_1(this, src, dest)
    class(index_map), intent(in) :: this
    complex(r8), intent(in) :: src(:)
    complex(r8), intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(src) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(dest) >= this%onp_size)
    call MPI_Scatterv(src, this%counts, this%displs, MPI_COMPLEX16, &
        dest, this%onp_size, MPI_COMPLEX16, this%root, this%comm, ierr)
  end subroutine

  module subroutine scat_dl_1(this, src, dest)
    class(index_map), intent(in) :: this
    logical, intent(in) :: src(:)
    logical, intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(src) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(dest) >= this%onp_size)
    call MPI_Scatterv(src, this%counts, this%displs, MPI_LOGICAL, &
        dest, this%onp_size, MPI_LOGICAL, this%root, this%comm, ierr)
  end subroutine

  !! This auxiliary subroutine performs the scatter operation for elements that
  !! are blocks of values of the specified MPI data type. The block size need
  !! only be specified on the root process. Note that the array dummy arguments
  !! are assumed-size, so if the actual arguments are not contiguous, copy-in/
  !! copy-out of contiguous temporaries will occur.

  subroutine scat_aux2(this, block_size, mpi_type, src, dest)

    type(index_map), intent(in) :: this
    integer, intent(inout) :: block_size
    integer, intent(in) :: mpi_type
    type(*), intent(in) :: src(*)
    type(*), intent(inout) :: dest(*)

    interface ! explicit interface needed to pass assumed-type arguments
      subroutine MPI_Scatterv(sendbuf, sendcounts, displs, sendtype, &
          recvbuf, recvcount, recvtype, root, comm, ierr)
        type(*), intent(in) :: sendbuf(*)
        type(*), intent(inout) :: recvbuf(*)
        integer, intent(in) :: sendcounts(*), displs(*), sendtype, recvcount, recvtype, root, comm
        integer, intent(out) :: ierr
      end subroutine
    end interface

    integer :: ierr
    integer :: block_type

    call MPI_Bcast(block_size, 1, MPI_INTEGER, this%root, this%comm, ierr)
    call MPI_Type_contiguous(block_size, mpi_type, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Scatterv(src, this%counts, this%displs, block_type, &
        dest, this%onp_size, block_type, this%root, this%comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine scat_aux2

!!!! RANK-2 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine scat_i4_2(this, src, dest)
    class(index_map), intent(in) :: this
    integer(i4), intent(in) :: src(:,:)
    integer(i4), intent(inout) :: dest(:,:)
    integer :: block_size
    ASSERT(size(src,2) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(dest,2) >= this%onp_size)
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(src(:,1))
    call scat_aux2(this, block_size, MPI_INTEGER4, src, dest)
  end subroutine

  module subroutine scat_i8_2(this, src, dest)
    class(index_map), intent(in) :: this
    integer(i8), intent(in) :: src(:,:)
    integer(i8), intent(inout) :: dest(:,:)
    integer :: block_size
    ASSERT(size(src,2) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(dest,2) >= this%onp_size)
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(src(:,1))
    call scat_aux2(this, block_size, MPI_INTEGER8, src, dest)
  end subroutine

  module subroutine scat_r4_2(this, src, dest)
    class(index_map), intent(in) :: this
    real(r4), intent(in) :: src(:,:)
    real(r4), intent(inout) :: dest(:,:)
    integer :: block_size
    ASSERT(size(src,2) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(dest,2) >= this%onp_size)
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(src(:,1))
    call scat_aux2(this, block_size, MPI_REAL4, src, dest)
  end subroutine

  module subroutine scat_r8_2(this, src, dest)
    class(index_map), intent(in) :: this
    real(r8), intent(in) :: src(:,:)
    real(r8), intent(inout) :: dest(:,:)
    integer :: block_size
    ASSERT(size(src,2) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(dest,2) >= this%onp_size)
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(src(:,1))
    call scat_aux2(this, block_size, MPI_REAL8, src, dest)
  end subroutine

  module subroutine scat_c4_2(this, src, dest)
    class(index_map), intent(in) :: this
    complex(r4), intent(in) :: src(:,:)
    complex(r4), intent(inout) :: dest(:,:)
    integer :: block_size
    ASSERT(size(src,2) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(dest,2) >= this%onp_size)
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(src(:,1))
    call scat_aux2(this, block_size, MPI_COMPLEX8, src, dest)
  end subroutine

  module subroutine scat_c8_2(this, src, dest)
    class(index_map), intent(in) :: this
    complex(r8), intent(in) :: src(:,:)
    complex(r8), intent(inout) :: dest(:,:)
    integer :: block_size
    ASSERT(size(src,2) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(dest,2) >= this%onp_size)
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(src(:,1))
    call scat_aux2(this, block_size, MPI_COMPLEX16, src, dest)
  end subroutine

  module subroutine scat_dl_2(this, src, dest)
    class(index_map), intent(in) :: this
    logical, intent(in) :: src(:,:)
    logical, intent(inout) :: dest(:,:)
    integer :: block_size
    ASSERT(size(src,2) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(dest,2) >= this%onp_size)
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(src(:,1))
    call scat_aux2(this, block_size, MPI_LOGICAL, src, dest)
  end subroutine

!!!! RANK-3 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine scat_i4_3(this, src, dest)
    class(index_map), intent(in) :: this
    integer(i4), intent(in) :: src(:,:,:)
    integer(i4), intent(inout) :: dest(:,:,:)
    integer :: block_size
    ASSERT(size(src,3) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(dest,3) >= this%onp_size)
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(src(:,:,1))
    call scat_aux2(this, block_size, MPI_INTEGER4, src, dest)
  end subroutine

  module subroutine scat_i8_3(this, src, dest)
    class(index_map), intent(in) :: this
    integer(i8), intent(in) :: src(:,:,:)
    integer(i8), intent(inout) :: dest(:,:,:)
    integer :: block_size
    ASSERT(size(src,3) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(dest,3) >= this%onp_size)
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(src(:,:,1))
    call scat_aux2(this, block_size, MPI_INTEGER8, src, dest)
  end subroutine

  module subroutine scat_r4_3(this, src, dest)
    class(index_map), intent(in) :: this
    real(r4), intent(in) :: src(:,:,:)
    real(r4), intent(inout) :: dest(:,:,:)
    integer :: block_size
    ASSERT(size(src,3) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(dest,3) >= this%onp_size)
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(src(:,:,1))
    call scat_aux2(this, block_size, MPI_REAL4, src, dest)
  end subroutine

  module subroutine scat_r8_3(this, src, dest)
    class(index_map), intent(in) :: this
    real(r8), intent(in) :: src(:,:,:)
    real(r8), intent(inout) :: dest(:,:,:)
    integer :: block_size
    ASSERT(size(src,3) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(dest,3) >= this%onp_size)
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(src(:,:,1))
    call scat_aux2(this, block_size, MPI_REAL8, src, dest)
  end subroutine

  module subroutine scat_c4_3(this, src, dest)
    class(index_map), intent(in) :: this
    complex(r4), intent(in) :: src(:,:,:)
    complex(r4), intent(inout) :: dest(:,:,:)
    integer :: block_size
    ASSERT(size(src,3) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(dest,3) >= this%onp_size)
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(src(:,:,1))
    call scat_aux2(this, block_size, MPI_COMPLEX8, src, dest)
  end subroutine

  module subroutine scat_c8_3(this, src, dest)
    class(index_map), intent(in) :: this
    complex(r8), intent(in) :: src(:,:,:)
    complex(r8), intent(inout) :: dest(:,:,:)
    integer :: block_size
    ASSERT(size(src,3) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(dest,3) >= this%onp_size)
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(src(:,:,1))
    call scat_aux2(this, block_size, MPI_COMPLEX16, src, dest)
  end subroutine

  module subroutine scat_dl_3(this, src, dest)
    class(index_map), intent(in) :: this
    logical, intent(in) :: src(:,:,:)
    logical, intent(inout) :: dest(:,:,:)
    integer :: block_size
    ASSERT(size(src,3) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(dest,3) >= this%onp_size)
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(src(:,:,1))
    call scat_aux2(this, block_size, MPI_LOGICAL, src, dest)
  end subroutine

end submodule scatter_impl
