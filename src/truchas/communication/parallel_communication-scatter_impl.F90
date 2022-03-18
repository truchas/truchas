!! Implementation of PARALLEL_COMMUNICATION SCATTER Procedures
!!
!! Copyright 2022 Neil N. Carlson <neil.n.carlson@gmail.com>
!! Use subject to the MIT license: https://opensource.org/licenses/MIT
!!

#include "f90_assert.fpp"

submodule(parallel_communication) scatter_impl
implicit none
contains

!!!! SCATTER SCALAR DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine scat_i4_0(src, dest)
    integer(i4), intent(in)  :: src(:)
    integer(i4), intent(out) :: dest
    integer :: ierr
    ASSERT(size(src) >= merge(npe,0,is_iop))
    call MPI_Scatter(src, 1, MPI_INTEGER4, dest, 1, MPI_INTEGER4, root, comm, ierr)
  end subroutine

  module subroutine scat_i8_0(src, dest)
    integer(i8), intent(in)  :: src(:)
    integer(i8), intent(out) :: dest
    integer :: ierr
    ASSERT(size(src) >= merge(npe,0,is_iop))
    call MPI_Scatter(src, 1, MPI_INTEGER8, dest, 1, MPI_INTEGER8, root, comm, ierr)
  end subroutine

  module subroutine scat_r4_0(src, dest)
    real(r4), intent(in)  :: src(:)
    real(r4), intent(out) :: dest
    integer :: ierr
    ASSERT(size(src) >= merge(npe,0,is_iop))
    call MPI_Scatter(src, 1, MPI_REAL4, dest, 1, MPI_REAL4, root, comm, ierr)
  end subroutine

  module subroutine scat_r8_0(src, dest)
    real(r8), intent(in)  :: src(:)
    real(r8), intent(out) :: dest
    integer :: ierr
    ASSERT(size(src) >= merge(npe,0,is_iop))
    call MPI_Scatter(src, 1, MPI_REAL8, dest, 1, MPI_REAL8, root, comm, ierr)
  end subroutine

  module subroutine scat_dl_0(src, dest)
    logical, intent(in)  :: src(:)
    logical, intent(out) :: dest
    integer :: ierr
    ASSERT(size(src) >= merge(npe,0,is_iop))
    call MPI_Scatter(src, 1, MPI_LOGICAL, dest, 1, MPI_LOGICAL, root, comm, ierr)
  end subroutine

  !! This auxiliary subroutine generates the COUNTS and DISPLS arrays for the
  !! send buffer on the IO process, given the size of the receive buffer OUTLEN
  !! on each of the processes. INLEN is the required size of the send buffer.

  subroutine make_counts_displs(outlen, counts, displs, inlen)

    integer, intent(in) :: outlen
    integer, allocatable :: counts(:), displs(:)
    integer, intent(out) :: inlen

    integer :: ierr, j

    allocate(counts(merge(npe,0,is_iop)))
    call MPI_Gather(outlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    inlen = sum(counts)
    call MPI_Bcast(inlen, 1, MPI_INTEGER, root, comm, ierr)
    if (inlen == 0) return ! nothing to do

    allocate(displs(merge(npe,0,is_iop)))
    if (is_iop) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

  end subroutine

  !! This auxiliary subroutine performs the scatter operation for unit-sized
  !! elements of the specified MPI data type. Note that the array dummy
  !! arguments are assumed-size, so if the actual arguments are not contiguous,
  !! copy-in/copy-out of contiguous temporaries will occur.

  subroutine scat_aux1(srclen, outlen, mpi_type, src, dest)

    integer, intent(in) :: srclen, outlen
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

    integer :: inlen, ierr
    integer, allocatable :: counts(:), displs(:)

    call make_counts_displs(outlen, counts, displs, inlen)
    if (inlen == 0) return ! nothing to do
    ASSERT(srclen >= merge(inlen,0,is_iop))

    call MPI_Scatterv(src, counts, displs, mpi_type, dest, outlen, mpi_type, root, comm, ierr)

  end subroutine scat_aux1

  !! This auxiliary subroutine performs the scatter operation for elements that
  !! are blocks of values of the specified MPI data type. The block size need
  !! only be specified on the root process. Note that the array ummy arguments
  !! are assumed-size, so if the actual arguments are not contiguous, copy-in/
  !! copy-out of contiguous temporaries will occur.

  subroutine scat_aux2(srclen, outlen, block_size, mpi_type, src, dest)

    integer, intent(in) :: srclen, outlen, block_size
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

    integer :: inlen, ierr
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    call make_counts_displs(outlen, counts, displs, inlen)
    if (inlen == 0) return ! nothing to do
    ASSERT(srclen >= merge(inlen,0,is_iop))

    call MPI_Bcast(block_size, 1, MPI_INTEGER, root, comm, ierr)
    call MPI_Type_contiguous(block_size, mpi_type, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Scatterv(src, counts, displs, block_type, dest, outlen, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine scat_aux2

!!!! SCATTER RANK-1 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine scat_i1_1(src, dest)
    integer(i1), intent(in)  :: src(:)
    integer(i1), intent(out) :: dest(:)
    call scat_aux1(size(src), size(dest), MPI_INTEGER1, src, dest)
  end subroutine scat_i1_1

  module subroutine scat_i4_1(src, dest)
    integer(i4), intent(in)  :: src(:)
    integer(i4), intent(out) :: dest(:)
    call scat_aux1(size(src), size(dest), MPI_INTEGER4, src, dest)
  end subroutine scat_i4_1

  module subroutine scat_i8_1(src, dest)
    integer(i8), intent(in)  :: src(:)
    integer(i8), intent(out) :: dest(:)
    call scat_aux1(size(src), size(dest), MPI_INTEGER8, src, dest)
  end subroutine scat_i8_1

  module subroutine scat_r4_1(src, dest)
    real(r4), intent(in)  :: src(:)
    real(r4), intent(out) :: dest(:)
    call scat_aux1(size(src), size(dest), MPI_REAL4, src, dest)
  end subroutine scat_r4_1

  module subroutine scat_r8_1(src, dest)
    real(r8), intent(in)  :: src(:)
    real(r8), intent(out) :: dest(:)
    call scat_aux1(size(src), size(dest), MPI_REAL8, src, dest)
  end subroutine scat_r8_1

  module subroutine scat_dl_1(src, dest)
    logical, intent(in)  :: src(:)
    logical, intent(out) :: dest(:)
    call scat_aux1(size(src), size(dest), MPI_LOGICAL, src, dest)
  end subroutine scat_dl_1

!!!! SCATTER RANK-2 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine scat_i1_2(src, dest)
    integer(i1), intent(in)  :: src(:,:)
    integer(i1), intent(out) :: dest(:,:)
    integer :: destlen, srclen, block_size
    destlen = size(dest, dim=2)
    srclen = size(src, dim=2)
    if (is_iop) block_size = size(src(:,1))
    call scat_aux2(srclen, destlen, block_size, MPI_INTEGER1, src, dest)
  end subroutine scat_i1_2

  module subroutine scat_i4_2(src, dest)
    integer(i4), intent(in)  :: src(:,:)
    integer(i4), intent(out) :: dest(:,:)
    integer :: destlen, srclen, block_size
    destlen = size(dest, dim=2)
    srclen = size(src, dim=2)
    if (is_iop) block_size = size(src(:,1))
    call scat_aux2(srclen, destlen, block_size, MPI_INTEGER4, src, dest)
  end subroutine scat_i4_2

  module subroutine scat_i8_2(src, dest)
    integer(i8), intent(in)  :: src(:,:)
    integer(i8), intent(out) :: dest(:,:)
    integer :: srclen, destlen, block_size
    destlen = size(dest, dim=2)
    srclen = size(src, dim=2)
    if (is_iop) block_size = size(src(:,1))
    call scat_aux2(srclen, destlen, block_size, MPI_INTEGER8, src, dest)
  end subroutine scat_i8_2

  module subroutine scat_r4_2(src, dest)
    real(r4), intent(in)  :: src(:,:)
    real(r4), intent(out) :: dest(:,:)
    integer :: srclen, destlen, block_size
    destlen = size(dest, dim=2)
    srclen = size(src, dim=2)
    if (is_iop) block_size = size(src(:,1))
    call scat_aux2(srclen, destlen, block_size, MPI_REAL4, src, dest)
  end subroutine scat_r4_2

  module subroutine scat_r8_2(src, dest)
    real(r8), intent(in)  :: src(:,:)
    real(r8), intent(out) :: dest(:,:)
    integer :: srclen, destlen, block_size
    destlen = size(dest, dim=2)
    srclen = size(src, dim=2)
    if (is_iop) block_size = size(src(:,1))
    call scat_aux2(srclen, destlen, block_size, MPI_REAL8, src, dest)
  end subroutine scat_r8_2

  module subroutine scat_dl_2(src, dest)
    logical, intent(in)  :: src(:,:)
    logical, intent(out) :: dest(:,:)
    integer :: srclen, destlen, block_size
    destlen = size(dest, dim=2)
    srclen = size(src, dim=2)
    if (is_iop) block_size = size(src(:,1))
    call scat_aux2(srclen, destlen, block_size, MPI_LOGICAL, src, dest)
  end subroutine scat_dl_2

!!!! SCATTER RANK-3 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine scat_i1_3(src, dest)
    integer(i1), intent(in)  :: src(:,:,:)
    integer(i1), intent(out) :: dest(:,:,:)
    integer :: srclen, destlen, block_size
    destlen = size(dest, dim=3)
    srclen = size(src, dim=3)
    if (is_iop) block_size = size(src(:,:,1))
    call scat_aux2(srclen, destlen, block_size, MPI_INTEGER1, src, dest)
  end subroutine scat_i1_3

  module subroutine scat_i4_3(src, dest)
    integer(i4), intent(in)  :: src(:,:,:)
    integer(i4), intent(out) :: dest(:,:,:)
    integer :: srclen, destlen, block_size
    destlen = size(dest, dim=3)
    srclen = size(src, dim=3)
    if (is_iop) block_size = size(src(:,:,1))
    call scat_aux2(srclen, destlen, block_size, MPI_INTEGER4, src, dest)
  end subroutine scat_i4_3

  module subroutine scat_i8_3(src, dest)
    integer(i8), intent(in)  :: src(:,:,:)
    integer(i8), intent(out) :: dest(:,:,:)
    integer :: srclen, destlen, block_size
    destlen = size(dest, dim=3)
    srclen = size(src, dim=3)
    if (is_iop) block_size = size(src(:,:,1))
    call scat_aux2(srclen, destlen, block_size, MPI_INTEGER8, src, dest)
  end subroutine scat_i8_3

  module subroutine scat_r4_3(src, dest)
    real(r4), intent(in)  :: src(:,:,:)
    real(r4), intent(out) :: dest(:,:,:)
    integer :: srclen, destlen, block_size
    destlen = size(dest, dim=3)
    srclen = size(src, dim=3)
    if (is_iop) block_size = size(src(:,:,1))
    call scat_aux2(srclen, destlen, block_size, MPI_REAL4, src, dest)
  end subroutine scat_r4_3

  module subroutine scat_r8_3(src, dest)
    real(r8), intent(in)  :: src(:,:,:)
    real(r8), intent(out) :: dest(:,:,:)
    integer :: srclen, destlen, block_size
    destlen = size(dest, dim=3)
    srclen = size(src, dim=3)
    if (is_iop) block_size = size(src(:,:,1))
    call scat_aux2(srclen, destlen, block_size, MPI_REAL8, src, dest)
  end subroutine scat_r8_3

  module subroutine scat_dl_3(src, dest)
    logical, intent(in)  :: src(:,:,:)
    logical, intent(out) :: dest(:,:,:)
    integer :: srclen, destlen, block_size
    destlen = size(dest, dim=3)
    srclen = size(src, dim=3)
    if (is_iop) block_size = size(src(:,:,1))
    call scat_aux2(srclen, destlen, block_size, MPI_LOGICAL, src, dest)
  end subroutine scat_dl_3

end submodule scatter_impl
