!! Implementation of PARALLEL_COMMUNICATION GATHER Procedures
!!
!! Copyright 2022 Neil N. Carlson <neil.n.carlson@gmail.com>
!! Use subject to the MIT license: https://opensource.org/licenses/MIT
!!

#include "f90_assert.fpp"

submodule(parallel_communication) gather_impl
implicit none

contains

!!!! GATHER SCALAR DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine gath_i1_0(src, dest)
    integer(i1), intent(in) :: src
    integer(i1), intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(dest) >= merge(npe,0,is_iop))
    call MPI_Gather(src, 1, MPI_INTEGER1, dest, 1, MPI_INTEGER1, root, comm, ierr)
  end subroutine

  module subroutine gath_i4_0(src, dest)
    integer(i4), intent(in) :: src
    integer(i4), intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(dest) >= merge(npe,0,is_iop))
    call MPI_Gather(src, 1, MPI_INTEGER4, dest, 1, MPI_INTEGER4, root, comm, ierr)
  end subroutine

  module subroutine gath_i8_0(src, dest)
    integer(i8), intent(in) :: src
    integer(i8), intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(dest) >= merge(npe,0,is_iop))
    call MPI_Gather(src, 1, MPI_INTEGER8, dest, 1, MPI_INTEGER8, root, comm, ierr)
  end subroutine

  module subroutine gath_r4_0(src, dest)
    real(r4), intent(in) :: src
    real(r4), intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(dest) >= merge(npe,0,is_iop))
    call MPI_Gather(src, 1, MPI_REAL4, dest, 1, MPI_REAL4, root, comm, ierr)
  end subroutine

  module subroutine gath_r8_0(src, dest)
    real(r8), intent(in) :: src
    real(r8), intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(dest) >= merge(npe,0,is_iop))
    call MPI_Gather(src, 1, MPI_REAL8, dest, 1, MPI_REAL8, root, comm, ierr)
  end subroutine

  module subroutine gath_dl_0(src, dest)
    logical, intent(in) :: src
    logical, intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(dest) >= merge(npe,0,is_iop))
    call MPI_Gather(src, 1, MPI_LOGICAL, dest, 1, MPI_LOGICAL, root, comm, ierr)
  end subroutine

  module subroutine gath_char_0(src, dest)
    character(*), intent(in) :: src
    character(*), intent(inout) :: dest(:)
    integer :: string_len, ierr
    ASSERT(size(dest) >= merge(npe,0,is_iop))
    string_len = len(src) ! must be same for all character args on all processes!
    call MPI_Gather(src, string_len, MPI_CHARACTER, dest, string_len, MPI_CHARACTER, root, comm, ierr)
  end subroutine

  !! This auxiliary subroutine generates the COUNTS and DISPLS arrays for the
  !! receive buffer on the IO process, given the size of the send buffer INLEN
  !! on each of the processes. OUTLEN is the required size of the receive buffer.

  subroutine make_counts_displs(inlen, counts, displs, outlen)

    integer, intent(in) :: inlen
    integer, allocatable :: counts(:), displs(:)
    integer, intent(out) :: outlen

    integer :: ierr, j

    allocate(counts(merge(npe,0,is_iop)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    outlen = sum(counts)
    call MPI_Bcast(outlen, 1, MPI_INTEGER, root, comm, ierr)
    if (outlen == 0) return ! nothing to do

    allocate(displs(merge(npe,0,is_iop)))
    if (is_iop) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

  end subroutine

  !! This auxiliary subroutine performs the gather operation for unit-sized
  !! elements of the specified MPI data type. Note that the array dummy
  !! arguments are assumed-size, so if the actual arguments are not contiguous,
  !! copy-in/copy-out of contiguous temporaries will occur.

  subroutine gath_aux1(inlen, destlen, mpi_type, src, dest)

    integer, intent(in) :: inlen, destlen
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

    integer :: outlen, ierr
    integer, allocatable :: counts(:), displs(:)

    call make_counts_displs(inlen, counts, displs, outlen)
    if (outlen == 0) return ! nothing to do
    ASSERT(destlen >= merge(outlen,0,is_iop))

    call MPI_Gatherv(src, inlen, mpi_type, dest, counts, displs, mpi_type, root, comm, ierr)

  end subroutine gath_aux1

  !! This auxiliary subroutine performs the gather operation for elements that
  !! are blocks of values of the specified MPI data type. The block size need
  !! only be specified on the root process. Note that the array ummy arguments
  !! are assumed-size, so if the actual arguments are not contiguous, copy-in/
  !! copy-out of contiguous temporaries will occur.

  subroutine gath_aux2(inlen, destlen, block_size, mpi_type, src, dest)

    integer, intent(in) :: inlen, destlen, block_size
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

    integer :: outlen, ierr
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    call make_counts_displs(inlen, counts, displs, outlen)
    if (outlen == 0) return ! nothing to do
    ASSERT(destlen >= merge(outlen,0,is_iop))

    call MPI_Bcast(block_size, 1, MPI_INTEGER, root, comm, ierr)
    call MPI_Type_contiguous(block_size, mpi_type, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Gatherv(src, inlen, block_type, dest, counts, displs, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine gath_aux2

!!!! GATHER RANK-1 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine gath_i1_1(src, dest)
    integer(i1), intent(in) :: src(:)
    integer(i1), intent(inout) :: dest(:)
    call gath_aux1(size(src), size(dest), MPI_INTEGER1, src, dest)
  end subroutine

  module subroutine gath_i4_1(src, dest)
    integer(i4), intent(in) :: src(:)
    integer(i4), intent(inout) :: dest(:)
    call gath_aux1(size(src), size(dest), MPI_INTEGER4, src, dest)
  end subroutine

  module subroutine gath_i8_1(src, dest)
    integer(i8), intent(in) :: src(:)
    integer(i8), intent(inout) :: dest(:)
    call gath_aux1(size(src), size(dest), MPI_INTEGER8, src, dest)
  end subroutine

  module subroutine gath_r4_1(src, dest)
    real(r4), intent(in) :: src(:)
    real(r4), intent(inout) :: dest(:)
    call gath_aux1(size(src), size(dest), MPI_REAL4, src, dest)
  end subroutine

  module subroutine gath_r8_1(src, dest)
    real(r8), intent(in) :: src(:)
    real(r8), intent(inout) :: dest(:)
    call gath_aux1(size(src), size(dest), MPI_REAL8, src, dest)
  end subroutine

  module subroutine gath_dl_1(src, dest)
    logical, intent(in) :: src(:)
    logical, intent(inout) :: dest(:)
    call gath_aux1(size(src), size(dest), MPI_LOGICAL, src, dest)
  end subroutine

  module subroutine gath_char_1(src, dest)
    character(*), intent(in) :: src(:)
    character(*), intent(inout) :: dest(:)
    integer :: block_size
    if (is_iop) block_size = len(src)
    call gath_aux2(size(src), size(dest), block_size, MPI_CHARACTER, src, dest)
  end subroutine

!!!! GATHER RANK-2 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine gath_i1_2(src, dest)
    integer(i1), intent(in) :: src(:,:)
    integer(i1), intent(inout) :: dest(:,:)
    integer :: srclen, destlen, block_size
    srclen = size(src, dim=2)
    destlen = size(dest, dim=2)
    if (is_iop) block_size = size(dest(:,1))
    call gath_aux2(srclen, destlen, block_size, MPI_INTEGER1, src, dest)
  end subroutine

  module subroutine gath_i4_2(src, dest)
    integer(i4), intent(in) :: src(:,:)
    integer(i4), intent(inout) :: dest(:,:)
    integer :: srclen, destlen, block_size
    srclen = size(src, dim=2)
    destlen = size(dest, dim=2)
    if (is_iop) block_size = size(dest(:,1))
    call gath_aux2(srclen, destlen, block_size, MPI_INTEGER4, src, dest)
  end subroutine

  module subroutine gath_i8_2(src, dest)
    integer(i8), intent(in) :: src(:,:)
    integer(i8), intent(inout) :: dest(:,:)
    integer :: srclen, destlen, block_size
    srclen = size(src, dim=2)
    destlen = size(dest, dim=2)
    if (is_iop) block_size = size(dest(:,1))
    call gath_aux2(srclen, destlen, block_size, MPI_INTEGER8, src, dest)
  end subroutine

  module subroutine gath_r4_2(src, dest)
    real(r4), intent(in) :: src(:,:)
    real(r4), intent(inout) :: dest(:,:)
    integer :: srclen, destlen, block_size
    srclen = size(src, dim=2)
    destlen = size(dest, dim=2)
    if (is_iop) block_size = size(dest(:,1))
    call gath_aux2(srclen, destlen, block_size, MPI_REAL4, src, dest)
  end subroutine

  module subroutine gath_r8_2(src, dest)
    real(r8), intent(in) :: src(:,:)
    real(r8), intent(inout) :: dest(:,:)
    integer :: srclen, destlen, block_size
    srclen = size(src, dim=2)
    destlen = size(dest, dim=2)
    if (is_iop) block_size = size(dest(:,1))
    call gath_aux2(srclen, destlen, block_size, MPI_REAL8, src, dest)
  end subroutine

  module subroutine gath_dl_2(src, dest)
    logical, intent(in) :: src(:,:)
    logical, intent(inout) :: dest(:,:)
    integer :: srclen, destlen, block_size
    srclen = size(src, dim=2)
    destlen = size(dest, dim=2)
    if (is_iop) block_size = size(dest(:,1))
    call gath_aux2(srclen, destlen, block_size, MPI_LOGICAL, src, dest)
  end subroutine

  module subroutine gath_char_2(src, dest)
    character(*), intent(in) :: src(:,:)
    character(*), intent(inout) :: dest(:,:)
    integer :: srclen, destlen, block_size
    srclen = size(src, dim=2)
    destlen = size(dest, dim=2)
    if (is_iop) block_size = len(dest)*size(dest(:,1))
    call gath_aux2(srclen, destlen, block_size, MPI_CHARACTER, src, dest)
  end subroutine

end submodule gather_impl
