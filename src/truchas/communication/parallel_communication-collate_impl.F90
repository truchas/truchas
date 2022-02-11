!! Implementation of PARALLEL_COMMUNICATION COLLATE Procedures
!!
!! Copyright 2022 Neil N. Carlson <neil.n.carlson@gmail.com>
!! Use subject to the MIT license: https://opensource.org/licenses/MIT
!!

#include "f90_assert.fpp"

submodule(parallel_communication) collate_impl
implicit none
contains

!!!! COLLATE SCALAR DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine coll_i1_0(scalarv_out, scalar_in)
    integer(int8), intent(inout) :: scalarv_out(:)
    integer(int8), intent(in) :: scalar_in
    integer :: ierr
    ASSERT(size(scalarv_out) >= merge(npe,0,is_iop))
    call MPI_Gather(scalar_in, 1, MPI_INTEGER1, scalarv_out, 1, MPI_INTEGER1, root, comm, ierr)
  end subroutine

  module subroutine coll_i4_0(scalarv_out, scalar_in)
    integer(int32), intent(inout) :: scalarv_out(:)
    integer(int32), intent(in) :: scalar_in
    integer :: ierr
    ASSERT(size(scalarv_out) >= merge(npe,0,is_iop))
    call MPI_Gather(scalar_in, 1, MPI_INTEGER4, scalarv_out, 1, MPI_INTEGER4, root, comm, ierr)
  end subroutine

  module subroutine coll_i8_0(scalarv_out, scalar_in)
    integer(int64), intent(inout) :: scalarv_out(:)
    integer(int64), intent(in) :: scalar_in
    integer :: ierr
    ASSERT(size(scalarv_out) >= merge(npe,0,is_iop))
    call MPI_Gather(scalar_in, 1, MPI_INTEGER8, scalarv_out, 1, MPI_INTEGER8, root, comm, ierr)
  end subroutine

  module subroutine coll_r4_0(scalarv_out, scalar_in)
    real(real32), intent(inout) :: scalarv_out(:)
    real(real32), intent(in) :: scalar_in
    integer :: ierr
    ASSERT(size(scalarv_out) >= merge(npe,0,is_iop))
    call MPI_Gather(scalar_in, 1, MPI_REAL4, scalarv_out, 1, MPI_REAL4, root, comm, ierr)
  end subroutine

  module subroutine coll_r8_0(scalarv_out, scalar_in)
    real(real64), intent(inout) :: scalarv_out(:)
    real(real64), intent(in) :: scalar_in
    integer :: ierr
    ASSERT(size(scalarv_out) >= merge(npe,0,is_iop))
    call MPI_Gather(scalar_in, 1, MPI_REAL8, scalarv_out, 1, MPI_REAL8, root, comm, ierr)
  end subroutine

  module subroutine coll_log_0(scalarv_out, scalar_in)
    logical, intent(inout) :: scalarv_out(:)
    logical, intent(in) :: scalar_in
    integer :: ierr
    ASSERT(size(scalarv_out) >= merge(npe,0,is_iop))
    call MPI_Gather(scalar_in, 1, MPI_LOGICAL, scalarv_out, 1, MPI_LOGICAL, root, comm, ierr)
  end subroutine

  module subroutine coll_char_0(scalarv_out, scalar_in)
    character(*), intent(inout) :: scalarv_out(:)
    character(*), intent(in) :: scalar_in
    integer :: string_len, ierr
    ASSERT(size(scalarv_out) >= merge(npe,0,is_iop))
    string_len = len(scalar_in) ! must be same for all character args on all processes!
    call MPI_Gather(scalar_in, string_len, MPI_CHARACTER, scalarv_out, string_len, MPI_CHARACTER, root, comm, ierr)
  end subroutine

!!!! COLLATE RANK-1 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine coll_i1_1(vector_out, vector_in)

    integer(int8), intent(inout) :: vector_out(:)
    integer(int8), intent(in) :: vector_in(:)

    integer :: inlen, outlen, ierr, j
    integer, allocatable :: counts(:), displs(:)

    inlen = size(vector_in, dim=1)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    outlen = sum(counts)
    call MPI_Bcast(outlen, 1, MPI_INTEGER, root, comm, ierr)
    if (outlen == 0) return ! nothing to do

    ASSERT(size(vector_out,dim=1) >= merge(outlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Gatherv(vector_in, inlen, MPI_INTEGER1, &
        vector_out, counts, displs, MPI_INTEGER1, root, comm, ierr)

  end subroutine

  module subroutine coll_i4_1(vector_out, vector_in)

    integer(int32), intent(inout) :: vector_out(:)
    integer(int32), intent(in) :: vector_in(:)

    integer :: inlen, outlen, ierr, j
    integer, allocatable :: counts(:), displs(:)

    inlen = size(vector_in, dim=1)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    outlen = sum(counts)
    call MPI_Bcast(outlen, 1, MPI_INTEGER, root, comm, ierr)
    if (outlen == 0) return ! nothing to do

    ASSERT(size(vector_out,dim=1) >= merge(outlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Gatherv(vector_in, inlen, MPI_INTEGER4, &
        vector_out, counts, displs, MPI_INTEGER4, root, comm, ierr)

  end subroutine

  module subroutine coll_i8_1(vector_out, vector_in)

    integer(int64), intent(inout) :: vector_out(:)
    integer(int64), intent(in) :: vector_in(:)

    integer :: inlen, outlen, ierr, j
    integer, allocatable :: counts(:), displs(:)

    inlen = size(vector_in, dim=1)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    outlen = sum(counts)
    call MPI_Bcast(outlen, 1, MPI_INTEGER, root, comm, ierr)
    if (outlen == 0) return ! nothing to do

    ASSERT(size(vector_out,dim=1) >= merge(outlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Gatherv(vector_in, inlen, MPI_INTEGER8, &
        vector_out, counts, displs, MPI_INTEGER8, root, comm, ierr)

  end subroutine

  module subroutine coll_r4_1(vector_out, vector_in)

    real(real32), intent(inout) :: vector_out(:)
    real(real32), intent(in) :: vector_in(:)

    integer :: inlen, outlen, ierr, j
    integer, allocatable :: counts(:), displs(:)

    inlen = size(vector_in, dim=1)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    outlen = sum(counts)
    call MPI_Bcast(outlen, 1, MPI_INTEGER, root, comm, ierr)
    if (outlen == 0) return ! nothing to do

    ASSERT(size(vector_out,dim=1) >= merge(outlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Gatherv(vector_in, inlen, MPI_REAL4, &
        vector_out, counts, displs, MPI_REAL4, root, comm, ierr)

  end subroutine

  module subroutine coll_r8_1(vector_out, vector_in)

    real(real64), intent(inout) :: vector_out(:)
    real(real64), intent(in) :: vector_in(:)

    integer :: inlen, outlen, ierr, j
    integer, allocatable :: counts(:), displs(:)

    inlen = size(vector_in, dim=1)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    outlen = sum(counts)
    call MPI_Bcast(outlen, 1, MPI_INTEGER, root, comm, ierr)
    if (outlen == 0) return ! nothing to do

    ASSERT(size(vector_out,dim=1) >= merge(outlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Gatherv(vector_in, inlen, MPI_REAL8, &
        vector_out, counts, displs, MPI_REAL8, root, comm, ierr)

  end subroutine

  module subroutine coll_log_1(vector_out, vector_in)

    logical, intent(inout) :: vector_out(:)
    logical, intent(in) :: vector_in(:)

    integer :: inlen, outlen, ierr, j
    integer, allocatable :: counts(:), displs(:)

    inlen = size(vector_in, dim=1)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    outlen = sum(counts)
    call MPI_Bcast(outlen, 1, MPI_INTEGER, root, comm, ierr)
    if (outlen == 0) return ! nothing to do

    ASSERT(size(vector_out,dim=1) >= merge(outlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Gatherv(vector_in, inlen, MPI_LOGICAL, &
        vector_out, counts, displs, MPI_LOGICAL, root, comm, ierr)

  end subroutine

  module subroutine coll_char_1(vector_out, vector_in)

    character(*), intent(inout) :: vector_out(:)
    character(*), intent(in) :: vector_in(:)

    integer :: inlen, outlen, ierr, j
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    inlen = size(vector_in, dim=1)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    outlen = sum(counts)
    call MPI_Bcast(outlen, 1, MPI_INTEGER, root, comm, ierr)
    if (outlen == 0) return ! nothing to do

    ASSERT(size(vector_out,dim=1) >= merge(outlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Type_contiguous(len(vector_in), MPI_CHARACTER, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Gatherv(vector_in, inlen, block_type, &
        vector_out, counts, displs, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine

!!!! COLLATE RANK-2 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine coll_i1_2(vector_out, vector_in)

    integer(int8), intent(inout) :: vector_out(:,:)
    integer(int8), intent(in) :: vector_in(:,:)

    integer :: inlen, outlen, ierr, j, block_size
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    inlen = size(vector_in, dim=2)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    outlen = sum(counts)
    call MPI_Bcast(outlen, 1, MPI_INTEGER, root, comm, ierr)
    if (outlen == 0) return ! nothing to do

    ASSERT(size(vector_out,dim=2) >= merge(outlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    if (is_iop) block_size = size(vector_out(:,1))
    call MPI_Bcast(block_size, 1, MPI_INTEGER, root, comm, ierr)
    call MPI_Type_contiguous(block_size, MPI_INTEGER1, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Gatherv(vector_in, inlen, block_type, &
        vector_out, counts, displs, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine

  module subroutine coll_i4_2(vector_out, vector_in)

    integer(int32), intent(inout) :: vector_out(:,:)
    integer(int32), intent(in) :: vector_in(:,:)

    integer :: inlen, outlen, ierr, j, block_size
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    inlen = size(vector_in, dim=2)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    outlen = sum(counts)
    call MPI_Bcast(outlen, 1, MPI_INTEGER, root, comm, ierr)
    if (outlen == 0) return ! nothing to do

    ASSERT(size(vector_out,dim=2) >= merge(outlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    if (is_iop) block_size = size(vector_out(:,1))
    call MPI_Bcast(block_size, 1, MPI_INTEGER, root, comm, ierr)
    call MPI_Type_contiguous(block_size, MPI_INTEGER4, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Gatherv(vector_in, inlen, block_type, &
        vector_out, counts, displs, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine

  module subroutine coll_i8_2(vector_out, vector_in)

    integer(int64), intent(inout) :: vector_out(:,:)
    integer(int64), intent(in) :: vector_in(:,:)

    integer :: inlen, outlen, ierr, j, block_size
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    inlen = size(vector_in, dim=2)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    outlen = sum(counts)
    call MPI_Bcast(outlen, 1, MPI_INTEGER, root, comm, ierr)
    if (outlen == 0) return ! nothing to do

    ASSERT(size(vector_out,dim=2) >= merge(outlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    if (is_iop) block_size = size(vector_out(:,1))
    call MPI_Bcast(block_size, 1, MPI_INTEGER, root, comm, ierr)
    call MPI_Type_contiguous(block_size, MPI_INTEGER8, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Gatherv(vector_in, inlen, block_type, &
        vector_out, counts, displs, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine

  module subroutine coll_r4_2(vector_out, vector_in)

    real(real32), intent(inout) :: vector_out(:,:)
    real(real32), intent(in) :: vector_in(:,:)

    integer :: inlen, outlen, ierr, j, block_size
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    inlen = size(vector_in, dim=2)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    outlen = sum(counts)
    call MPI_Bcast(outlen, 1, MPI_INTEGER, root, comm, ierr)
    if (outlen == 0) return ! nothing to do

    ASSERT(size(vector_out,dim=2) >= merge(outlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    if (is_iop) block_size = size(vector_out(:,1))
    call MPI_Bcast(block_size, 1, MPI_INTEGER, root, comm, ierr)
    call MPI_Type_contiguous(block_size, MPI_REAL4, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Gatherv(vector_in, inlen, block_type, &
        vector_out, counts, displs, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine

  module subroutine coll_r8_2(vector_out, vector_in)

    real(real64), intent(inout) :: vector_out(:,:)
    real(real64), intent(in) :: vector_in(:,:)

    integer :: inlen, outlen, ierr, j, block_size
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    inlen = size(vector_in, dim=2)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    outlen = sum(counts)
    call MPI_Bcast(outlen, 1, MPI_INTEGER, root, comm, ierr)
    if (outlen == 0) return ! nothing to do

    ASSERT(size(vector_out,dim=2) >= merge(outlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    if (is_iop) block_size = size(vector_out(:,1))
    call MPI_Bcast(block_size, 1, MPI_INTEGER, root, comm, ierr)
    call MPI_Type_contiguous(block_size, MPI_REAL8, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Gatherv(vector_in, inlen, block_type, &
        vector_out, counts, displs, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine

  module subroutine coll_log_2(vector_out, vector_in)

    logical, intent(inout) :: vector_out(:,:)
    logical, intent(in) :: vector_in(:,:)

    integer :: inlen, outlen, ierr, j, block_size
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    inlen = size(vector_in, dim=2)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    outlen = sum(counts)
    call MPI_Bcast(outlen, 1, MPI_INTEGER, root, comm, ierr)
    if (outlen == 0) return ! nothing to do

    ASSERT(size(vector_out,dim=2) >= merge(outlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    if (is_iop) block_size = size(vector_out(:,1))
    call MPI_Bcast(block_size, 1, MPI_INTEGER, root, comm, ierr)
    call MPI_Type_contiguous(block_size, MPI_LOGICAL, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Gatherv(vector_in, inlen, block_type, &
        vector_out, counts, displs, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine

  module subroutine coll_char_2(vector_out, vector_in)

    character(*), intent(inout) :: vector_out(:,:)
    character(*), intent(in) :: vector_in(:,:)

    integer :: inlen, outlen, ierr, j, block_size
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    inlen = size(vector_in, dim=2)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    outlen = sum(counts)
    call MPI_Bcast(outlen, 1, MPI_INTEGER, root, comm, ierr)
    if (outlen == 0) return ! nothing to do

    ASSERT(size(vector_out,dim=2) >= merge(outlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    if (is_iop) block_size = len(vector_out)*size(vector_out(:,1))
    call MPI_Bcast(block_size, 1, MPI_INTEGER, root, comm, ierr)
    call MPI_Type_contiguous(block_size, MPI_CHARACTER, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Gatherv(vector_in, inlen, block_type, &
        vector_out, counts, displs, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine

end submodule collate_impl
