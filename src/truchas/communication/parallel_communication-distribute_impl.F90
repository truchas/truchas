!! Implementation of PARALLEL_COMMUNICATION DISTRIBUTE Procedures
!!
!! Copyright 2022 Neil N. Carlson <neil.n.carlson@gmail.com>
!! Use subject to the MIT license: https://opensource.org/licenses/MIT
!!

#include "f90_assert.fpp"

submodule(parallel_communication) distribute_impl
implicit none
contains

!!!! DISTRIBUTE SCALAR DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine dist_i4_0(scalar_out, scalarv_in)
    integer(int32), intent(out) :: scalar_out
    integer(int32), intent(in)  :: scalarv_in(:)
    integer :: ierr
    ASSERT(size(scalarv_in) >= merge(npe,0,is_iop))
    call MPI_Scatter(scalarv_in, 1, MPI_INTEGER4, scalar_out, 1, MPI_INTEGER4, root, comm, ierr)
  end subroutine

  module subroutine dist_i8_0(scalar_out, scalarv_in)
    integer(int64), intent(out) :: scalar_out
    integer(int64), intent(in)  :: scalarv_in(:)
    integer :: ierr
    ASSERT(size(scalarv_in) >= merge(npe,0,is_iop))
    call MPI_Scatter(scalarv_in, 1, MPI_INTEGER8, scalar_out, 1, MPI_INTEGER8, root, comm, ierr)
  end subroutine

  module subroutine dist_r4_0(scalar_out, scalarv_in)
    real(real32), intent(out) :: scalar_out
    real(real32), intent(in)  :: scalarv_in(:)
    integer :: ierr
    ASSERT(size(scalarv_in) >= merge(npe,0,is_iop))
    call MPI_Scatter(scalarv_in, 1, MPI_REAL4, scalar_out, 1, MPI_REAL4, root, comm, ierr)
  end subroutine

  module subroutine dist_r8_0(scalar_out, scalarv_in)
    real(real64), intent(out) :: scalar_out
    real(real64), intent(in)  :: scalarv_in(:)
    integer :: ierr
    ASSERT(size(scalarv_in) >= merge(npe,0,is_iop))
    call MPI_Scatter(scalarv_in, 1, MPI_REAL8, scalar_out, 1, MPI_REAL8, root, comm, ierr)
  end subroutine

  module subroutine dist_log_0(scalar_out, scalarv_in)
    logical, intent(out) :: scalar_out
    logical, intent(in)  :: scalarv_in(:)
    integer :: ierr
    ASSERT(size(scalarv_in) >= merge(npe,0,is_iop))
    call MPI_Scatter(scalarv_in, 1, MPI_LOGICAL, scalar_out, 1, MPI_LOGICAL, root, comm, ierr)
  end subroutine

!!!! DISTRIBUTE RANK-1 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine dist_i1_1(vector_out, vector_in)

    integer(int8), intent(out) :: vector_out(:)
    integer(int8), intent(in)  :: vector_in(:)

    integer :: inlen, outlen, ierr, j
    integer, allocatable :: counts(:), displs(:)

    outlen = size(vector_out)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(outlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    inlen = sum(counts)
    call MPI_Bcast(inlen, 1, MPI_INTEGER, root, comm, ierr)
    if (inlen == 0) return ! nothing to do

    ASSERT(size(vector_in) >= merge(inlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Scatterv(vector_in, counts, displs, MPI_INTEGER1, &
        vector_out, outlen, MPI_INTEGER1, root, comm, ierr)

  end subroutine dist_i1_1


  module subroutine dist_i4_1(vector_out, vector_in)

    integer(int32), intent(out) :: vector_out(:)
    integer(int32), intent(in)  :: vector_in(:)

    integer :: inlen, outlen, ierr, j
    integer, allocatable :: counts(:), displs(:)

    outlen = size(vector_out)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(outlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    inlen = sum(counts)
    call MPI_Bcast(inlen, 1, MPI_INTEGER, root, comm, ierr)
    if (inlen == 0) return ! nothing to do

    ASSERT(size(vector_in) >= merge(inlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Scatterv(vector_in, counts, displs, MPI_INTEGER4, &
        vector_out, outlen, MPI_INTEGER4, root, comm, ierr)

  end subroutine dist_i4_1


  module subroutine dist_i8_1(vector_out, vector_in)

    integer(int64), intent(out) :: vector_out(:)
    integer(int64), intent(in)  :: vector_in(:)

    integer :: inlen, outlen, ierr, j
    integer, allocatable :: counts(:), displs(:)

    outlen = size(vector_out)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(outlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    inlen = sum(counts)
    call MPI_Bcast(inlen, 1, MPI_INTEGER, root, comm, ierr)
    if (inlen == 0) return ! nothing to do

    ASSERT(size(vector_in) >= merge(inlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Scatterv(vector_in, counts, displs, MPI_INTEGER8, &
        vector_out, outlen, MPI_INTEGER8, root, comm, ierr)

  end subroutine dist_i8_1


  module subroutine dist_r4_1(vector_out, vector_in)

    real(real32), intent(out) :: vector_out(:)
    real(real32), intent(in)  :: vector_in(:)

    integer :: inlen, outlen, ierr, j
    integer, allocatable :: counts(:), displs(:)

    outlen = size(vector_out)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(outlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    inlen = sum(counts)
    call MPI_Bcast(inlen, 1, MPI_INTEGER, root, comm, ierr)
    if (inlen == 0) return ! nothing to do

    ASSERT(size(vector_in) >= merge(inlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Scatterv(vector_in, counts, displs, MPI_REAL4, &
        vector_out, outlen, MPI_REAL4, root, comm, ierr)

  end subroutine dist_r4_1


  module subroutine dist_r8_1(vector_out, vector_in)

    real(real64), intent(out) :: vector_out(:)
    real(real64), intent(in)  :: vector_in(:)

    integer :: inlen, outlen, ierr, j
    integer, allocatable :: counts(:), displs(:)

    outlen = size(vector_out)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(outlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    inlen = sum(counts)
    call MPI_Bcast(inlen, 1, MPI_INTEGER, root, comm, ierr)
    if (inlen == 0) return ! nothing to do

    ASSERT(size(vector_in) >= merge(inlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Scatterv(vector_in, counts, displs, MPI_REAL8, &
        vector_out, outlen, MPI_REAL8, root, comm, ierr)

  end subroutine dist_r8_1


  module subroutine dist_log_1(vector_out, vector_in)

    logical, intent(out) :: vector_out(:)
    logical, intent(in)  :: vector_in(:)

    integer :: inlen, outlen, ierr, j
    integer, allocatable :: counts(:), displs(:)

    outlen = size(vector_out)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(outlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    inlen = sum(counts)
    call MPI_Bcast(inlen, 1, MPI_INTEGER, root, comm, ierr)
    if (inlen == 0) return ! nothing to do

    ASSERT(size(vector_in) >= merge(inlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Scatterv(vector_in, counts, displs, MPI_LOGICAL, &
        vector_out, size(vector_out), MPI_LOGICAL, root, comm, ierr)

  end subroutine dist_log_1

!!!! DISTRIBUTE RANK-2 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine dist_i1_2(vector_out, vector_in)

    integer(int8), intent(out) :: vector_out(:,:)
    integer(int8), intent(in)  :: vector_in(:,:)

    integer :: inlen, outlen, ierr, j, block_size
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    outlen = size(vector_out, dim=2)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(outlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    inlen = sum(counts)
    call MPI_Bcast(inlen, 1, MPI_INTEGER, root, comm, ierr)
    if (inlen == 0) return ! nothing to do

    ASSERT(size(vector_in,dim=2) >= merge(inlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    if (is_iop) block_size = size(vector_in(:,1))
    call MPI_Bcast(block_size, 1, MPI_INTEGER, root, comm, ierr)
    call MPI_Type_contiguous(block_size, MPI_INTEGER1, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Scatterv(vector_in, counts, displs, block_type, &
        vector_out, outlen, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine dist_i1_2


  module subroutine dist_i4_2(vector_out, vector_in)

    integer(int32), intent(out) :: vector_out(:,:)
    integer(int32), intent(in)  :: vector_in(:,:)

    integer :: inlen, outlen, ierr, j, block_size
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    outlen = size(vector_out, dim=2)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(outlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    inlen = sum(counts)
    call MPI_Bcast(inlen, 1, MPI_INTEGER, root, comm, ierr)
    if (inlen == 0) return ! nothing to do

    ASSERT(size(vector_in,dim=2) >= merge(inlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    if (is_iop) block_size = size(vector_in(:,1))
    call MPI_Bcast(block_size, 1, MPI_INTEGER, root, comm, ierr)
    call MPI_Type_contiguous(block_size, MPI_INTEGER4, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Scatterv(vector_in, counts, displs, block_type, &
        vector_out, outlen, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine dist_i4_2


  module subroutine dist_i8_2(vector_out, vector_in)

    integer(int64), intent(out) :: vector_out(:,:)
    integer(int64), intent(in)  :: vector_in(:,:)

    integer :: inlen, outlen, ierr, j, block_size
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    outlen = size(vector_out, dim=2)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(outlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    inlen = sum(counts)
    call MPI_Bcast(inlen, 1, MPI_INTEGER, root, comm, ierr)
    if (inlen == 0) return ! nothing to do

    ASSERT(size(vector_in,dim=2) >= merge(inlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    if (is_iop) block_size = size(vector_in(:,1))
    call MPI_Bcast(block_size, 1, MPI_INTEGER, root, comm, ierr)
    call MPI_Type_contiguous(block_size, MPI_INTEGER8, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Scatterv(vector_in, counts, displs, block_type, &
        vector_out, outlen, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine dist_i8_2


  module subroutine dist_r4_2(vector_out, vector_in)

    real(real32), intent(out) :: vector_out(:,:)
    real(real32), intent(in)  :: vector_in(:,:)

    integer :: inlen, outlen, ierr, j, block_size
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    outlen = size(vector_out, dim=2)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(outlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    inlen = sum(counts)
    call MPI_Bcast(inlen, 1, MPI_INTEGER, root, comm, ierr)
    if (inlen == 0) return ! nothing to do

    ASSERT(size(vector_in,dim=2) >= merge(inlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    if (is_iop) block_size = size(vector_in(:,1))
    call MPI_Bcast(block_size, 1, MPI_INTEGER, root, comm, ierr)
    call MPI_Type_contiguous(block_size, MPI_REAL4, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Scatterv(vector_in, counts, displs, block_type, &
        vector_out, outlen, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine dist_r4_2


  module subroutine dist_r8_2(vector_out, vector_in)

    real(real64), intent(out) :: vector_out(:,:)
    real(real64), intent(in)  :: vector_in(:,:)

    integer :: inlen, outlen, ierr, j, block_size
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    outlen = size(vector_out, dim=2)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(outlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    inlen = sum(counts)
    call MPI_Bcast(inlen, 1, MPI_INTEGER, root, comm, ierr)
    if (inlen == 0) return ! nothing to do

    ASSERT(size(vector_in,dim=2) >= merge(inlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    if (is_iop) block_size = size(vector_in(:,1))
    call MPI_Bcast(block_size, 1, MPI_INTEGER, root, comm, ierr)
    call MPI_Type_contiguous(block_size, MPI_REAL8, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Scatterv(vector_in, counts, displs, block_type, &
        vector_out, outlen, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine dist_r8_2


  module subroutine dist_log_2(vector_out, vector_in)

    logical, intent(out) :: vector_out(:,:)
    logical, intent(in)  :: vector_in(:,:)

    integer :: inlen, outlen, ierr, j, block_size
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    outlen = size(vector_out, dim=2)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(outlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    inlen = sum(counts)
    call MPI_Bcast(inlen, 1, MPI_INTEGER, root, comm, ierr)
    if (inlen == 0) return ! nothing to do

    ASSERT(size(vector_in,dim=2) >= merge(inlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    if (is_iop) block_size = size(vector_in(:,1))
    call MPI_Bcast(block_size, 1, MPI_INTEGER, root, comm, ierr)
    call MPI_Type_contiguous(block_size, MPI_LOGICAL, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Scatterv(vector_in, counts, displs, block_type, &
        vector_out, outlen, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine dist_log_2

!!!! DISTRIBUTE RANK-3 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine dist_i1_3(vector_out, vector_in)

    integer(int8), intent(out) :: vector_out(:,:,:)
    integer(int8), intent(in)  :: vector_in(:,:,:)

    integer :: inlen, outlen, ierr, j, block_size
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    outlen = size(vector_out, dim=3)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(outlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    inlen = sum(counts)
    call MPI_Bcast(inlen, 1, MPI_INTEGER, root, comm, ierr)
    if (inlen == 0) return ! nothing to do

    ASSERT(size(vector_in,dim=3) >= merge(inlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    if (is_iop) block_size = size(vector_in(:,:,1))
    call MPI_Bcast(block_size, 1, MPI_INTEGER, root, comm, ierr)
    call MPI_Type_contiguous(block_size, MPI_INTEGER1, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Scatterv(vector_in, counts, displs, block_type, &
        vector_out, outlen, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine dist_i1_3

  module subroutine dist_i4_3(vector_out, vector_in)

    integer(int32), intent(out) :: vector_out(:,:,:)
    integer(int32), intent(in)  :: vector_in(:,:,:)

    integer :: inlen, outlen, ierr, j, block_size
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    outlen = size(vector_out, dim=3)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(outlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    inlen = sum(counts)
    call MPI_Bcast(inlen, 1, MPI_INTEGER, root, comm, ierr)
    if (inlen == 0) return ! nothing to do

#ifndef NDEBUG
    ASSERT(size(vector_in,dim=3) >= merge(inlen,0,is_iop))
#endif

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    if (is_iop) block_size = size(vector_in(:,:,1))
    call MPI_Bcast(block_size, 1, MPI_INTEGER, root, comm, ierr)
    call MPI_Type_contiguous(block_size, MPI_INTEGER4, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Scatterv(vector_in, counts, displs, block_type, &
        vector_out, outlen, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine dist_i4_3

  module subroutine dist_i8_3(vector_out, vector_in)

    integer(int64), intent(out) :: vector_out(:,:,:)
    integer(int64), intent(in)  :: vector_in(:,:,:)

    integer :: inlen, outlen, ierr, j, block_size
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    outlen = size(vector_out, dim=3)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(outlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    inlen = sum(counts)
    call MPI_Bcast(inlen, 1, MPI_INTEGER, root, comm, ierr)
    if (inlen == 0) return ! nothing to do

    ASSERT(size(vector_in,dim=3) >= merge(inlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    if (is_iop) block_size = size(vector_in(:,:,1))
    call MPI_Bcast(block_size, 1, MPI_INTEGER, root, comm, ierr)
    call MPI_Type_contiguous(block_size, MPI_INTEGER8, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Scatterv(vector_in, counts, displs, block_type, &
        vector_out, outlen, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine dist_i8_3

  module subroutine dist_r4_3(vector_out, vector_in)

    real(real32), intent(out) :: vector_out(:,:,:)
    real(real32), intent(in)  :: vector_in(:,:,:)

    integer :: inlen, outlen, ierr, j, block_size
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    outlen = size(vector_out, dim=3)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(outlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    inlen = sum(counts)
    call MPI_Bcast(inlen, 1, MPI_INTEGER, root, comm, ierr)
    if (inlen == 0) return ! nothing to do

    ASSERT(size(vector_in,dim=3) >= merge(inlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    if (is_iop) block_size = size(vector_in(:,:,1))
    call MPI_Bcast(block_size, 1, MPI_INTEGER, root, comm, ierr)
    call MPI_Type_contiguous(block_size, MPI_REAL4, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Scatterv(vector_in, counts, displs, block_type, &
        vector_out, outlen, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine dist_r4_3

  module subroutine dist_r8_3(vector_out, vector_in)

    real(real64), intent(out) :: vector_out(:,:,:)
    real(real64), intent(in)  :: vector_in(:,:,:)

    integer :: inlen, outlen, ierr, j, block_size
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    outlen = size(vector_out, dim=3)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(outlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    inlen = sum(counts)
    call MPI_Bcast(inlen, 1, MPI_INTEGER, root, comm, ierr)
    if (inlen == 0) return ! nothing to do

    ASSERT(size(vector_in,dim=3) >= merge(inlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    if (is_iop) block_size = size(vector_in(:,:,1))
    call MPI_Bcast(block_size, 1, MPI_INTEGER, root, comm, ierr)
    call MPI_Type_contiguous(block_size, MPI_REAL8, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Scatterv(vector_in, counts, displs, block_type, &
        vector_out, outlen, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine dist_r8_3

  module subroutine dist_log_3(vector_out, vector_in)

    logical, intent(out) :: vector_out(:,:,:)
    logical, intent(in)  :: vector_in(:,:,:)

    integer :: inlen, outlen, ierr, j, block_size
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    outlen = size(vector_out, dim=3)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(outlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    inlen = sum(counts)
    call MPI_Bcast(inlen, 1, MPI_INTEGER, root, comm, ierr)
    if (inlen == 0) return ! nothing to do

    ASSERT(size(vector_in,dim=3) >= merge(inlen,0,is_iop))

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    if (is_iop) block_size = size(vector_in(:,:,1))
    call MPI_Bcast(block_size, 1, MPI_INTEGER, root, comm, ierr)
    call MPI_Type_contiguous(block_size, MPI_LOGICAL, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Scatterv(vector_in, counts, displs, block_type, &
        vector_out, outlen, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine dist_log_3

end submodule distribute_impl
