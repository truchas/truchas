!! VF_DATA_TYPE
!!
!! A container for the raw distributed view factor data with an initialization
!! method that reads the data from an enclosure radiation disk file.
!!
!! David Neill-Asanza <dhna@lanl.gov>, July 2019
!! Neil N. Carlson <nnc@lanl.gov>, refactoring June 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! The container stores the view factor matrix PHI as a distributed CSR matrix,
!! and the ambient view factors as a distibuted vector. The CSR matrix is
!! block-row partitioned and distributed: NPATCH is the number of local rows,
!! NPATCH_TOT is the total number of rows (and columns). Here a "patch" may be
!! a single enclosure face or a collection of such faces -- the distinction is
!! irrelevant to this container. The public (read-only) data components are:
!!
!!    NPATCH, NPATCH_TOT
!!    AMB_VF(:) -- distributed patch-based vector of ambient view factors.
!!      If there are no ambient view factors this array is not allocated.
!!
!! The public type-bound methods are:
!!
!!    INIT(FILE) initializes the object, reading the view factor data from a
!!    radiation enclosure disk file connected to the ENCL_RAD_FILE object FILE.
!!    The rows are equidistributed across the processors.
!!
!!    MATVEC(X, Y) computes the matrix-vector product Y = PHI*X, for face-based
!!    vectors X and Y.
!!

!TODO: global_x component to reduce allocations

#include "f90_assert.fpp"

module vf_data_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use rad_encl_file_type
  use parallel_communication
  implicit none
  private

  type, public :: vf_data
    private
    integer, public :: npatch, npatch_tot
    integer, allocatable :: ia(:), ja(:)
    real, allocatable :: vf(:)
    real, allocatable, public :: amb_vf(:)
    logical, public :: has_ambient
  contains
    procedure :: init
    procedure :: matvec
  end type

contains

  subroutine init(this, file)

    class(vf_data), intent(out) :: this
    type(rad_encl_file), intent(in) :: file

    integer :: nface_tot

    if (is_IOP) call file%get_patch_dims(nface_tot, this%npatch_tot)
    call broadcast(this%npatch_tot)

    !! Block-row partition of the view factor matrix
    this%npatch = this%npatch_tot/nPE
    if (this_PE <= modulo(this%npatch_tot,nPE)) this%npatch = 1 + this%npatch

    call load_view_factors(this, file, this%npatch, this%npatch_tot)

  end subroutine init

  !! Computes the matrix-vector product Y = PHI*X.
  subroutine matvec(this, x, y)

    class(vf_data), intent(in) :: this
    real(r8), intent(in)  :: x(:)
    real(r8), intent(out) :: y(:)

    real(r8) :: tmp, global_x(this%npatch_tot)
    integer :: i, j

    ASSERT(size(x) == this%npatch)
    ASSERT(size(y) == this%npatch)

    call collate(global_x, x)
    call broadcast(global_x)

    do i = 1, this%npatch
      tmp = 0
      do j = this%ia(i), this%ia(i+1)-1
        tmp = tmp + this%vf(j) * global_x(this%ja(j))
      end do
      y(i) = tmp
    end do

  end subroutine matvec

  !! This auxiliary subroutine reads the view factor data from FILE and
  !! distributes it in a block pattern. NROWS is the number of rows owned by
  !! this process, and NROWS_TOT is the total number of rows (and columns)
  !! of the view factor matrix and the size of the ambient view factor vector.

  subroutine load_view_factors(this, file, nrows, nrows_tot)

    use,intrinsic :: iso_fortran_env, only: i8 => int64

    class(vf_data), intent(inout) :: this
    type(rad_encl_file), intent(in) :: file
    integer, intent(in) :: nrows      ! number of rows on this process
    integer, intent(in) :: nrows_tot  ! total number of rows (and columns)

    integer(i8) :: start
    integer :: j, n
    integer :: vf_bsize(nPE), lengths(nPE), idum0(0)
    real :: rdum0(0)
    integer, allocatable :: ibuf(:)
    real, allocatable :: rbuf(:)

    !! Read and distribute the ambient view factors.
    if (is_IOP) this%has_ambient = file%has_ambient()
    call broadcast(this%has_ambient)

    if (this%has_ambient) then
      allocate(rbuf(merge(nrows_tot, 0, is_IOP)))
      if (is_IOP) call file%get_ambient(rbuf)
      allocate(this%amb_vf(nrows))
      call distribute(this%amb_vf, rbuf)
      deallocate(rbuf)
    end if

    !! Read and distribute the VF matrix nonzero row counts.
    allocate(ibuf(merge(nrows_tot, 0, is_IOP)))
    if (is_IOP) call file%get_vf_rowcount(ibuf)
    allocate(this%ia(nrows+1))
    call distribute(this%ia(2:), ibuf)
    deallocate(ibuf)

    !! Convert the row counts into the local IA indexing array.
    this%ia(1) = 1
    do j = 2, size(this%ia)
      this%ia(j) = this%ia(j) + this%ia(j-1)
    end do

    !! Determine the sizes of the distributed VF matrix.
    n = this%ia(nrows+1) - this%ia(1)
    allocate(this%vf(n), this%ja(n))
    call collate(vf_bsize, n)

    !! Read the VF matrix in process-sized blocks, sending them to the owning
    !! processes as we go.  PGSLib only provides the distribute collective to
    !! do this, and so we distribute with 0-sized destination arrays for all
    !! processes but the receiving one.  We also use the optional LENGTHS
    !! argument to distribute (only referenced on the IO process) that gives
    !! the number of items to distribute to each process, resulting in the
    !! actual sizes of the array arguments being ignored except to check that
    !! they are sufficiently large.  This allows us to simplify the calls.
    !! The code is structured for future use of MPI_Isend/Irecv, abandoning
    !! the usual 'single code path' pattern.

    if (is_IOP) then
      call file%get_vf_rows(this%vf, this%ja, start=1_i8)
      if (nPE > 1) then
        n = maxval(vf_bsize(2:))
        allocate(ibuf(n), rbuf(n))
        lengths = 0
        start = 1 + vf_bsize(1)
        do n = 2, nPE
          call file%get_vf_rows(rbuf(:vf_bsize(n)), ibuf(:vf_bsize(n)), start)
          lengths(n) = vf_bsize(n)
          call distribute(rdum0, rbuf, lengths)
          call distribute(idum0, ibuf, lengths)
          lengths(n) = 0
          start = start + vf_bsize(n)
        end do
        deallocate(ibuf, rbuf)
      end if
    else  ! everybody else just participates in the distribute calls.
      do n = 2, nPE
        if (n == this_PE) then
          call distribute(this%vf, rdum0, lengths)
          call distribute(this%ja, idum0, lengths)
        else
          call distribute(rdum0, rdum0, lengths)
          call distribute(idum0, idum0, lengths)
        end if
      end do
    end if

  end subroutine load_view_factors

end module vf_data_type
