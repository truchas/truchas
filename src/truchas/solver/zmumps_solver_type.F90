!!
!! ZMUMPS_SOLVER_TYPE
!!
!! This module provides an interface to the double-precision complex direct
!! multifrontal solver in the MUMPS library.
!!
!! Zach Jibben <zjibben@lanl.gov>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module zmumps_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use index_map_type
  use zmumps_serial_solver_type
  implicit none
#ifdef USE_MUMPS
  include 'zmumps_struc.h'
#endif
  private

  type, public :: zmumps_solver
    private
    type(zmumps_serial_solver) :: serial
#ifdef USE_MUMPS
    type(zmumps_struc), pointer :: mumps => null()
#endif
    type(index_map), pointer :: imap => null() ! unowned reference
    integer :: nrow, bsize
    integer :: verbosity = 0
  contains
    procedure :: init
    procedure :: setup
    procedure :: solve
    final :: zmumps_solver_delete
  end type zmumps_solver

contains

  subroutine zmumps_solver_delete(this)
    type(zmumps_solver), intent(inout) :: this
#ifdef USE_MUMPS
    if (associated(this%mumps)) then
      this%mumps%job = -2
      call zmumps(this%mumps)
      if (associated(this%mumps%irn_loc)) deallocate(this%mumps%irn_loc)
      if (associated(this%mumps%jcn_loc)) deallocate(this%mumps%jcn_loc)
      if (associated(this%mumps%a_loc)) deallocate(this%mumps%a_loc)
      if (associated(this%mumps%irhs_loc)) deallocate(this%mumps%irhs_loc)
      if (associated(this%mumps%rhs_loc)) deallocate(this%mumps%rhs_loc)
      if (associated(this%mumps%isol_loc)) deallocate(this%mumps%isol_loc)
      if (associated(this%mumps%sol_loc)) deallocate(this%mumps%sol_loc)
      deallocate(this%mumps)
    end if
#endif
  end subroutine zmumps_solver_delete


  ! Symmetry must be set at initialization.
  !   symmetry = 0 --> unsymmetric
  !   symmetry = 1 --> symmetric positive-definite
  !   symmetry = 2 --> general symmetric
  subroutine init(this, symmetry, stat, verbosity)

    use parallel_communication, only: comm
    use truchas_logging_services

    class(zmumps_solver), intent(out) :: this
    integer, intent(in) :: symmetry
    integer, intent(out) :: stat
    integer, intent(in), optional :: verbosity

#ifndef USE_MUMPS
    call tls_fatal("MUMPS requested, but Truchas wasn't compiled with MUMPS support.")
#else
    if (present(verbosity)) this%verbosity = verbosity

    if (this%verbosity > 2) then
      call this%serial%init(symmetry, stat)
      return
    end if

    allocate(this%mumps)
    this%mumps%irn_loc => null()
    this%mumps%jcn_loc => null()
    this%mumps%a_loc => null()
    this%mumps%irhs_loc => null()
    this%mumps%rhs_loc => null()
    this%mumps%isol_loc => null()

    ! initialize mumps instance
    this%mumps%comm = comm
    this%mumps%par = 1
    this%mumps%sym = symmetry
    this%mumps%job = -1

    call zmumps(this%mumps)
    stat = this%mumps%infog(1)

    select case (this%verbosity)
    case(0)
      this%mumps%icntl(2) = 0
      this%mumps%icntl(3) = 0
      this%mumps%icntl(4) = 1
      this%mumps%icntl(11) = 0
    case(1)
      this%mumps%icntl(4) = 2
      this%mumps%icntl(11) = 2
    case(2)
      this%mumps%icntl(4) = 4
      this%mumps%icntl(11) = 1
    case default
      INSIST(.false.)
    end select
#endif

  end subroutine init


  ! Provide the matrix and perform analysis & factorization.
  subroutine setup(this, A, stat)

    use parallel_communication, only: global_sum
    use complex_pcsr_matrix_type

    class(zmumps_solver), intent(inout) :: this
    type(complex_pcsr_matrix), intent(in) :: A
    integer, intent(out) :: stat

    integer :: i, nnz_loc

#ifdef USE_MUMPS
    if (this%verbosity > 2) then
      call this%serial%setup(A, stat)
      return
    end if

    INSIST(associated(this%mumps))

    if (.not.associated(this%imap)) this%imap => A%graph%row_imap

    ! general input settings
    this%mumps%icntl(5) = 0 ! assembled format
    this%mumps%icntl(18) = 3 ! parallel distribution

    ! settings
    ! this%mumps%icntl(24) = 1 ! null pivot rows should be computed
    ! this%mumps%icntl(25) = 0 ! don't compute the null space basis
    this%mumps%icntl(35) = 0 ! enable BLR feature (0 = off, 1 = auto, 2 = on)
    this%mumps%icntl(20) = 11 ! distributed RHS
    this%mumps%icntl(21) = 1 ! distributed solution
    this%mumps%icntl(14) = 30 ! increase work space array storage by 30%

    ! matrix, RHS, and sol specification
    nnz_loc = A%graph%xadj(A%nrow_onP+1)-1
    this%mumps%N = A%graph%row_imap%global_size
    this%mumps%NNZ = global_sum(nnz_loc)
    this%mumps%NNZ_loc = nnz_loc
    this%mumps%nrhs = 1
    this%mumps%nloc_rhs = A%nrow_onp
    this%mumps%lrhs_loc = A%nrow_onp
    this%mumps%lsol_loc = A%nrow_onp

    allocate(this%mumps%IRN_loc(nnz_loc), this%mumps%JCN_loc(nnz_loc), this%mumps%A_loc(nnz_loc), &
        this%mumps%irhs_loc(this%mumps%nloc_rhs), this%mumps%rhs_loc(this%mumps%nloc_rhs))

    do i = 1, A%nrow_onP
      this%mumps%IRN_loc(A%graph%xadj(i):A%graph%xadj(i+1)-1) = A%graph%row_imap%global_index(i)
    end do
    this%mumps%JCN_loc(:) = A%graph%col_imap%global_index(A%graph%adjncy(:nnz_loc))
    this%mumps%A_loc(:) = A%values(:nnz_loc)

    ! perform analysis & factorization
    this%mumps%job = 4
    call zmumps(this%mumps)
    stat = this%mumps%infog(1)

    ! rhs & solution specification
    do i = 1, this%mumps%nloc_rhs
      this%mumps%irhs_loc(i) = i + A%graph%row_imap%first_gid - 1
    end do
    !this%mumps%irhs_loc(:) = [ (i, i = A%graph%row_imap%first_gid, A%graph%row_imap%last_gid) ]
    this%mumps%lsol_loc = this%mumps%info(23)
    allocate(this%mumps%isol_loc(this%mumps%lsol_loc), this%mumps%sol_loc(this%mumps%lsol_loc))
#endif

  end subroutine setup


  subroutine solve(this, b, x, stat)

    class(zmumps_solver), intent(inout) :: this
    complex(r8), intent(in) :: b(:)
    complex(r8), intent(out) :: x(:)
    integer, intent(out) :: stat

    integer :: xi, ig, il
    character(64) :: msg

#ifdef USE_MUMPS
    if (this%verbosity > 2) then
      call this%serial%solve(b, x, stat)
      return
    end if

    INSIST(associated(this%mumps))
    ASSERT(size(b) >= this%mumps%nloc_rhs)
    !ASSERT(size(x) >= this%mumps%lsol_loc)

    this%mumps%rhs_loc(:) = b(:this%mumps%nloc_rhs)

    this%mumps%job = 3
    call zmumps(this%mumps)
    stat = this%mumps%infog(1)
    if (stat /= 0) return

    call gather_distributed(this%imap, this%mumps%isol_loc, this%mumps%sol_loc, x)

    ! block
    !   use truchas_logging_services
    !   write (msg,"(a,es13.3)") "MUMPS solve complete: ", this%mumps%rinfog(5)
    !   call tls_info(msg)
    ! end block
#endif

  end subroutine solve


  ! Gather unusually-distributed data onto the appropriate MPI ranks.
  subroutine gather_distributed(imap, data_in_gid, data_in, data_out)

    use mpi
    use parallel_communication, only: comm, this_pe, npe
    use sort_utilities, only: insertion_sort

    type(index_map), intent(in) :: imap
    integer, intent(in) :: data_in_gid(:)
    complex(r8), intent(in) :: data_in(:)
    complex(r8), intent(out) :: data_out(:)

    integer :: i, ierr
    integer, dimension(npe) :: first_gid, sendcounts, sdispls, recvcounts, rdispls
    integer :: all_sendcounts(npe,npe)
    integer :: owned_rank(size(data_in)), tmp(size(data_in)), isendbuf(size(data_in))
    complex(r8) :: rsendbuf(size(data_in))
    complex(r8), allocatable :: rrecvbuf(:)
    integer, allocatable :: irecvbuf(:)

    ASSERT(size(data_in) == size(data_in_gid))

    call MPI_Allgather(imap%first_gid, 1, MPI_INTEGER, first_gid, 1, MPI_INTEGER, comm, ierr)
    ASSERT(ierr == 0)

    ! set up send buffers
    sendcounts = 0
    do i = 1, size(data_in)
      ! TODO: bisectional search likely faster than findloc
      owned_rank(i) = findloc(first_gid <= data_in_gid(i), .true., back=.true., dim=1)
      INSIST(owned_rank(i) > 0 .and. owned_rank(i) <= npe)
      sendcounts(owned_rank(i)) = sendcounts(owned_rank(i)) + 1
      tmp(i) = owned_rank(i)
      rsendbuf(i) = data_in(i)
      isendbuf(i) = data_in_gid(i)
    end do
    sdispls = inclusive_scan(sendcounts)
    call insertion_sort(rsendbuf, tmp) ! TODO: quicksort faster at high processor counts
    call insertion_sort(isendbuf, owned_rank)

    ! set up recv buffers
    call MPI_Alltoall(sendcounts, 1, MPI_INTEGER, recvcounts, 1, MPI_INTEGER, comm, ierr)
    ASSERT(ierr == 0)
    rdispls = inclusive_scan(recvcounts)
    allocate(rrecvbuf(sum(recvcounts)), irecvbuf(sum(recvcounts)))

    ! communicate data
    ! TODO: use non-blocking calls
    call MPI_Alltoallv(rsendbuf, sendcounts, sdispls, MPI_COMPLEX16, &
        rrecvbuf, recvcounts, rdispls, MPI_COMPLEX16, comm, ierr)
    ASSERT(ierr == 0)
    call MPI_Alltoallv(isendbuf, sendcounts, sdispls, MPI_INTEGER, &
        irecvbuf, recvcounts, rdispls, MPI_INTEGER, comm, ierr)
    ASSERT(ierr == 0)

    ! populate output data
    do i = 1, size(rrecvbuf)
      ASSERT(irecvbuf(i) - imap%first_gid + 1 <= imap%onp_size)
      data_out(irecvbuf(i) - imap%first_gid + 1) = rrecvbuf(i)
    end do

  contains

    function inclusive_scan(data_in)
      integer, intent(in) :: data_in(:)
      integer :: inclusive_scan(size(data_in))
      integer :: i
      inclusive_scan(1) = 0
      do i = 2, size(data_in)
        inclusive_scan(i) = inclusive_scan(i-1) + data_in(i-1)
      end do
    end function inclusive_scan

  end subroutine gather_distributed

end module zmumps_solver_type
