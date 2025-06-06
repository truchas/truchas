!!
!! DMUMPS_SERIAL_SOLVER_TYPE
!!
!! This module provides an interface to the serial version of the
!! double-precision direct multifrontal solver in the MUMPS library. The serial
!! version allows for more output metrics.
!!
!! Zach Jibben <zjibben@lanl.gov>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module dmumps_serial_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use index_map_type
  use,intrinsic :: iso_fortran_env, only: int64, real64
  implicit none
  include 'dmumps_struc.h'
  private

  type, public :: dmumps_serial_solver
    private
    type(dmumps_struc), pointer :: mumps => null()
    type(index_map), pointer :: imap => null() ! potentially unowned reference
    logical :: imap_owned = .false.
    integer :: nrow, bsize
  contains
    procedure :: init
    procedure :: setup
    procedure :: solve
    final :: mumps_solver_delete
  end type dmumps_serial_solver

contains

  subroutine mumps_solver_delete(this)
    type(dmumps_serial_solver), intent(inout) :: this
    if (associated(this%mumps)) then
      this%mumps%job = -2
      call dmumps(this%mumps)
      if (associated(this%mumps%irn)) deallocate(this%mumps%irn)
      if (associated(this%mumps%jcn)) deallocate(this%mumps%jcn)
      if (associated(this%mumps%A)) deallocate(this%mumps%A)
      if (associated(this%mumps%rhs)) deallocate(this%mumps%rhs)
      deallocate(this%mumps)
    end if
    if (this%imap_owned) deallocate(this%imap)
  end subroutine mumps_solver_delete


  ! Symmetry must be set at initialization.
  !   symmetry = 0 --> unsymmetric
  !   symmetry = 1 --> symmetric positive-definite
  !   symmetry = 2 --> general symmetric
  subroutine init(this, symmetry, stat)

    use parallel_communication, only: comm, npe
    use truchas_logging_services

    class(dmumps_serial_solver), intent(out) :: this
    integer, intent(in) :: symmetry
    integer, intent(out) :: stat

    allocate(this%mumps)

    INSIST(npe == 1)

    ! initialize mumps instance
    this%mumps%comm = comm
    this%mumps%par = 1
    this%mumps%sym = symmetry
    this%mumps%job = -1
    this%mumps%icntl(4) = 4 ! verbosity
    call dmumps(this%mumps)
    stat = this%mumps%infog(1)

  end subroutine init


  ! Provide the matrix and perform analysis & factorization.
  subroutine setup(this, A, stat)

    use parallel_communication, only: global_sum
    use pcsr_matrix_type

    class(dmumps_serial_solver), intent(inout) :: this
    type(pcsr_matrix), intent(in) :: A
    integer, intent(out) :: stat

    integer :: i, ii

    INSIST(associated(this%mumps))

    if (.not.associated(this%imap)) this%imap => A%graph%row_imap

    this%mumps%icntl(11) = 1 ! compute norms

    ! general input settings
    this%mumps%icntl(5) = 0 ! assembled format
    this%mumps%icntl(18) = 0 ! serial

    ! settings
    this%mumps%icntl(24) = 1 ! null pivot rows should be computed
    this%mumps%icntl(25) = 0 ! don't compute the null space basis
    this%mumps%icntl(35) = 0 ! enable BLR feature (0 = off, 1 = auto, 2 = on)
    this%mumps%icntl(20) = 0 ! serial RHS
    this%mumps%icntl(21) = 0 ! serial solution
    !this%mumps%icntl(14) = 150 ! allow 150% workspace increase

    ! matrix, RHS, and sol specification
    this%mumps%N = A%graph%row_imap%global_size
    this%mumps%NNZ = global_sum(size(A%values))
    this%mumps%nrhs = 1
    this%mumps%lrhs = A%nrow_onp

    associate(nnz => this%mumps%NNZ)
      allocate(this%mumps%irn(nnz), this%mumps%jcn(nnz), this%mumps%A(nnz), &
          this%mumps%rhs(this%mumps%lrhs))
    end associate

    do i = 1, A%nrow_onP
      this%mumps%irn(A%graph%xadj(i):A%graph%xadj(i+1)-1) = A%graph%row_imap%global_index(i)
    end do
    this%mumps%jcn(:) = A%graph%col_imap%global_index(A%graph%adjncy)
    this%mumps%A(:) = A%values

    ! perform analysis & factorization
    this%mumps%job = 4
    print *, "START SETUP"
    call dmumps(this%mumps)
    print *, "DONE SETUP"
    stat = this%mumps%infog(1)
    print *, '1 setup pcsr', stat, this%mumps%info(23), A%nrow_onP, A%nrow

  end subroutine setup


  subroutine solve(this, b, x, stat)

    class(dmumps_serial_solver), intent(inout) :: this
    real(r8), intent(in) :: b(:)
    real(r8), intent(out) :: x(:)
    integer, intent(out) :: stat

    integer :: xi, ig, il

    INSIST(associated(this%mumps))
    ASSERT(size(b) >= this%mumps%nloc_rhs)
    ASSERT(size(x) >= this%mumps%lsol_loc)

    this%mumps%rhs(:) = b(:this%mumps%lrhs)

    this%mumps%job = 3
    print *, "START SOLVE"
    call dmumps(this%mumps)
    print *, "DONE SOLVE"
    stat = this%mumps%infog(1)

    do xi = 1, this%mumps%lrhs
      x(xi) = this%mumps%rhs(xi)
    end do

  end subroutine solve

end module dmumps_serial_solver_type
