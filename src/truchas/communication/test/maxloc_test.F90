!! Unit Tests for PARALLEL_COMMUNICATION Ad Hoc Maxloc/Minloc Procedures
!!
!! Copyright 2022 Neil N. Carlson <neil.n.carlson@gmail.com>
!! Use subject to the MIT license: https://opensource.org/licenses/MIT
!!

program maxloc_test

  use parallel_communication
  use,intrinsic :: iso_fortran_env
  use mpi
  implicit none

  integer :: ierr, status, global_status

  call MPI_Init(ierr)
  call init_parallel_communication

  if (npe /= 3) then
    if (is_IOP) write(error_unit,'(a)') 'test must be run with 3 ranks'
    call MPI_Finalize(ierr)
    stop 1
  end if

  status = 0
  call test_a0
  call test_a1
  call test_a2
  call test_b0
  call test_b1
  call test_b2
  call test_b3
  call test_b4

  call MPI_Allreduce(status, global_status, 1, MPI_INTEGER, MPI_MAX, comm, ierr)

  call MPI_Finalize(ierr)

  if (global_status /= 0) stop 1

contains

  subroutine write_fail (errmsg)
    use,intrinsic :: iso_fortran_env, only: error_unit
    character(*), intent(in) :: errmsg
    status = 1
    write(error_unit,'(a)') 'Failed: ' // errmsg
  end subroutine

  ! Corner case: 0-sized array on every process should return 0
  subroutine test_a0
    real(real64), allocatable :: a(:)
    allocate(a(0))
    if (any(global_maxloc(a) /= 0)) call write_fail('test_a0')
  end subroutine

  ! Corner case: Initial 0-sized array and all values of remaining arrays equal
  ! to the value returned by maxloc applied to a 0-sized array should return 1
  subroutine test_a1
    real(real64), allocatable :: a(:)
    allocate(a(2*(this_pe-1)))
    a = maxloc([real(real64) :: ], dim=1)
    if (any(global_maxloc(a) /= 1)) call write_fail('test_a1')
  end subroutine

  ! Generic case. Should return the global index of the *first* element that
  ! attains the maximum value
  subroutine test_a2
    real(real64), allocatable :: a(:)
    select case (this_pe)
    case (1)
      allocate(a(1))
      a(1) = 1
    case (2)
      allocate(a(3))
      a(1) = 2
      a(2) = 3
      a(3) = 5
    case (3)
      allocate(a(2))
      a(1) = 5
      a(2) = 4
    end select
    if (any(global_maxloc(a) /= 4)) call write_fail('test_a2')
  end subroutine

  ! Corner case: 0-sized array on every process should return 0
  subroutine test_b0
    real(real64), allocatable :: a(:)
    integer :: pid, lindex
    allocate(a(0))
    call global_maxloc_sub(a, pid, lindex)
    if (pid /= 0 .or. lindex /= 0) call write_fail('test_b0')
  end subroutine

  ! Corner case: Initial 0-sized array and all values of remaining arrays equal
  ! to the value returned by maxloc applied to a 0-sized array should return 2
  ! for the pid and 1 for lindex
  subroutine test_b1
    real(real64), allocatable :: a(:)
    integer :: pid, lindex
    allocate(a(2*(this_pe-1)))
    a = maxloc([real(real64) :: ], dim=1)
    call global_maxloc_sub(a, pid, lindex)
    if (pid /= 2 .or. lindex /= 1) call write_fail('test_b1')
  end subroutine

  ! Generic case. Should return the global index of the *first* element that
  ! attains the maximum value
  subroutine test_b2
    real(real64), allocatable :: a(:)
    integer :: pid, lindex
    select case (this_pe)
    case (1)
      allocate(a(1))
      a(1) = 1
    case (2)
      allocate(a(3))
      a(1) = 2
      a(2) = 3
      a(3) = 5
    case (3)
      allocate(a(2))
      a(1) = 5
      a(2) = 4
    end select
    call global_maxloc_sub(a, pid, lindex)
    if (pid /= 2 .or. lindex /= 3) call write_fail('test_b2')
  end subroutine

  ! Generic case with mask. Should return the global index of the
  ! *second* element that attains the maximum value as the mask
  ! excludes the first such element.
  subroutine test_b3
    real(real64), allocatable :: a(:)
    logical, allocatable :: mask(:)
    integer :: pid, lindex
    select case (this_pe)
    case (1)
      allocate(a(1), mask(1))
      a(1) = 1
      mask = .true.
    case (2)
      allocate(a(3), mask(3))
      a(1) = 2
      a(2) = 3
      a(3) = 5
      mask = .true.
      mask(3) = .false.
    case (3)
      allocate(a(2), mask(2))
      a(1) = 5
      a(2) = 4
      mask = .true.
    end select
    call global_maxloc_sub(a, pid, lindex, mask)
    if (pid /= 3 .or. lindex /= 1) call write_fail('test_b3')
  end subroutine

  ! Corner case with mask. All elements are ignored. Should return the same
  ! results as for all 0-sized arrays: pid = 0, lindex = 0
  subroutine test_b4
    real(real64), allocatable :: a(:)
    logical, allocatable :: mask(:)
    integer :: pid, lindex
    select case (this_pe)
    case (1)
      allocate(a(1), mask(1))
      a(1) = 1
    case (2)
      allocate(a(3), mask(3))
      a(1) = 2
      a(2) = 3
      a(3) = 5
    case (3)
      allocate(a(2), mask(2))
      a(1) = 5
      a(2) = 4
    end select
    mask = .false.
    call global_maxloc_sub(a, pid, lindex, mask)
    if (pid /= 0 .or. lindex /= 0) call write_fail('test_b4')
  end subroutine

end program
