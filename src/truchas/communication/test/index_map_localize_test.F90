!! Unit Tests for INDEX_MAP Localize Procedures
!!
!! Copyright 2022 Neil N. Carlson <neil.n.carlson@gmail.com>
!! Use subject to the MIT license: https://opensource.org/licenses/MIT
!!

program main

  use,intrinsic :: iso_fortran_env, only: output_unit
  use index_map_type
  use mpi
  implicit none

  logical :: is_root
  integer :: my_rank, nproc, ierr, status

  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  is_root = (my_rank == 0)

  status = 0
  call test_basic
  call test_rank1
  call test_rank2
  call test_rank1_offp
  call test_rank1_domain_offp
  call test_struct

  call MPI_Finalize(ierr)
  if (status /= 0) error stop 1

contains

  subroutine write_result(pass, name)
    logical, value :: pass
    character(*), intent(in) :: name
    call MPI_Allreduce(MPI_IN_PLACE, pass, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    if (pass) then
      if (is_root) write(output_unit,'(a)') 'Passed: ' //  name
    else
      status = 1
      if (is_root) write(output_unit,'(a)') 'FAILED: ' //  name
    end if
  end subroutine

  ! No off-process references

  subroutine test_basic

    type(index_map) :: domain, range
    integer, allocatable :: g_index(:), l_index(:)
    integer, allocatable :: g_input(:), l_input(:), g_output(:), l_output(:)
    integer :: j

    call domain%init(2)
    call range%init(4)

    g_index = [(2*j, j=1,merge(domain%global_size,0,is_root))]
    call domain%localize_index_array(g_index, range, l_index)
    call write_result(.not.allocated(range%offp_index), 'test_basic_no-offp')

    g_input = [(j, j=1,merge(range%global_size,0,is_root))]
    allocate(l_input(range%local_size))
    call range%scatter(g_input, l_input)

    g_output = g_input(g_index)
    allocate(l_output(domain%local_size))
    call domain%scatter(g_output, l_output)

    call write_result(all(l_output == l_input(l_index)), 'test_basic')

  end subroutine

  ! Rank-1 indexing array with off-process references
  ! No off-process in domain or range.
  ! Distributed index_map initialization.

  subroutine test_rank1

    type(index_map) :: domain, range
    integer, allocatable :: g_index(:), l_index(:)
    integer, allocatable :: g_input(:), l_input(:), g_output(:), l_output(:)
    integer :: j

    call domain%init(1+my_rank)
    call range%init(2*(nproc-my_rank))

    g_index = [(2*j, j=1,merge(domain%global_size,0,is_root))]
    call domain%localize_index_array(g_index, range, l_index)

    g_input = [(j, j=1,merge(range%global_size,0,is_root))]
    allocate(l_input(range%local_size))
    call range%scatter(g_input, l_input)
    call range%gather_offp(l_input)

    g_output = g_input(g_index)
    allocate(l_output(domain%local_size))
    call domain%scatter(g_output, l_output)

    call write_result(all(l_output == l_input(l_index)), 'test_rank1')

  end subroutine

  ! Rank-2 indexing array with off-process references
  ! No off-process in domain or range.
  ! Distributed index_map initialization.

  subroutine test_rank2

    type(index_map) :: domain, range
    integer, allocatable :: g_index(:,:), l_index(:,:)
    integer, allocatable :: g_input(:), l_input(:), g_output(:,:), l_output(:,:)
    integer :: j

    call domain%init(1+my_rank)
    call range%init(2*(nproc-my_rank))

    allocate(g_index(2,merge(domain%global_size,0,is_root)))
    g_index(1,:) = [(2*j, j=1,merge(domain%global_size,0,is_root))]
    g_index(2,:) = [(2*j-1, j=1,merge(domain%global_size,0,is_root))]
    call domain%localize_index_array(g_index, range, l_index)

    g_input = [(j, j=1,merge(range%global_size,0,is_root))]
    allocate(l_input(range%local_size))
    call range%scatter(g_input, l_input)
    call range%gather_offp(l_input)

    allocate(g_output, mold=g_index)
    g_output(1,:) = g_input(g_index(1,:))
    g_output(2,:) = g_input(g_index(2,:))
    allocate(l_output, mold=l_index)
    call domain%scatter(g_output, l_output)

    call write_result(all(l_output(1,:) == l_input(l_index(1,:)) .and. &
                          l_output(2,:) == l_input(l_index(2,:))), 'test_rank2')

  end subroutine

  ! Rank-1 indexing array with off-process references
  ! No off-process in domain, but range has pre-existing off-process
  ! Exercises distributed index_map initialization

  subroutine test_rank1_offp

    type(index_map) :: domain, range
    integer, allocatable :: g_index(:), l_index(:)
    integer, allocatable :: g_input(:), l_input(:), g_output(:), l_output(:)
    integer, allocatable :: offp_index(:)
    integer :: j, bsize

    bsize = 3
    offp_index = [1+modulo((my_rank+1)*bsize, nproc*bsize)]
    call range%init(bsize, offp_index)
    call domain%init(1)

    g_index = [(1+modulo(j,nproc)*bsize, j=1,merge(nproc,0,is_root))]
    call domain%localize_index_array(g_index, range, l_index)

    g_input = [(j, j=1,merge(range%global_size,0,is_root))]
    allocate(l_input(range%local_size))
    call range%scatter(g_input, l_input)
    call range%gather_offp(l_input)

    g_output = g_input(g_index)
    allocate(l_output(domain%local_size))
    call domain%scatter(g_output, l_output)

    call write_result(all(l_output == l_input(l_index)), 'test_rank1_offp')

  end subroutine

  ! Rank-1 indexing array with off-process references
  ! Domain has off-process elements; no off-process in range.
  ! Distributed and root index_map initialization.

  subroutine test_rank1_domain_offp

    type(index_map) :: domain, range
    integer, allocatable :: g_index(:), l_index(:)
    integer, allocatable :: g_input(:), l_input(:), g_output(:), l_output(:)
    integer, allocatable :: bsizes(:), offp_counts(:), offp_indices(:)
    integer :: j, n

    if (is_root) then
      allocate(bsizes(nproc), offp_counts(nproc), offp_indices(nproc))
      bsizes = [(j, j=1,nproc)]
      offp_counts = 1
      n = (nproc*(nproc+1))/2
      offp_indices = [(1+modulo((j*(j+1))/2,n), j=1,nproc)]
    else
      allocate(bsizes(0), offp_counts(0), offp_indices(0))
    end if

    call domain%init(bsizes, offp_counts, offp_indices)
    call range%init(2*(nproc-my_rank))

    g_index = [(2*j, j=1,merge(domain%global_size,0,is_root))]
    call domain%localize_index_array(g_index, range, l_index)

    g_input = [(j, j=1,merge(range%global_size,0,is_root))]
    allocate(l_input(range%local_size))
    call range%scatter(g_input, l_input)
    call range%gather_offp(l_input)

    g_output = g_input(g_index)
    allocate(l_output(domain%local_size))
    call domain%scatter(g_output, l_output)
    call domain%gather_offp(l_output)

    call write_result(all(l_output == l_input(l_index)), 'test_rank1_domain_offp')

  end subroutine

  ! Indexing array structure with off-process reference and
  ! domain with off-process

  subroutine test_struct

    type(index_map) :: domain, range, temp
    integer, allocatable :: g_count(:), g_index(:), l_count(:), l_index(:)
    integer, allocatable :: g_input(:), l_input(:), g_output(:), l_output(:)

    integer :: j, n, m

    m = (nproc*(nproc+1))/2
    n = 1 + modulo(((my_rank+1)*(my_rank+2))/2, m)
    call domain%init(my_rank+1, offp_index=[n])
    call range%init(nproc-my_rank)

    g_count = [(modulo(j,3), j=1,merge(domain%global_size,0,is_root))]
    g_index = [(1+modulo(j,range%global_size), j=1,sum(g_count))]
    call domain%localize_index_array(g_count, g_index, range, l_count, l_index)

    g_input = [(j, j=1,merge(range%global_size,0,is_root))]
    allocate(l_input(range%local_size))
    call range%scatter(g_input, l_input)
    call range%gather_offp(l_input)

    g_output = g_input(g_index)
    call temp%init(domain, g_count)
    allocate(l_output(temp%local_size))
    call temp%scatter(g_output, l_output)
    call temp%gather_offp(l_output)

    call write_result(all(l_output == l_input(l_index)), 'test_struct')

  end subroutine

end program
