!! Unit Tests for INDEX_MAP Collate Procedures
!!
!! Copyright 2022 Neil N. Carlson <neil.n.carlson@gmail.com>
!! Use subject to the MIT license: https://opensource.org/licenses/MIT
!!

program main

  use mpi
  use index_map_type
  use,intrinsic :: iso_fortran_env
  implicit none

  integer :: ierr, my_rank, nproc, status
  logical :: is_root

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  is_root = (my_rank == 0)

  status = 0
  call coll_rank1
  call coll_rank1_zero
  call coll_rank2
  call coll_rank2_zero
  call coll_array_section
  call coll_rank3
  call coll_rank3_zero
  call coll_log_rank1
  call coll_log_rank2
  call coll_log_rank3

  call MPI_Finalize(ierr)
  if (status /= 0) error stop 1

contains

  subroutine write_result(pass, name)
    logical, value :: pass
    character(*), intent(in) :: name
#ifdef AVOID_MPI_IN_PLACE
    block
      logical :: pass_loc
      pass_loc = pass
      call MPI_Allreduce(pass_loc, pass, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    end block
#else
    call MPI_Allreduce(MPI_IN_PLACE, pass, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
#endif
    if (pass) then
      if (is_root) write(output_unit,'(a)') 'Passed: ' //  name
    else
      status = 1
      if (is_root) write(output_unit,'(a)') 'FAILED: ' //  name
    end if
  end subroutine

  ! generic rank-1 array case
  subroutine coll_rank1
    type(index_map) :: imap
    integer, allocatable :: asrc(:), adest(:)
    integer :: j
    call imap%init(1+my_rank)
    if (is_root) then
      adest = [(j, j=1, imap%global_size)]
    else
      allocate(adest(0))
    end if
    asrc = [(j, j=imap%first_gid,imap%last_gid)]
    block
      integer(int32), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call imap%gather(src, dest)
      call write_result(all(dest == adest), 'coll_rank1_int32')
    end block
    block
      integer(int64), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call imap%gather(src, dest)
      call write_result(all(dest == adest), 'coll_rank1_int64')
    end block
    block
      real(real32), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call imap%gather(src, dest)
      call write_result(all(dest == adest), 'coll_rank1_real32')
    end block
    block
      real(real64), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call imap%gather(src, dest)
      call write_result(all(dest == adest), 'coll_rank1_real64')
    end block
    block
      complex(real32), allocatable :: src(:), dest(:)
      src = (1,-1)*asrc
      allocate(dest(size(adest)))
      call imap%gather(src, dest)
      call write_result(all(dest == (1,-1)*adest), 'coll_rank1_complex32')
    end block
    block
      complex(real64), allocatable :: src(:), dest(:)
      src = (1,-1)*asrc
      allocate(dest(size(adest)))
      call imap%gather(src, dest)
      call write_result(all(dest == (1,-1)*adest), 'coll_rank1_complex64')
    end block
  end subroutine

  ! Rank-1 array case with a 0-sized vector
  subroutine coll_rank1_zero
    type(index_map) :: imap
    integer, allocatable :: asrc(:), adest(:)
    integer :: j
    call imap%init(modulo(my_rank+2,nproc))
    if (is_root) then
      adest = [(j, j=1, imap%global_size)]
    else
      allocate(adest(0))
    end if
    asrc = [(j, j=imap%first_gid,imap%last_gid)]
    block
      integer(int32), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call imap%gather(src, dest)
      call write_result(all(dest == adest), 'coll_rank1_zero_int32')
    end block
    block
      integer(int64), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call imap%gather(src, dest)
      call write_result(all(dest == adest), 'coll_rank1_zero_int64')
    end block
    block
      real(real32), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call imap%gather(src, dest)
      call write_result(all(dest == adest), 'coll_rank1_zero_real32')
    end block
    block
      real(real64), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call imap%gather(src, dest)
      call write_result(all(dest == adest), 'coll_rank1_zero_real64')
    end block
    block
      complex(real32), allocatable :: src(:), dest(:)
      src = (1,-1)*asrc
      allocate(dest(size(adest)))
      call imap%gather(src, dest)
      call write_result(all(dest == (1,-1)*adest), 'coll_rank1_zero_complex32')
    end block
    block
      complex(real64), allocatable :: src(:), dest(:)
      src = (1,-1)*asrc
      allocate(dest(size(adest)))
      call imap%gather(src, dest)
      call write_result(all(dest == (1,-1)*adest), 'coll_rank1_zero_complex64')
    end block
  end subroutine

  ! generic rank-2 array case
  subroutine coll_rank2
    type(index_map) :: imap
    integer, allocatable :: asrc(:,:), adest(:,:)
    integer :: j
    call imap%init(1+my_rank)
    if (is_root) then
      allocate(adest(2,imap%global_size))
      adest(1,:) = [(j, j=1, imap%global_size)]
      adest(2,:) = -adest(1,:)
    else
      allocate(adest(2,0))
    end if
    allocate(asrc(2,imap%onp_size))
    asrc(1,:) = [(j, j=imap%first_gid,imap%last_gid)]
    asrc(2,:) = -asrc(1,:)
    block
      integer(int32), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call imap%gather(src, dest)
      call write_result(all(dest == adest), 'coll_rank2_int32')
    end block
    block
      integer(int64), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call imap%gather(src, dest)
      call write_result(all(dest == adest), 'coll_rank2_int64')
    end block
    block
      real(real32), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call imap%gather(src, dest)
      call write_result(all(dest == adest), 'coll_rank2_real32')
    end block
    block
      real(real64), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call imap%gather(src, dest)
      call write_result(all(dest == adest), 'coll_rank2_real64')
    end block
    block
      complex(real32), allocatable :: src(:,:), dest(:,:)
      src = (1,-1)*asrc
      allocate(dest(2,size(adest,2)))
      call imap%gather(src, dest)
      call write_result(all(dest == (1,-1)*adest), 'coll_rank2_complex32')
    end block
    block
      complex(real64), allocatable :: src(:,:), dest(:,:)
      src = (1,-1)*asrc
      allocate(dest(2,size(adest,2)))
      call imap%gather(src, dest)
      call write_result(all(dest == (1,-1)*adest), 'coll_rank2_complex64')
    end block
  end subroutine

  ! Rank-2 array case with a 0-sized vector
  subroutine coll_rank2_zero
    type(index_map) :: imap
    integer, allocatable :: asrc(:,:), adest(:,:)
    integer :: j
    call imap%init(modulo(my_rank+nproc-1,nproc))
    if (is_root) then
      allocate(adest(2,imap%global_size))
      adest(1,:) = [(j, j=1, imap%global_size)]
      adest(2,:) = -adest(1,:)
    else
      allocate(adest(2,0))
    end if
    allocate(asrc(2,imap%onp_size))
    asrc(1,:) = [(j, j=imap%first_gid,imap%last_gid)]
    asrc(2,:) = -asrc(1,:)
    block
      integer(int32), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call imap%gather(src, dest)
      call write_result(all(dest == adest), 'coll_rank2_zero_int32')
    end block
    block
      integer(int64), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call imap%gather(src, dest)
      call write_result(all(dest == adest), 'coll_rank2_zero_int64')
    end block
    block
      real(real32), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call imap%gather(src, dest)
      call write_result(all(dest == adest), 'coll_rank2_zero_real32')
    end block
    block
      real(real64), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call imap%gather(src, dest)
      call write_result(all(dest == adest), 'coll_rank2_zero_real64')
    end block
    block
      complex(real32), allocatable :: src(:,:), dest(:,:)
      src = (1,-1)*asrc
      allocate(dest(2,size(adest,2)))
      call imap%gather(src, dest)
      call write_result(all(dest == (1,-1)*adest), 'coll_rank2_zero_complex32')
    end block
    block
      complex(real64), allocatable :: src(:,:), dest(:,:)
      src = (1,-1)*asrc
      allocate(dest(2,size(adest,2)))
      call imap%gather(src, dest)
      call write_result(all(dest == (1,-1)*adest), 'coll_rank2_zero_complex64')
    end block
  end subroutine

  ! rank-2 array section case
  subroutine coll_array_section
    type(index_map) :: imap
    integer, allocatable :: asrc(:,:), adest(:,:)
    integer :: j
    call imap%init(1+my_rank)
    if (is_root) then
      allocate(adest(3,2*imap%global_size), source=0)
      adest(1,1::2) = [(j, j=1, imap%global_size)]
      adest(3,1::2) = -adest(1,1::2)
    else
      allocate(adest(3,0))
    end if
    allocate(asrc(3,2*imap%onp_size), source=127)
    asrc(1,1::2) = [(j, j=imap%first_gid,imap%last_gid)]
    asrc(3,1::2) = -asrc(1,1::2)
    block
      integer(int32), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(3,size(adest,2)), source=0_int32)
      call imap%gather(src(1::2,1::2), dest(1::2,1::2))
      call write_result(all(dest == adest), 'coll_array_section_int32')
    end block
    block
      integer(int64), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(3,size(adest,2)), source=0_int64)
      call imap%gather(src(1::2,1::2), dest(1::2,1::2))
      call write_result(all(dest == adest), 'coll_array_section_int64')
    end block
    block
      real(real32), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(3,size(adest,2)), source=0.0_real32)
      call imap%gather(src(1::2,1::2), dest(1::2,1::2))
      call write_result(all(dest == adest), 'coll_array_section_real32')
    end block
    block
      real(real64), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(3,size(adest,2)), source=0.0_real64)
      call imap%gather(src(1::2,1::2), dest(1::2,1::2))
      call write_result(all(dest == adest), 'coll_array_section_real64')
    end block
    block
      complex(real32), allocatable :: src(:,:), dest(:,:)
      src = (1,-1)*asrc
      allocate(dest(3,size(adest,2)), source=(0.0_real32,0.0_real32))
      call imap%gather(src(1::2,1::2), dest(1::2,1::2))
      call write_result(all(dest == (1,-1)*adest), 'coll_array_section_complex32')
    end block
    block
      complex(real64), allocatable :: src(:,:), dest(:,:)
      src = (1,-1)*asrc
      allocate(dest(3,size(adest,2)), source=(0.0_real64,0.0_real64))
      call imap%gather(src(1::2,1::2), dest(1::2,1::2))
      call write_result(all(dest == (1,-1)*adest), 'coll_array_section_complex64')
    end block
  end subroutine

  ! generic rank-3 array case
  subroutine coll_rank3
    type(index_map) :: imap
    integer, allocatable :: asrc(:,:,:), adest(:,:,:)
    integer :: j
    call imap%init(1+my_rank)
    if (is_root) then
      allocate(adest(2,2,imap%global_size))
      adest(1,1,:) = [(j, j=1, imap%global_size)]
      adest(2,1,:) = -adest(1,1,:)
      adest(:,2,:) = 2*adest(:,1,:)
    else
      allocate(adest(2,2,0))
    end if
    allocate(asrc(2,2,imap%onp_size))
    asrc(1,1,:) = [(j, j=imap%first_gid,imap%last_gid)]
    asrc(2,1,:) = -asrc(1,1,:)
    asrc(:,2,:) = 2*asrc(:,1,:)
    block
      integer(int32), allocatable :: src(:,:,:), dest(:,:,:)
      src = asrc
      allocate(dest(2,2,size(adest,3)))
      call imap%gather(src, dest)
      call write_result(all(dest == adest), 'coll_rank3_int32')
    end block
    block
      integer(int64), allocatable :: src(:,:,:), dest(:,:,:)
      src = asrc
      allocate(dest(2,2,size(adest,3)))
      call imap%gather(src, dest)
      call write_result(all(dest == adest), 'coll_rank3_int64')
    end block
    block
      real(real32), allocatable :: src(:,:,:), dest(:,:,:)
      src = asrc
      allocate(dest(2,2,size(adest,3)))
      call imap%gather(src, dest)
      call write_result(all(dest == adest), 'coll_rank3_real32')
    end block
    block
      real(real64), allocatable :: src(:,:,:), dest(:,:,:)
      src = asrc
      allocate(dest(2,2,size(adest,3)))
      call imap%gather(src, dest)
      call write_result(all(dest == adest), 'coll_rank3_real64')
    end block
    block
      complex(real32), allocatable :: src(:,:,:), dest(:,:,:)
      src = (1,-1)*asrc
      allocate(dest(2,2,size(adest,3)))
      call imap%gather(src, dest)
      call write_result(all(dest == (1,-1)*adest), 'coll_rank3_complex32')
    end block
    block
      complex(real64), allocatable :: src(:,:,:), dest(:,:,:)
      src = (1,-1)*asrc
      allocate(dest(2,2,size(adest,3)))
      call imap%gather(src, dest)
      call write_result(all(dest == (1,-1)*adest), 'coll_rank3_complex64')
    end block
  end subroutine

  ! Rank-3 array case with a 0-sized vector
  subroutine coll_rank3_zero
    type(index_map) :: imap
    integer, allocatable :: asrc(:,:,:), adest(:,:,:)
    integer :: j
    call imap%init(modulo(my_rank+nproc-1,nproc))
    if (is_root) then
      allocate(adest(2,2,imap%global_size))
      adest(1,1,:) = [(j, j=1, imap%global_size)]
      adest(2,1,:) = -adest(1,1,:)
      adest(:,2,:) = 2*adest(:,1,:)
    else
      allocate(adest(2,2,0))
    end if
    allocate(asrc(2,2,imap%onp_size))
    asrc(1,1,:) = [(j, j=imap%first_gid,imap%last_gid)]
    asrc(2,1,:) = -asrc(1,1,:)
    asrc(:,2,:) = 2*asrc(:,1,:)
    block
      integer(int32), allocatable :: src(:,:,:), dest(:,:,:)
      src = asrc
      allocate(dest(2,2,size(adest,3)))
      call imap%gather(src, dest)
      call write_result(all(dest == adest), 'coll_rank3_zero_int32')
    end block
    block
      integer(int64), allocatable :: src(:,:,:), dest(:,:,:)
      src = asrc
      allocate(dest(2,2,size(adest,3)))
      call imap%gather(src, dest)
      call write_result(all(dest == adest), 'coll_rank3_zero_int64')
    end block
    block
      real(real32), allocatable :: src(:,:,:), dest(:,:,:)
      src = asrc
      allocate(dest(2,2,size(adest,3)))
      call imap%gather(src, dest)
      call write_result(all(dest == adest), 'coll_rank3_zero_real32')
    end block
    block
      real(real64), allocatable :: src(:,:,:), dest(:,:,:)
      src = asrc
      allocate(dest(2,2,size(adest,3)))
      call imap%gather(src, dest)
      call write_result(all(dest == adest), 'coll_rank3_zero_real64')
    end block
    block
      complex(real32), allocatable :: src(:,:,:), dest(:,:,:)
      src = (1,-1)*asrc
      allocate(dest(2,2,size(adest,3)))
      call imap%gather(src, dest)
      call write_result(all(dest == (1,-1)*adest), 'coll_rank3_zero_complex32')
    end block
    block
      complex(real64), allocatable :: src(:,:,:), dest(:,:,:)
      src = (1,-1)*asrc
      allocate(dest(2,2,size(adest,3)))
      call imap%gather(src, dest)
      call write_result(all(dest == (1,-1)*adest), 'coll_rank3_zero_complex64')
    end block
  end subroutine

  subroutine coll_log_rank1
    type(index_map) :: imap
    logical, allocatable :: src(:), dest(:)
    call imap%init(1+my_rank)
    if (is_root) then
      allocate(dest(imap%global_size), source=.false.)
    else
      allocate(dest(0))
    end if
    allocate(src(imap%onp_size), source=.true.)
    call imap%gather(src, dest)
    call write_result(all(dest), 'coll_log_rank1')
  end subroutine

  subroutine coll_log_rank2
    type(index_map) :: imap
    logical, allocatable :: src(:,:), dest(:,:)
    call imap%init(1+my_rank)
    if (is_root) then
      allocate(dest(2,imap%global_size))
      dest(1,:) = .false.
      dest(2,:) = .true.
    else
      allocate(dest(2,0))
    end if
    allocate(src(2,imap%onp_size))
    src(1,:) = .true.
    src(2,:) = .false.
    call imap%gather(src, dest)
    call write_result(all(dest(1,:) .and. .not.dest(2,:)), 'coll_log_rank2')
  end subroutine

  subroutine coll_log_rank3
    type(index_map) :: imap
    logical, allocatable :: src(:,:,:), dest(:,:,:)
    call imap%init(1+my_rank)
    if (is_root) then
      allocate(dest(2,2,imap%global_size))
      dest(1,1,:) = .false.
      dest(2,1,:) = .true.
      dest(:,2,:) = .not.dest(:,1,:)
    else
      allocate(dest(2,2,0))
    end if
    allocate(src(2,2,imap%onp_size))
    src(1,1,:) = .true.
    src(2,1,:) = .false.
    src(:,2,:) = .not.src(:,1,:)
    call imap%gather(src, dest)
    call write_result(all(dest(1,1,:) .and. .not.dest(2,1,:) .and. .not.dest(1,2,:) .and. dest(2,2,:)), 'coll_log_rank3')
  end subroutine

end program
