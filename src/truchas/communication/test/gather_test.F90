!! Unit Tests for PARALLEL_COMMUNICATION Gather Procedures
!!
!! Copyright 2022 Neil N. Carlson <neil.n.carlson@gmail.com>
!! Use subject to the MIT license: https://opensource.org/licenses/MIT
!!

program main

  use mpi
  use parallel_communication
  use,intrinsic :: iso_fortran_env
  implicit none

  integer :: ierr, status, global_status
  logical :: pass ! local to tests

  call MPI_Init(ierr)
  call init_parallel_communication

  status = 0
  call gath_scalar
  call gath_rank1
  call gath_rank1_zero
  call gath_rank2
  call gath_rank2_zero
  call gath_array_section
  call gath_log_scalar
  call gath_log_rank1
  call gath_log_rank2
  call gath_char_scalar
  call gath_char_rank1
  call gath_char_rank2

  call MPI_Allreduce(status, global_status, 1, MPI_INTEGER, MPI_MAX, comm, ierr)

  call MPI_Finalize(ierr)

  if (global_status /= 0) stop 1

contains

  subroutine write_result(pass, name)
    logical, intent(in) :: pass
    character(*), intent(in) :: name
    if (global_all(pass)) then
      if (is_IOP) write(output_unit,'(a)') 'Passed: ' // name
    else
      status = 1
      if (is_IOP) write(output_unit,'(a)') 'FAILED: ' // name
    end if
  end subroutine

  ! Distribute a scalar value to each process
  subroutine gath_scalar
    integer, allocatable :: adest(:)
    integer :: n, asrc
    asrc = this_pe
    if (is_IOP) then
      adest = [(n, n=1,npe)]
    else
      allocate(adest(0))
    end if
    block
      integer(int8) :: dest(size(adest)), src
      src = asrc
      dest = 0
      call gather(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'gath_scalar_int8')
    end block
    block
      integer(int32) :: dest(size(adest)), src
      src = asrc
      dest = 0
      call gather(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'gath_scalar_int32')
    end block
    block
      integer(int64) :: dest(size(adest)), src
      src = asrc
      dest = 0
      call gather(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'gath_scalar_int64')
    end block
    block
      real(real32) :: dest(size(adest)), src
      src = asrc
      dest = 0
      call gather(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'gath_scalar_real32')
    end block
    block
      real(real64) :: dest(size(adest)), src
      src = asrc
      dest = 0
      call gather(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'gath_scalar_real64')
    end block
    block
      complex(real32) :: dest(size(adest)), src
      src = (1,-1)*asrc
      dest = 0
      call gather(src, dest)
      pass = all(dest == (1,-1)*adest)
      call write_result(pass, 'gath_scalar_complex32')
    end block
    block
      complex(real64) :: dest(size(adest)), src
      src = (1,-1)*asrc
      dest = 0
      call gather(src, dest)
      pass = all(dest == (1,-1)*adest)
      call write_result(pass, 'gath_scalar_complex64')
    end block
  end subroutine

  ! generic rank-1 vector case
  subroutine gath_rank1
    integer, allocatable :: asrc(:), adest(:)
    integer :: j, n
    n = (this_pe*(this_pe-1))/2
    asrc = [(n+j, j=1,this_pe)]
    if (is_IOP) then
      adest = [(j, j=1, (npe*(npe+1))/2)]
    else
      allocate(adest(0))
    end if
    block
      integer(int8), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call gather(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'gath_rank1_int8')
    end block
    block
      integer(int32), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call gather(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'gath_rank1_int32')
    end block
    block
      integer(int64), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call gather(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'gath_rank1_int64')
    end block
    block
      real(real32), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call gather(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'gath_rank1_real32')
    end block
    block
      real(real64), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call gather(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'gath_rank1_real64')
    end block
    block
      complex(real32), allocatable :: src(:), dest(:)
      src = (1,-1)*asrc
      allocate(dest(size(adest)))
      call gather(src, dest)
      pass = all(dest == (1,-1)*adest)
      call write_result(pass, 'gath_rank1_complex32')
    end block
    block
      complex(real64), allocatable :: src(:), dest(:)
      src = (1,-1)*asrc
      allocate(dest(size(adest)))
      call gather(src, dest)
      pass = all(dest == (1,-1)*adest)
      call write_result(pass, 'gath_rank1_complex64')
    end block
  end subroutine

  ! Rank-1 vector case with a 0-sized vector
  subroutine gath_rank1_zero
    integer, allocatable :: asrc(:), adest(:)
    integer :: j, n
    n = ((this_pe-1)*(this_pe-2))/2
    asrc = [(n+j, j=1,this_pe-1)]
    if (is_IOP) then
      adest = [(j, j=1, (npe*(npe-1))/2)]
    else
      allocate(adest(0))
    end if
    block
      integer(int8), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call gather(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'gath_rank1_zero_int8')
    end block
    block
      integer(int32), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call gather(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'gath_rank1_zero_int32')
    end block
    block
      integer(int64), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call gather(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'gath_rank1_zero_int64')
    end block
    block
      real(real32), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call gather(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'gath_rank1_zero_real32')
    end block
    block
      real(real64), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call gather(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'gath_rank1_zero_real64')
    end block
    block
      complex(real32), allocatable :: src(:), dest(:)
      src = (1,-1)*asrc
      allocate(dest(size(adest)))
      call gather(src, dest)
      pass = all(dest == (1,-1)*adest)
      call write_result(pass, 'gath_rank1_zero_complex32')
    end block
    block
      complex(real64), allocatable :: src(:), dest(:)
      src = (1,-1)*asrc
      allocate(dest(size(adest)))
      call gather(src, dest)
      pass = all(dest == (1,-1)*adest)
      call write_result(pass, 'gath_rank1_zero_complex64')
    end block
  end subroutine

  ! generic rank-2 vector case
  subroutine gath_rank2
    integer, allocatable :: asrc(:,:), adest(:,:)
    integer :: j, n
    n = (this_pe*(this_pe-1))/2
    allocate(asrc(2,this_pe))
    asrc(1,:) = [(n+j, j=1,this_pe)]
    asrc(2,:) = -asrc(1,:)
    if (is_IOP) then
      allocate(adest(2,(npe*(npe+1))/2))
      adest(1,:) = [(j, j=1, (npe*(npe+1))/2)]
      adest(2,:) = -adest(1,:)
    else
      allocate(adest(2,0))
    end if
    block
      integer(int8), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call gather(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'gath_rank2_int8')
    end block
    block
      integer(int32), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call gather(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'gath_rank2_int32')
    end block
    block
      integer(int64), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call gather(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'gath_rank2_int64')
    end block
    block
      real(real32), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call gather(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'gath_rank2_real32')
    end block
    block
      real(real64), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call gather(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'gath_rank2_real64')
    end block
    block
      complex(real32), allocatable :: src(:,:), dest(:,:)
      src = (1,-1)*asrc
      allocate(dest(2,size(adest,2)))
      call gather(src, dest)
      pass = all(dest == (1,-1)*adest)
      call write_result(pass, 'gath_rank2_complex32')
    end block
    block
      complex(real64), allocatable :: src(:,:), dest(:,:)
      src = (1,-1)*asrc
      allocate(dest(2,size(adest,2)))
      call gather(src, dest)
      pass = all(dest == (1,-1)*adest)
      call write_result(pass, 'gath_rank2_complex64')
    end block
  end subroutine

  ! Rank-1 vector case with a 0-sized vector
  subroutine gath_rank2_zero
    integer, allocatable :: asrc(:,:), adest(:,:)
    integer :: j, n
    n = ((this_pe-1)*(this_pe-2))/2
    allocate(asrc(2,this_pe-1))
    asrc(1,:) = [(n+j, j=1,this_pe-1)]
    asrc(2,:) = -asrc(1,:)
    if (is_IOP) then
      n = (npe*(npe-1))/2
      allocate(adest(2,n))
      adest(1,:) = [(j, j=1,n)]
      adest(2,:) = -adest(1,:)
    else
      allocate(adest(2,0))
    end if
    block
      integer(int8), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call gather(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'gath_rank2_zero_int8')
    end block
    block
      integer(int32), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call gather(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'gath_rank2_zero_int32')
    end block
    block
      integer(int64), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call gather(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'gath_rank2_zero_int64')
    end block
    block
      real(real32), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call gather(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'gath_rank2_zero_real32')
    end block
    block
      real(real64), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call gather(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'gath_rank2_zero_real64')
    end block
    block
      complex(real32), allocatable :: src(:,:), dest(:,:)
      src = (1,-1)*asrc
      allocate(dest(2,size(adest,2)))
      call gather(src, dest)
      pass = all(dest == (1,-1)*adest)
      call write_result(pass, 'gath_rank2_zero_complex32')
    end block
    block
      complex(real64), allocatable :: src(:,:), dest(:,:)
      src = (1,-1)*asrc
      allocate(dest(2,size(adest,2)))
      call gather(src, dest)
      pass = all(dest == (1,-1)*adest)
      call write_result(pass, 'gath_rank2_zero_complex64')
    end block
  end subroutine

  ! generic rank-2 vector case
  subroutine gath_array_section
    integer, allocatable :: asrc(:,:), adest(:,:)
    integer :: j, n
    n = (this_pe*(this_pe-1))/2
    allocate(asrc(3,2*this_pe), source=-1)
    asrc(1,1::2) = [(n+j, j=1,this_pe)]
    asrc(3,1::2) = -asrc(1,1::2)
    if (is_IOP) then
      n = (npe*(npe+1))/2
      allocate(adest(3,2*n), source=0)
      adest(1,1::2) = [(j, j=1,n)]
      adest(3,1::2) = -adest(1,1::2)
    else
      allocate(adest(3,0))
    end if
    block
      integer(int8), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(3,size(adest,2)), source=0_int8)
      call gather(src(1::2,1::2), dest(1::2,1::2))
      pass = all(dest == adest)
      call write_result(pass, 'gath_array_section_int8')
    end block
    block
      integer(int32), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(3,size(adest,2)), source=0_int32)
      call gather(src(1::2,1::2), dest(1::2,1::2))
      pass = all(dest == adest)
      call write_result(pass, 'gath_array_section_int32')
    end block
    block
      integer(int64), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(3,size(adest,2)), source=0_int64)
      call gather(src(1::2,1::2), dest(1::2,1::2))
      pass = all(dest == adest)
      call write_result(pass, 'gath_array_section_int64')
    end block
    block
      real(real32), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(3,size(adest,2)), source=0.0_real32)
      call gather(src(1::2,1::2), dest(1::2,1::2))
      pass = all(dest == adest)
      call write_result(pass, 'gath_array_section_real32')
    end block
    block
      real(real64), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(3,size(adest,2)), source=0.0_real64)
      call gather(src(1::2,1::2), dest(1::2,1::2))
      pass = all(dest == adest)
      call write_result(pass, 'gath_array_section_real64')
    end block
    block
      complex(real32), allocatable :: src(:,:), dest(:,:)
      src =(1,-1)* asrc
      allocate(dest(3,size(adest,2)), source=cmplx(0,kind=real32))
      call gather(src(1::2,1::2), dest(1::2,1::2))
      pass = all(dest == (1,-1)*adest)
      call write_result(pass, 'gath_array_section_complex32')
    end block
    block
      complex(real64), allocatable :: src(:,:), dest(:,:)
      src = (1,-1)*asrc
      allocate(dest(3,size(adest,2)), source=cmplx(0,kind=real64))
      call gather(src(1::2,1::2), dest(1::2,1::2))
      pass = all(dest == (1,-1)*adest)
      call write_result(pass, 'gath_array_section_complex64')
    end block
  end subroutine

  subroutine gath_log_scalar
    logical, allocatable :: dest(:)
    logical :: src
    src = .true.
    if (is_IOP) then
      allocate(dest(npe), source=.false.)
    else
      allocate(dest(0))
    end if
    call gather(src, dest)
    pass = all(dest)
    call write_result(pass, 'gath_log_scalar')
  end subroutine

  subroutine gath_log_rank1
    logical, allocatable :: src(:), dest(:)
    integer :: n
    allocate(src(this_pe), source=.true.)
    if (is_IOP) then
      n = (npe*(npe+1))/2
      allocate(dest(n), source=.false.)
    else
      allocate(dest(0))
    end if
    call gather(src, dest)
    pass = all(dest)
    call write_result(pass, 'gath_log_rank1')
  end subroutine

  subroutine gath_log_rank2
    logical, allocatable :: src(:,:), dest(:,:)
    integer :: n
    allocate(src(2,this_pe))
    src(1,:) = .true.
    src(2,:) = .false.
    if (is_IOP) then
      n = (npe*(npe+1))/2
      allocate(dest(2,n))
      dest(1,:) = .false.
      dest(2,:) = .true.
    else
      allocate(dest(2,0))
    end if
    call gather(src, dest)
    pass = all(dest(1,:)) .and. all(.not.dest(2,:))
    call write_result(pass, 'gath_log_rank2')
  end subroutine

  subroutine gath_char_scalar
    character(3), allocatable :: adest(:), dest(:)
    character(3) :: src
    integer :: n
    write(src,'(i3)') this_pe
    if (is_IOP) then
      allocate(dest(npe), adest(npe))
      do n = 1, npe
        write(adest(n),'(i3)') n
      end do
    else
      allocate(dest(0), adest(0))
    end if
    dest = ''
    call gather(src, dest)
    pass = all(dest == adest)
    call write_result(pass, 'gath_char_scalar')
  end subroutine

  subroutine gath_char_rank1
    character(3), allocatable :: src(:), adest(:), dest(:)
    integer :: j, n
    n = (this_pe*(this_pe-1))/2
    allocate(src(this_pe))
    do j = 1, this_pe
      write(src(j),'(i3)') n+j
    end do
    if (is_IOP) then
      n = (npe*(npe+1))/2
      allocate(dest(n), adest(n))
      do j = 1, n
        write(adest(j),'(i3)') j
      end do
    else
      allocate(dest(0), adest(0))
    end if
    dest = ''
    call gather(src, dest)
    pass = all(dest == adest)
    call write_result(pass, 'gath_char_rank1')
  end subroutine

  subroutine gath_char_rank2
    character(3), allocatable :: src(:,:), adest(:,:), dest(:,:)
    integer :: j, n
    n = (this_pe*(this_pe-1))/2
    allocate(src(2,this_pe))
    do j = 1, this_pe
      write(src(1,j),'(i3)') n+j
      write(src(2,j),'(i3.3)') n+j
    end do
    if (is_IOP) then
      n = (npe*(npe+1))/2
      allocate(dest(2,n), adest(2,n))
      do j = 1, n
        write(adest(1,j),'(i3)') j
        write(adest(2,j),'(i3.3)') j
      end do
    else
      allocate(dest(2,0), adest(2,0))
    end if
    dest = ''
    call gather(src, dest)
    pass = all(dest == adest)
    call write_result(pass, 'gath_char_rank2')
  end subroutine

end program
