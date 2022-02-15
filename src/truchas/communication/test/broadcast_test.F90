!! Unit Tests for PARALLEL_COMMUNICATION Broadcast Procedures
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
  call bcast_scalar
  call bcast_rank1
  call bcast_rank2
  call bcast_rank3
  call bcast_section
  call bcast_log_scalar
  call bcast_log_rank1
  call bcast_log_rank2
  call bcast_log_rank3
  call bcast_char_scalar
  call bcast_char_rank1
  call bcast_char_rank2
  call bcast_char_rank3

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

  ! scalar argument case
  subroutine bcast_scalar
    integer :: a
    a = 42
    block
      integer(int8) :: x
      x = merge(a, 0, is_iop)
      call broadcast(x)
      call write_result((x == a), 'bcast_scalar_int8')
    end block
    block
      integer(int32) :: x
      x = merge(a, 0, is_iop)
      call broadcast(x)
      call write_result((x == a), 'bcast_scalar_int32')
    end block
    block
      integer(int64) :: x
      x = merge(a, 0, is_iop)
      call broadcast(x)
      call write_result((x == a), 'bcast_scalar_int64')
    end block
    block
      real(real32) :: x
      x = merge(a, 0, is_iop)
      call broadcast(x)
      call write_result((x == a), 'bcast_scalar_real32')
    end block
    block
      real(real64) :: x
      x = merge(a, 0, is_iop)
      call broadcast(x)
      call write_result((x == a), 'bcast_scalar_real64')
    end block
  end subroutine

  ! rank-1 array case
  subroutine bcast_rank1
    integer :: a(2)
    a = [42,17]
    block
      integer(int8) :: x(size(a))
      x = merge(a, 0, is_iop)
      call broadcast(x)
      call write_result(all(x == a), 'bcast_rank1_int8')
    end block
    block
      integer(int32) :: x(size(a))
      x = merge(a, 0, is_iop)
      call broadcast(x)
      call write_result(all(x == a), 'bcast_rank1_int32')
    end block
    block
      integer(int64) :: x(size(a))
      x = merge(a, 0, is_iop)
      call broadcast(x)
      call write_result(all(x == a), 'bcast_rank1_int64')
    end block
    block
      real(real32) :: x(size(a))
      x = merge(a, 0, is_iop)
      call broadcast(x)
      call write_result(all(x == a), 'bcast_rank1_real32')
    end block
    block
      real(real64) :: x(size(a))
      x = merge(a, 0, is_iop)
      call broadcast(x)
      call write_result(all(x == a), 'bcast_rank1_real64')
    end block
  end subroutine

  ! rank-2 array case
  subroutine bcast_rank2
    integer :: a(2,2)
    a = reshape([4,9,16,25],shape=[2,2])
    block
      integer(int8) :: x(2,2)
      x = merge(a, 0, is_iop)
      call broadcast(x)
      call write_result(all(x == a), 'bcast_rank2_int8')
    end block
    block
      integer(int32) :: x(2,2)
      x = merge(a, 0, is_iop)
      call broadcast(x)
      call write_result(all(x == a), 'bcast_rank2_int32')
    end block
    block
      integer(int64) :: x(2,2)
      x = merge(a, 0, is_iop)
      call broadcast(x)
      call write_result(all(x == a), 'bcast_rank2_int64')
    end block
    block
      real(real32) :: x(2,2)
      x = merge(a, 0, is_iop)
      call broadcast(x)
      call write_result(all(x == a), 'bcast_rank2_real32')
    end block
    block
      real(real64) :: x(2,2)
      x = merge(a, 0, is_iop)
      call broadcast(x)
      call write_result(all(x == a), 'bcast_rank2_real64')
    end block
  end subroutine

  ! rank-3 array case
  subroutine bcast_rank3
    integer :: a(2,1,2)
    a = reshape([4,9,16,25],shape=[2,1,2])
    block
      integer(int8) :: x(2,1,2)
      x = merge(a, 0, is_iop)
      call broadcast(x)
      call write_result(all(x == int(a,int8)), 'bcast_rank3_int8')
    end block
    block
      integer(int32) :: x(2,1,2)
      x = merge(a, 0, is_iop)
      call broadcast(x)
      call write_result(all(x == int(a,int32)), 'bcast_rank3_int32')
    end block
    block
      integer(int64) :: x(2,1,2)
      x = merge(a, 0, is_iop)
      call broadcast(x)
      call write_result(all(x == int(a,int64)), 'bcast_rank3_int64')
    end block
    block
      real(real32) :: x(2,1,2)
      x = merge(a, 0, is_iop)
      call broadcast(x)
      call write_result(all(x == real(a,real32)), 'bcast_rank3_real32')
    end block
    block
      real(real64) :: x(2,1,2)
      x = merge(a, 0, is_iop)
      call broadcast(x)
      call write_result(all(x == real(a,real64)), 'bcast_rank3_real64')
    end block
  end subroutine

  ! rank-2 array section case
  subroutine bcast_section
    integer :: a(3,3), xx(3,3)
    if (is_IOP) then
      a = -1; a(1::2,1::2) = reshape([4,9,16,25],shape=[2,2])
      xx = a
    else
      a = 0; a(1::2,1::2) = reshape([4,9,16,25],shape=[2,2])
      xx = 0
    end if
    block
      integer(int8) :: x(3,3)
      x = xx
      call broadcast(x(1::2,1::2))
      call write_result(all(x == a), 'bcast_section_int8')
    end block
    block
      integer(int32) :: x(3,3)
      x = xx
      call broadcast(x(1::2,1::2))
      call write_result(all(x == a), 'bcast_section_int32')
    end block
    block
      integer(int64) :: x(3,3)
      x = xx
      call broadcast(x(1::2,1::2))
      call write_result(all(x == a), 'bcast_section_int64')
    end block
    block
      real(real32) :: x(3,3)
      x = xx
      call broadcast(x(1::2,1::2))
      call write_result(all(x == a), 'bcast_section_real32')
    end block
    block
      real(real64) :: x(3,3)
      x = xx
      call broadcast(x(1::2,1::2))
      call write_result(all(x == a), 'bcast_section_real64')
    end block
  end subroutine

  ! scalar logical argument
  subroutine bcast_log_scalar
    logical :: a
    a = .false.
    if (is_IOP) a = .true.
    call broadcast(a)
    pass = a
    a = .true.
    if (is_IOP) a = .false.
    call broadcast(a)
    pass = pass .and. .not.a
    call write_result(pass, 'bcast_log_scalar')
  end subroutine

  ! rank-1 array logical argument
  subroutine bcast_log_rank1
    logical :: a(2), b(2)
    b = [.false., .true.]
    if (is_IOP) then
      a = b
    else
      a = .not.b
    end if
    call broadcast(a)
    pass = all(a.eqv.b)
    call write_result(pass, 'bcast_log_rank1')
  end subroutine

  ! rank-2 array logical argument
  subroutine bcast_log_rank2
    logical :: a(2,2), b(2,2)
    b = reshape([.false., .true., .true., .false.], shape=[2,2])
    if (is_IOP) then
      a = b
    else
      a = .not.b
    end if
    call broadcast(a)
    pass = all(a.eqv.b)
    call write_result(pass, 'bcast_log_rank2')
  end subroutine

  ! rank-3 array logical argument
  subroutine bcast_log_rank3
    logical :: a(2,1,2), b(2,1,2)
    b = reshape([.false., .true., .true., .false.], shape=[2,1,2])
    if (is_IOP) then
      a = b
    else
      a = .not.b
    end if
    call broadcast(a)
    pass = all(a.eqv.b)
    call write_result(pass, 'bcast_log_rank3')
  end subroutine

  ! scalar character argument
  subroutine bcast_char_scalar
    character(5) :: a, x
    a = 'fubar'
    x = ''
    if (is_IOP) x = a
    call broadcast(x)
    pass = (x == a)
    call write_result(pass, 'bcast_char_scalar')
  end subroutine

  ! NB: I'm a bit surprised the following character array broadcasts work.
  ! Each element is an odd number of bytes, and I half expected padding to
  ! be included between elements so that the stride between elements to be
  ! some multiple of 2, 4, or 8 bytes. On the MPI C side it is expecting
  ! contiguous bytes. So the storage for the array really is contiguous or
  ! there is some copy-in/copy-out happening at the interface to MPI.

  ! rank-1 character array argument
  subroutine bcast_char_rank1
    character(5) :: a(2), x(2)
    a = ['fubar', 'hello']
    x = ''
    if (is_IOP) x = a
    call broadcast(x)
    pass = all(x == a)
    call write_result(pass, 'bcast_char_rank1')
  end subroutine

  ! rank-2 character array argument
  subroutine bcast_char_rank2
    character(1) :: a(2,2), x(2,2)
    a = reshape(['a','b','c','d'], shape=[2,2])
    x = ''
    if (is_IOP) x = a
    call broadcast(x)
    pass = all(x == a)
    call write_result(pass, 'bcast_char_rank2')
  end subroutine

  ! rank-3 character array argument
  subroutine bcast_char_rank3
    character(1) :: a(2,2,2), x(2,2,2)
    a = reshape(['a','b','c','d', 'e', 'f', 'g', 'h'], shape=[2,2,2])
    x = ''
    if (is_IOP) x = a
    call broadcast(x)
    pass = all(x == a)
    call write_result(pass, 'bcast_char_rank3')
  end subroutine

end program
