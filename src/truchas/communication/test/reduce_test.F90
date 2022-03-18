!! Unit Tests for PARALLEL_COMMUNICATION Reduction Procedures
!!
!! Copyright 2022 Neil N. Carlson <neil.n.carlson@gmail.com>
!! Use subject to the MIT license: https://opensource.org/licenses/MIT
!!

program reduce_test

  use mpi
  use parallel_communication
  use,intrinsic :: iso_fortran_env
  implicit none

  integer :: ierr, status, global_status
  logical :: pass ! local to tests

  call MPI_Init(ierr)
  call init_parallel_communication

  status = 0
  call all_scalar !NB: tests functionality used by write_result!
  call all_rank1
  call any_scalar
  call any_rank1
  call count_scalar
  call count_rank1
  call sum_scalar
  call sum_rank1
  call sum_rank1_mask
  call minval_scalar
  call minval_rank1
  call minval_rank1_mask
  call maxval_scalar
  call maxval_rank1
  call maxval_rank1_mask
  call dot_prod
  call all_rank1_all_zero
  call any_rank1_all_zero
  call count_rank1_all_zero
  call sum_rank1_all_zero
  call minval_rank1_all_zero
  call maxval_rank1_all_zero
  call dot_prod_all_zero
  if (npe > 1) then
    call all_rank1_zero
    call any_rank1_zero
    call count_rank1_zero
    call sum_rank1_zero
    call minval_rank1_zero
    call maxval_rank1_zero
    call dot_prod_zero
  end if

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

  ! scalar argument
  subroutine all_scalar
    logical :: pass
    logical :: a
    a = .true.
    pass = global_all(a)
    if (this_pe == npe) a = .false.
    pass = pass .and. .not.global_all(a)
    call write_result(pass, 'all_scalar')
  end subroutine

  ! vector argument generic case
  subroutine all_rank1
    logical :: pass
    logical, allocatable :: a(:)
    allocate(a(this_pe), source=.true.)
    pass = global_all(a)
    if (this_pe == npe) a(1) = .false.
    pass = pass .and. .not.global_all(a)
    call write_result(pass, 'all_rank1')
  end subroutine

  ! vector argument with 0-sized array
  subroutine all_rank1_zero
    logical :: pass
    logical, allocatable :: a(:)
    allocate(a(this_pe-1), source=.true.)
    pass = global_all(a)
    if (this_pe == npe) a(1) = .false.
    pass = pass .and. .not.global_all(a)
    call write_result(pass, 'all_rank1_zero')
  end subroutine

  ! vector argument generic case
  subroutine all_rank2
    logical :: pass
    logical, allocatable :: a(:,:)
    allocate(a(2,this_pe), source=.true.)
    pass = global_all(a)
    if (this_pe == npe) a(1,1) = .false.
    pass = pass .and. .not.global_all(a)
    call write_result(pass, 'all_rank2')
  end subroutine

  ! vector argument with 0-sized array
  subroutine all_rank2_zero
    logical :: pass
    logical, allocatable :: a(:,:)
    allocate(a(2,this_pe-1), source=.true.)
    pass = global_all(a)
    if (this_pe == npe) a(1,1) = .false.
    pass = pass .and. .not.global_all(a)
    call write_result(pass, 'all_rank2_zero')
  end subroutine

  ! vector argument corner case: all 0-sized arrays
  subroutine all_rank1_all_zero
    logical :: pass
    logical :: a(0)
    pass = global_all(a)
    call write_result(pass, 'all_rank1_all_zero')
  end subroutine

  ! scalar argument
  subroutine any_scalar
    logical :: pass
    logical :: a
    a = .false.
    pass = .not.global_any(a)
    if (this_pe == npe) a = .true.
    pass = pass .and. global_any(a)
    call write_result(pass, 'any_scalar')
  end subroutine

  ! vector argument generic case
  subroutine any_rank1
    logical :: pass
    logical, allocatable :: a(:)
    allocate(a(this_pe), source=.false.)
    pass = .not.global_any(a)
    if (this_pe == npe) a(npe) = .true.
    pass = pass .and. global_any(a)
    call write_result(pass, 'any_rank1')
  end subroutine

  ! vector argument with 0-sized array
  subroutine any_rank1_zero
    logical :: pass
    logical, allocatable :: a(:)
    allocate(a(this_pe-1), source=.false.)
    pass = .not.global_any(a)
    if (this_pe == npe) a(npe-1) = .true.
    pass = pass .and. global_any(a)
    call write_result(pass, 'any_rank1_zero')
  end subroutine

  ! vector argument corner case: all 0-sized arrays
  subroutine any_rank1_all_zero
    logical :: pass
    logical :: a(0)
    pass = .not.global_any(a)
    call write_result(pass, 'any_all_zero')
  end subroutine

  ! scalar argument case
  subroutine count_scalar
    logical :: a
    a = .true.
    if (this_pe == 1) a = .false.
    pass = (global_count(a) == npe-1)
    call write_result(pass, 'count_scalar')
  end subroutine

  ! generic vector argument case
  subroutine count_rank1
    logical, allocatable :: a(:)
    allocate(a(this_pe), source=.true.)
    a(1) = .false.
    pass = (global_count(a) == npe*(npe-1)/2)
    call write_result(pass, 'count_rank1')
  end subroutine

  ! vector arguments with a 0-sized vector
  subroutine count_rank1_zero
    logical, allocatable :: a(:)
    allocate(a(this_pe-1), source=.true.)
    pass = (global_count(a) == npe*(npe-1)/2)
    call write_result(pass, 'count_rank1_zero')
  end subroutine

  ! vector arguments corner case: all 0-sized vectors
  subroutine count_rank1_all_zero
    logical :: a(0)
    pass = (global_count(a) == 0)
    call write_result(pass, 'count_rank1_all_zero')
  end subroutine

  ! scalar argument case
  subroutine sum_scalar
    integer :: a, x
    x = this_pe
    a = npe*(npe+1)/2
    call write_result((global_sum(int(x,int32)) == a), 'sum_scalar_int32')
    call write_result((global_sum(int(x,int64)) == a), 'sum_scalar_int64')
    call write_result((global_sum(real(x,real32)) == a), 'sum_scalar_real32')
    call write_result((global_sum(real(x,real64)) == a), 'sum_scalar_real64')
  end subroutine

  ! generic vector argument case
  subroutine sum_rank1
    integer :: a, x(this_pe)
    x = 1
    a = npe*(npe+1)/2
    call write_result((global_sum(int(x,int32)) == a), 'sum_rank1_int32')
    call write_result((global_sum(int(x,int64)) == a), 'sum_rank1_int64')
    call write_result((global_sum(real(x,real32)) == a), 'sum_rank1_real32')
    call write_result((global_sum(real(x,real64)) == a), 'sum_rank1_real64')
  end subroutine

  ! generic vector argument case with mask
  subroutine sum_rank1_mask
    integer :: x(this_pe)
    logical :: mask(this_pe)
    mask = .false.
    if (is_iop) mask(1) = .true.
    x = 1
    call write_result((global_sum(int(x,int32),mask) == 1), 'sum_rank1_mask_int32')
    call write_result((global_sum(int(x,int64),mask) == 1), 'sum_rank1_mask_int64')
    call write_result((global_sum(real(x,real32),mask) == 1), 'sum_rank1_mask_real32')
    call write_result((global_sum(real(x,real64),mask) == 1), 'sum_rank1_mask_real64')
  end subroutine

  ! vector arguments with a 0-sized vector
  subroutine sum_rank1_zero
    integer :: a, x(this_pe-1)
    x = 1
    a = npe*(npe-1)/2
    call write_result((global_sum(int(x,int32)) == a), 'sum_rank1_zero_int32')
    call write_result((global_sum(int(x,int64)) == a), 'sum_rank1_zero_int64')
    call write_result((global_sum(real(x,real32)) == a), 'sum_rank1_zero_real32')
    call write_result((global_sum(real(x,real64)) == a), 'sum_rank1_zero_real64')
  end subroutine

  ! vector arguments corner case: all 0-sized vectors
  subroutine sum_rank1_all_zero
    call write_result((global_sum([integer(int32)::]) == 0), 'sum_rank1_all_zero_int32')
    call write_result((global_sum([integer(int64)::]) == 0), 'sum_rank1_all_zero_int64')
    call write_result((global_sum([real(real32)::]) == 0), 'sum_rank1_all_zero_real32')
    call write_result((global_sum([real(real64)::]) == 0), 'sum_rank1_all_zero_real64')
  end subroutine

  ! scalar argument case
  subroutine minval_scalar
    integer :: a, x
    x = -this_pe
    a = -npe
    call write_result((global_minval(int(x,int32)) == a), 'minval_scalar_int32')
    call write_result((global_minval(int(x,int64)) == a), 'minval_scalar_int64')
    call write_result((global_minval(real(x,real32)) == a), 'minval_scalar_real32')
    call write_result((global_minval(real(x,real64)) == a), 'minval_scalar_real64')
  end subroutine

  ! generic vector argument case
  subroutine minval_rank1
    integer :: a, x(this_pe), j
    x = [(-j, j=1,size(x))]
    a = -npe
    call write_result((global_minval(int(x,int32)) == a), 'minval_rank1_int32')
    call write_result((global_minval(int(x,int64)) == a), 'minval_rank1_int64')
    call write_result((global_minval(real(x,real32)) == a), 'minval_rank1_real32')
    call write_result((global_minval(real(x,real64)) == a), 'minval_rank1_real64')
  end subroutine

  ! generic vector argument case with mask
  subroutine minval_rank1_mask
    integer :: x(this_pe), j
    logical :: mask(this_pe)
    mask = .false.
    if (is_iop) mask(1) = .true.
    x = [(-j, j=1,size(x))]
    call write_result((global_minval(int(x,int32),mask) == -1), 'minval_rank1_mask_int32')
    call write_result((global_minval(int(x,int64),mask) == -1), 'minval_rank1_mask_int64')
    call write_result((global_minval(real(x,real32),mask) == -1), 'minval_rank1_mask_real32')
    call write_result((global_minval(real(x,real64),mask) == -1), 'minval_rank1_mask_real64')
  end subroutine

  ! vector arguments with 0-sized vector
  subroutine minval_rank1_zero
    integer :: a, x(this_pe-1), j
    x = [(-j, j=1,size(x))]
    a = -npe+1
    call write_result((global_minval(int(x,int32)) == a), 'minval_rank1_zero_int32')
    call write_result((global_minval(int(x,int64)) == a), 'minval_rank1_zero_int64')
    call write_result((global_minval(real(x,real32)) == a), 'minval_rank1_zero_real32')
    call write_result((global_minval(real(x,real64)) == a), 'minval_rank1_zero_real64')
  end subroutine

  ! vector arguments corner case: all 0-sized vectors
  subroutine minval_rank1_all_zero
    call write_result((global_minval([integer(int32)::]) == minval([integer(int32)::])), 'minval_rank1_all_zero_int32')
    call write_result((global_minval([integer(int64)::]) == minval([integer(int64)::])), 'minval_rank1_all_zero_int64')
    call write_result((global_minval([real(real32)::])   == minval([real(real32)::]  )), 'minval_rank1_all_zero_real32')
    call write_result((global_minval([real(real64)::])   == minval([real(real64)::]  )), 'minval_rank1_all_zero_real64')
  end subroutine

  ! scalar argument case
  subroutine maxval_scalar
    integer :: a, x
    x = this_pe
    a = npe
    call write_result((global_maxval(int(x,int32)) == a), 'maxval_scalar_int32')
    call write_result((global_maxval(int(x,int64)) == a), 'maxval_scalar_int64')
    call write_result((global_maxval(real(x,real32)) == a), 'maxval_scalar_real32')
    call write_result((global_maxval(real(x,real64)) == a), 'maxval_scalar_real64')
  end subroutine

  ! generic vector argument case
  subroutine maxval_rank1
    integer :: a, x(this_pe), j
    x = [(j, j=1,size(x))]
    a = npe
    call write_result((global_maxval(int(x,int32)) == a), 'maxval_rank1_int32')
    call write_result((global_maxval(int(x,int64)) == a), 'maxval_rank1_int64')
    call write_result((global_maxval(real(x,real32)) == a), 'maxval_rank1_real32')
    call write_result((global_maxval(real(x,real64)) == a), 'maxval_rank1_real64')
  end subroutine

  ! generic vector argument case with mask
  subroutine maxval_rank1_mask
    integer :: x(this_pe), j
    logical :: mask(this_pe)
    mask = .false.
    if (is_iop) mask(1) = .true.
    x = [(j, j=1,size(x))]
    call write_result((global_maxval(int(x,int32),mask) == 1), 'maxval_rank1_mask_int32')
    call write_result((global_maxval(int(x,int64),mask) == 1), 'maxval_rank1_mask_int64')
    call write_result((global_maxval(real(x,real32),mask) == 1), 'maxval_rank1_mask_real32')
    call write_result((global_maxval(real(x,real64),mask) == 1), 'maxval_rank1_mask_real64')
  end subroutine

  ! vector arguments with 0-sized vector
  subroutine maxval_rank1_zero
    integer :: a, x(this_pe-1), j
    x = [(j, j=1,size(x))]
    a = npe-1
    call write_result((global_maxval(int(x,int32)) == a), 'maxval_rank1_zero_int32')
    call write_result((global_maxval(int(x,int64)) == a), 'maxval_rank1_zero_int64')
    call write_result((global_maxval(real(x,real32)) == a), 'maxval_rank1_zero_real32')
    call write_result((global_maxval(real(x,real64)) == a), 'maxval_rank1_zero_real64')
  end subroutine

  ! vector arguments corner case: all 0-sized vectors
  subroutine maxval_rank1_all_zero
    call write_result((global_maxval([integer(int32)::]) == maxval([integer(int32)::])), 'maxval_rank1_all_zero_int32')
    call write_result((global_maxval([integer(int64)::]) == maxval([integer(int64)::])), 'maxval_rank1_all_zero_int64')
    call write_result((global_maxval([real(real32)::])   == maxval([real(real32)::]  )), 'maxval_rank1_all_zero_real32')
    call write_result((global_maxval([real(real64)::])   == maxval([real(real64)::]  )), 'maxval_rank1_all_zero_real64')
  end subroutine

  ! generic case
  subroutine dot_prod
    integer :: a, x(this_pe), y(this_pe), j
    a = (this_pe*(this_pe-1))/2
    x = [(a+j, j=1,size(x))]
    y = 1
    a = (npe*(npe+1))/2
    a = (a*(a+1))/2
    call write_result((global_dot_product(real(x,real32),real(y,real32)) == a), 'dot_prod_real32')
    call write_result((global_dot_product(real(x,real64),real(y,real64)) == a), 'dot_prod_real64')
  end subroutine

  ! case with 0-sized vector
  subroutine dot_prod_zero
    integer :: a, x(this_pe-1), y(this_pe-1), j
    a = ((this_pe-1)*(this_pe-2))/2
    x = [(a+j, j=1,size(x))]
    y = 1
    a = (npe*(npe-1))/2
    a = (a*(a+1))/2
    call write_result((global_dot_product(real(x,real32),real(y,real32)) == a), 'dot_prod_zero_real32')
    call write_result((global_dot_product(real(x,real64),real(y,real64)) == a), 'dot_prod_zero_real64')
  end subroutine

  ! corner case with all 0-sized vectors
  subroutine dot_prod_all_zero
    call write_result((global_dot_product([real(real32)::],[real(real32)::]) == 0), 'dot_prod_all_zero_real32')
    call write_result((global_dot_product([real(real64)::],[real(real64)::]) == 0), 'dot_prod_all_zero_real64')
  end subroutine

end program
