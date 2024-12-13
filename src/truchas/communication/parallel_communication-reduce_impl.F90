!! Implementation of PARALLEL_COMMUNICATION Reduction Procedures
!!
!! Copyright 2022 Neil N. Carlson <neil.n.carlson@gmail.com>
!! Use subject to the MIT license: https://opensource.org/licenses/MIT
!!

#include "f90_assert.fpp"

submodule(parallel_communication) reduce_impl
implicit none
contains

!!!! GLOBAL ANY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module function any_0(mask) result(global)
    logical, intent(in) :: mask
    logical :: global
    integer :: ierr
    call MPI_Allreduce(mask, global, 1, MPI_LOGICAL, MPI_LOR, comm, ierr)
  end function

  module function any_1(mask) result(global)
    logical, intent(in) :: mask(:)
    logical :: global, local
    integer :: ierr
    local = any(mask)
    call MPI_Allreduce(local, global, 1, MPI_LOGICAL, MPI_LOR, comm, ierr)
  end function

  module function any_2(mask) result(global)
    logical, intent(in) :: mask(:,:)
    logical :: global, local
    integer :: ierr
    local = any(mask)
    call MPI_Allreduce(local, global, 1, MPI_LOGICAL, MPI_LOR, comm, ierr)
  end function

!!!! GLOBAL ALL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module function all_0(mask) result(global)
    logical, intent(in) :: mask
    logical :: global
    integer :: ierr
    call MPI_Allreduce(mask, global, 1, MPI_LOGICAL, MPI_LAND, comm, ierr)
  end function

  module function all_1(mask) result(global)
    logical, intent(in) :: mask(:)
    logical :: global, local
    integer :: ierr
    local = all(mask)
    call MPI_Allreduce(local, global, 1, MPI_LOGICAL, MPI_LAND, comm, ierr)
  end function

  module function all_2(mask) result(global)
    logical, intent(in) :: mask(:,:)
    logical :: global, local
    integer :: ierr
    local = all(mask)
    call MPI_Allreduce(local, global, 1, MPI_LOGICAL, MPI_LAND, comm, ierr)
  end function

!!!! GLOBAL COUNT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module function count_0(mask) result(global)
    logical, intent(in) :: mask
    integer :: global, local
    integer :: ierr
    local = merge(1, 0, mask)
    call MPI_Allreduce(local, global, 1, MPI_INTEGER, MPI_SUM, comm, ierr)
  end function

  module function count_1(mask) result(global)
    logical, intent(in) :: mask(:)
    integer :: global, local
    integer :: ierr
    local = count(mask)
    call MPI_Allreduce(local, global, 1, MPI_INTEGER, MPI_SUM, comm, ierr)
  end function

  module function count_2(mask) result(global)
    logical, intent(in) :: mask(:,:)
    integer :: global, local
    integer :: ierr
    local = count(mask)
    call MPI_Allreduce(local, global, 1, MPI_INTEGER, MPI_SUM, comm, ierr)
  end function

!!!! GLOBAL SUM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module function sum_i4_0(a) result(s)
    integer(i4), intent(in) :: a
    integer(i4) :: s
    integer :: ierr
    call MPI_Allreduce(a, s, 1, MPI_INTEGER4, MPI_SUM, comm, ierr)
  end function

  module function sum_i8_0(a) result(s)
    integer(i8), intent(in) :: a
    integer(i8) :: s
    integer :: ierr
    call MPI_Allreduce(a, s, 1, MPI_INTEGER8, MPI_SUM, comm, ierr)
  end function

  module function sum_r4_0(a) result(s)
    real(r4), intent(in) :: a
    real(r4) :: s
    integer :: ierr
    call MPI_Allreduce(a, s, 1, MPI_REAL4, MPI_SUM, comm, ierr)
  end function

  module function sum_r8_0(a) result(s)
    real(r8), intent(in) :: a
    real(r8) :: s
    integer :: ierr
    call MPI_Allreduce(a, s, 1, MPI_REAL8, MPI_SUM, comm, ierr)
  end function

  module function sum_c4_0(a) result(s)
    complex(r4), intent(in) :: a
    complex(r4) :: s
    integer :: ierr
    call MPI_Allreduce(a, s, 1, MPI_COMPLEX8, MPI_SUM, comm, ierr)
  end function

  module function sum_c8_0(a) result(s)
    complex(r8), intent(in) :: a
    complex(r8) :: s
    integer :: ierr
    call MPI_Allreduce(a, s, 1, MPI_COMPLEX16, MPI_SUM, comm, ierr)
  end function

  module function sum_i4_1(a, mask) result(s)
    integer(i4), intent(in) :: a(:)
    logical, intent(in), optional :: mask(:)
    integer(i4) :: s
    integer :: ierr
    call MPI_Allreduce(sum(a,mask), s, 1, MPI_INTEGER4, MPI_SUM, comm, ierr)
  end function

  module function sum_i8_1(a, mask) result(s)
    integer(i8), intent(in) :: a(:)
    logical, intent(in), optional :: mask(:)
    integer(i8) :: s
    integer :: ierr
    call MPI_Allreduce(sum(a,mask), s, 1, MPI_INTEGER8, MPI_SUM, comm, ierr)
  end function

  module function sum_r4_1(a, mask) result(s)
    real(r4), intent(in) :: a(:)
    logical, intent(in), optional :: mask(:)
    real(r4) :: s
    integer :: ierr
    call MPI_Allreduce(sum(a,mask), s, 1, MPI_REAL4, MPI_SUM, comm, ierr)
  end function

  module function sum_r8_1(a, mask) result(s)
    real(r8), intent(in) :: a(:)
    logical, intent(in), optional :: mask(:)
    real(r8) :: s
    integer :: ierr
    call MPI_Allreduce(sum(a,mask), s, 1, MPI_REAL8, MPI_SUM, comm, ierr)
  end function

  module function sum_c4_1(a, mask) result(s)
    complex(r4), intent(in) :: a(:)
    logical, intent(in), optional :: mask(:)
    complex(r4) :: s
    integer :: ierr
    call MPI_Allreduce(sum(a,mask), s, 1, MPI_COMPLEX8, MPI_SUM, comm, ierr)
  end function

  module function sum_c8_1(a, mask) result(s)
    complex(r8), intent(in) :: a(:)
    logical, intent(in), optional :: mask(:)
    complex(r8) :: s
    integer :: ierr
    call MPI_Allreduce(sum(a,mask), s, 1, MPI_COMPLEX16, MPI_SUM, comm, ierr)
  end function

!!!! GLOBAL MINVAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module function minval_i4_0(a) result(s)
    integer(i4), intent(in) :: a
    integer(i4) :: s
    integer :: ierr
    call MPI_Allreduce(a, s, 1, MPI_INTEGER4, MPI_MIN, comm, ierr)
  end function

  module function minval_i8_0(a) result(s)
    integer(i8), intent(in) :: a
    integer(i8) :: s
    integer :: ierr
    call MPI_Allreduce(a, s, 1, MPI_INTEGER8, MPI_MIN, comm, ierr)
  end function

  module function minval_r4_0(a) result(s)
    real(r4), intent(in) :: a
    real(r4) :: s
    integer :: ierr
    call MPI_Allreduce(a, s, 1, MPI_REAL4, MPI_MIN, comm, ierr)
  end function

  module function minval_r8_0(a) result(s)
    real(r8), intent(in) :: a
    real(r8) :: s
    integer :: ierr
    call MPI_Allreduce(a, s, 1, MPI_REAL8, MPI_MIN, comm, ierr)
  end function

  module function minval_i4_1(a, mask) result(s)
    integer(i4), intent(in) :: a(:)
    logical, intent(in), optional :: mask(:)
    integer(i4) :: s
    integer :: ierr
    call MPI_Allreduce(minval(a,mask), s, 1, MPI_INTEGER4, MPI_MIN, comm, ierr)
  end function

  module function minval_i8_1(a, mask) result(s)
    integer(i8), intent(in) :: a(:)
    logical, intent(in), optional :: mask(:)
    integer(i8) :: s
    integer :: ierr
    call MPI_Allreduce(minval(a,mask), s, 1, MPI_INTEGER8, MPI_MIN, comm, ierr)
  end function

  module function minval_r4_1(a, mask) result(s)
    real(r4), intent(in) :: a(:)
    logical, intent(in), optional :: mask(:)
    real(r4) :: s
    integer :: ierr
    call MPI_Allreduce(minval(a,mask), s, 1, MPI_REAL4, MPI_MIN, comm, ierr)
  end function

  module function minval_r8_1(a, mask) result(s)
    real(r8), intent(in) :: a(:)
    logical, intent(in), optional :: mask(:)
    real(r8) :: s
    integer :: ierr
    call MPI_Allreduce(minval(a,mask), s, 1, MPI_REAL8, MPI_MIN, comm, ierr)
  end function

!!!! GLOBAL MAXVAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module function maxval_i4_0(a) result(s)
    integer(i4), intent(in) :: a
    integer(i4) :: s
    integer :: ierr
    call MPI_Allreduce(a, s, 1, MPI_INTEGER4, MPI_MAX, comm, ierr)
  end function

  module function maxval_i8_0(a) result(s)
    integer(i8), intent(in) :: a
    integer(i8) :: s
    integer :: ierr
    call MPI_Allreduce(a, s, 1, MPI_INTEGER8, MPI_MAX, comm, ierr)
  end function

  module function maxval_r4_0(a) result(s)
    real(r4), intent(in) :: a
    real(r4) :: s
    integer :: ierr
    call MPI_Allreduce(a, s, 1, MPI_REAL4, MPI_MAX, comm, ierr)
  end function

  module function maxval_r8_0(a) result(s)
    real(r8), intent(in) :: a
    real(r8) :: s
    integer :: ierr
    call MPI_Allreduce(a, s, 1, MPI_REAL8, MPI_MAX, comm, ierr)
  end function

  module function maxval_i4_1(a, mask) result(s)
    integer(i4), intent(in) :: a(:)
    logical, intent(in), optional :: mask(:)
    integer(i4) :: s
    integer :: ierr
    call MPI_Allreduce(maxval(a,mask), s, 1, MPI_INTEGER4, MPI_MAX, comm, ierr)
  end function

  module function maxval_i8_1(a, mask) result(s)
    integer(i8), intent(in) :: a(:)
    logical, intent(in), optional :: mask(:)
    integer(i8) :: s
    integer :: ierr
    call MPI_Allreduce(maxval(a,mask), s, 1, MPI_INTEGER8, MPI_MAX, comm, ierr)
  end function

  module function maxval_r4_1(a, mask) result(s)
    real(r4), intent(in) :: a(:)
    logical, intent(in), optional :: mask(:)
    real(r4) :: s
    integer :: ierr
    call MPI_Allreduce(maxval(a,mask), s, 1, MPI_REAL4, MPI_MAX, comm, ierr)
  end function

  module function maxval_r8_1(a, mask) result(s)
    real(r8), intent(in) :: a(:)
    logical, intent(in), optional :: mask(:)
    real(r8) :: s
    integer :: ierr
    call MPI_Allreduce(maxval(a,mask), s, 1, MPI_REAL8, MPI_MAX, comm, ierr)
  end function

!!!! GLOBAL DOT PRODUCT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module function dot_prod_r4(a, b) result(dp)
    real(r4), intent(in) :: a(:), b(:)
    real(r4) :: dp
    integer :: ierr
    call MPI_Allreduce(dot_product(a,b), dp, 1, MPI_REAL4, MPI_SUM, comm, ierr)
  end function

  module function dot_prod_r8(a, b) result(dp)
    real(r8), intent(in) :: a(:), b(:)
    real(r8) :: dp
    integer :: ierr
    call MPI_Allreduce(dot_product(a,b), dp, 1, MPI_REAL8, MPI_SUM, comm, ierr)
  end function

  module function dot_prod_c4(a, b) result(dp)
    complex(r4), intent(in) :: a(:), b(:)
    complex(r4) :: dp
    integer :: ierr
    call MPI_Allreduce(dot_product(a,b), dp, 1, MPI_COMPLEX8, MPI_SUM, comm, ierr)
  end function

  module function dot_prod_c8(a, b) result(dp)
    complex(r8), intent(in) :: a(:), b(:)
    complex(r8) :: dp
    integer :: ierr
    call MPI_Allreduce(dot_product(a,b), dp, 1, MPI_COMPLEX16, MPI_SUM, comm, ierr)
  end function

!!!! GLOBAL NORM2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module function norm2_r4(a) result(nrm2)
    real(r4), intent(in) :: a(:)
    real(r4) :: nrm2
    integer :: ierr
    call MPI_Allreduce(norm2(a)**2, nrm2, 1, MPI_REAL4, MPI_SUM, comm, ierr)
    nrm2 = sqrt(nrm2)
  end function

  module function norm2_r8(a) result(nrm2)
    real(r8), intent(in) :: a(:)
    real(r8) :: nrm2
    integer :: ierr
    call MPI_Allreduce(norm2(a)**2, nrm2, 1, MPI_REAL8, MPI_SUM, comm, ierr)
    nrm2 = sqrt(nrm2)
  end function

  module function norm2_c4(a) result(nrm2)
    complex(r4), intent(in) :: a(:)
    real(r4) :: nrm2, tmp
    integer :: ierr
    tmp = norm2(a%re)**2 + norm2(a%im)**2
    call MPI_Allreduce(tmp, nrm2, 1, MPI_REAL4, MPI_SUM, comm, ierr)
    nrm2 = sqrt(nrm2)
  end function

  module function norm2_c8(a) result(nrm2)
    complex(r8), intent(in) :: a(:)
    real(r8) :: nrm2, tmp
    integer :: ierr
    tmp = norm2(a%re)**2 + norm2(a%im)**2
    call MPI_Allreduce(tmp, nrm2, 1, MPI_REAL8, MPI_SUM, comm, ierr)
    nrm2 = sqrt(nrm2)
  end function

end submodule reduce_impl
