!! Implementation of PARALLEL_COMMUNICATION BROADCAST Procedures
!!
!! Copyright 2022 Neil N. Carlson <neil.n.carlson@gmail.com>
!! Use subject to the MIT license: https://opensource.org/licenses/MIT
!!

#include "f90_assert.fpp"

submodule(parallel_communication) broadcast_impl
implicit none
contains

!!!! BROADCAST SCALAR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine bcast_i1_0(scalar)
    integer(i1), intent(inout) :: scalar
    integer :: ierr
    call MPI_Bcast(scalar, 1, MPI_INTEGER1, root, comm, ierr)
  end subroutine

  module subroutine bcast_i4_0(scalar)
    integer(i4), intent(inout) :: scalar
    integer :: ierr
    call MPI_Bcast(scalar, 1, MPI_INTEGER4, root, comm, ierr)
  end subroutine

  module subroutine bcast_i8_0(scalar)
    integer(i8), intent(inout) :: scalar
    integer :: ierr
    call MPI_Bcast(scalar, 1, MPI_INTEGER8, root, comm, ierr)
  end subroutine

  module subroutine bcast_r4_0(scalar)
    real(r4), intent(inout) :: scalar
    integer :: ierr
    call MPI_Bcast(scalar, 1, MPI_REAL4, root, comm, ierr)
  end subroutine

  module subroutine bcast_r8_0(scalar)
    real(r8), intent(inout) :: scalar
    integer :: ierr
    call MPI_Bcast(scalar, 1, MPI_REAL8, root, comm, ierr)
  end subroutine

  module subroutine bcast_c4_0(scalar)
    complex(r4), intent(inout) :: scalar
    integer :: ierr
    call MPI_Bcast(scalar, 1, MPI_COMPLEX8, root, comm, ierr)
  end subroutine

  module subroutine bcast_c8_0(scalar)
    complex(r8), intent(inout) :: scalar
    integer :: ierr
    call MPI_Bcast(scalar, 1, MPI_COMPLEX16, root, comm, ierr)
  end subroutine

  module subroutine bcast_dl_0(scalar)
    logical, intent(inout) :: scalar
    integer :: ierr
    call MPI_Bcast(scalar, 1, MPI_LOGICAL, root, comm, ierr)
  end subroutine

  module subroutine bcast_char_0(scalar)
    character(*), intent(inout) :: scalar
    integer :: ierr
    call MPI_Bcast(scalar, len(scalar), MPI_CHARACTER, root, comm, ierr)
  end subroutine

!!!! BROADCAST RANK-1 ARRAY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine bcast_i1_1(vector)
    integer(i1), intent(inout) :: vector(:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_INTEGER1, root, comm, ierr)
  end subroutine

  module subroutine bcast_i4_1(vector)
    integer(i4), intent(inout) :: vector(:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_INTEGER4, root, comm, ierr)
  end subroutine

  module subroutine bcast_i8_1(vector)
    integer(i8), intent(inout) :: vector(:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_INTEGER8, root, comm, ierr)
  end subroutine

  module subroutine bcast_r4_1(vector)
    real(r4), intent(inout) :: vector(:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_REAL4, root, comm, ierr)
  end subroutine

  module subroutine bcast_r8_1(vector)
    real(r8), intent(inout) :: vector(:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_REAL8, root, comm, ierr)
  end subroutine

  module subroutine bcast_c4_1(vector)
    complex(r4), intent(inout) :: vector(:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_COMPLEX8, root, comm, ierr)
  end subroutine

  module subroutine bcast_c8_1(vector)
    complex(r8), intent(inout) :: vector(:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_COMPLEX16, root, comm, ierr)
  end subroutine

  module subroutine bcast_dl_1(vector)
    logical, intent(inout) :: vector(:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_LOGICAL, root, comm, ierr)
  end subroutine

  module subroutine bcast_char_1(vector)
    character(*), intent(inout) :: vector(:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector)*len(vector), MPI_CHARACTER, root, comm, ierr)
  end subroutine

!!!! BROADCAST RANK-2 ARRAY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine bcast_i1_2(vector)
    integer(i1), intent(inout) :: vector(:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_INTEGER1, root, comm, ierr)
  end subroutine

  module subroutine bcast_i4_2(vector)
    integer(i4), intent(inout) :: vector(:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_INTEGER4, root, comm, ierr)
  end subroutine

  module subroutine bcast_i8_2(vector)
    integer(i8), intent(inout) :: vector(:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_INTEGER8, root, comm, ierr)
  end subroutine

  module subroutine bcast_r4_2(vector)
    real(r4), intent(inout) :: vector(:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_REAL4, root, comm, ierr)
  end subroutine

  module subroutine bcast_r8_2(vector)
    real(r8), intent(inout) :: vector(:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_REAL8, root, comm, ierr)
  end subroutine

  module subroutine bcast_c4_2(vector)
    complex(r4), intent(inout) :: vector(:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_COMPLEX8, root, comm, ierr)
  end subroutine

  module subroutine bcast_c8_2(vector)
    complex(r8), intent(inout) :: vector(:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_COMPLEX16, root, comm, ierr)
  end subroutine

  module subroutine bcast_dl_2(vector)
    logical, intent(inout) :: vector(:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_LOGICAL, root, comm, ierr)
  end subroutine

  module subroutine bcast_char_2(vector)
    character(*), intent(inout) :: vector(:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector)*len(vector), MPI_CHARACTER, root, comm, ierr)
  end subroutine

!!!! BROADCAST RANK-3 ARRAY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine bcast_i1_3(vector)
    integer(i1), intent(inout) :: vector(:,:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_INTEGER1, root, comm, ierr)
  end subroutine

  module subroutine bcast_i4_3(vector)
    integer(i4), intent(inout) :: vector(:,:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_INTEGER4, root, comm, ierr)
  end subroutine

  module subroutine bcast_i8_3(vector)
    integer(i8), intent(inout) :: vector(:,:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_INTEGER8, root, comm, ierr)
  end subroutine

  module subroutine bcast_r4_3(vector)
    real(r4), intent(inout) :: vector(:,:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_REAL4, root, comm, ierr)
  end subroutine

  module subroutine bcast_r8_3(vector)
    real(r8), intent(inout) :: vector(:,:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_REAL8, root, comm, ierr)
  end subroutine

  module subroutine bcast_c4_3(vector)
    complex(r4), intent(inout) :: vector(:,:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_COMPLEX8, root, comm, ierr)
  end subroutine

  module subroutine bcast_c8_3(vector)
    complex(r8), intent(inout) :: vector(:,:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_COMPLEX16, root, comm, ierr)
  end subroutine

  module subroutine bcast_dl_3(vector)
    logical, intent(inout) :: vector(:,:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_LOGICAL, root, comm, ierr)
  end subroutine

  module subroutine bcast_char_3(vector)
    character(*), intent(inout) :: vector(:,:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector)*len(vector), MPI_CHARACTER, root, comm, ierr)
  end subroutine

end submodule broadcast_impl
