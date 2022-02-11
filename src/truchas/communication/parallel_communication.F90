!! PARALLEL_COMMUNICATION
!!
!! This is a reimplementation of the original module that replaces internal
!! use of PGSLib with direct use of MPI. All interfaces are unchanged.
!!
!! * Retains numbering processes (PEs) starting with 1.
!! * The IO process is PE 1 (MPI rank 0). Was user-choice with PGSLib.
!! * User may initialize MPI; if not, module will initialize. Whoever
!!   initializes MPI is responsible for finalizing MPI.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright 2022 Neil N. Carlson <neil.n.carlson@gmail.com>
!! Use subject to the MIT license: https://opensource.org/licenses/MIT
!!

#include "f90_assert.fpp"

module parallel_communication

  use mpi
  use,intrinsic :: iso_fortran_env, only: int8, int32, int64, real32, real64
  implicit none
  private

  public :: init_parallel_communication, halt_parallel_communication, abort_parallel_communication
  public :: broadcast, distribute, collate
  public :: global_any, global_all, global_count
  public :: global_sum, global_minval, global_maxval, global_dot_product
  public :: global_minloc, global_maxloc, global_maxloc_sub

  integer, parameter :: root = 0
  integer, parameter, public :: comm = MPI_COMM_WORLD
  integer, parameter, public :: io_pe = root + 1

  integer, public, protected :: npe = 1
  integer, public, protected :: this_pe = io_pe
  logical, public, protected :: is_iop = .true.

  interface broadcast
    module subroutine bcast_i1_0(scalar)
      integer(int8), intent(inout) :: scalar
    end subroutine
    module subroutine bcast_i4_0(scalar)
      integer(int32), intent(inout) :: scalar
    end subroutine
    module subroutine bcast_i8_0(scalar)
      integer(int64), intent(inout) :: scalar
    end subroutine
    module subroutine bcast_r4_0(scalar)
      real(real32), intent(inout) :: scalar
    end subroutine
    module subroutine bcast_r8_0(scalar)
      real(real64), intent(inout) :: scalar
    end subroutine
    module subroutine bcast_log_0(scalar)
      logical, intent(inout) :: scalar
    end subroutine
    module subroutine bcast_char_0(scalar)
      character(*), intent(inout) :: scalar
    end subroutine
    module subroutine bcast_i1_1(vector)
      integer(int8), intent(inout) :: vector(:)
    end subroutine
    module subroutine bcast_i4_1(vector)
      integer(int32), intent(inout) :: vector(:)
    end subroutine
    module subroutine bcast_i8_1(vector)
      integer(int64), intent(inout) :: vector(:)
    end subroutine
    module subroutine bcast_r4_1(vector)
      real(real32), intent(inout) :: vector(:)
    end subroutine
    module subroutine bcast_r8_1(vector)
      real(real64), intent(inout) :: vector(:)
    end subroutine
    module subroutine bcast_log_1(vector)
      logical, intent(inout) :: vector(:)
    end subroutine
    module subroutine bcast_char_1(vector)
      character(*), intent(inout) :: vector(:)
    end subroutine
    module subroutine bcast_i1_2(vector)
      integer(int8), intent(inout) :: vector(:,:)
    end subroutine
    module subroutine bcast_i4_2(vector)
      integer(int32), intent(inout) :: vector(:,:)
    end subroutine
    module subroutine bcast_i8_2(vector)
      integer(int64), intent(inout) :: vector(:,:)
    end subroutine
    module subroutine bcast_r4_2(vector)
      real(real32), intent(inout) :: vector(:,:)
    end subroutine
    module subroutine bcast_r8_2(vector)
      real(real64), intent(inout) :: vector(:,:)
    end subroutine
    module subroutine bcast_log_2(vector)
      logical, intent(inout) :: vector(:,:)
    end subroutine
    module subroutine bcast_char_2(vector)
      character(*), intent(inout) :: vector(:,:)
    end subroutine
    module subroutine bcast_i1_3(vector)
      integer(int8), intent(inout) :: vector(:,:,:)
    end subroutine
    module subroutine bcast_i4_3(vector)
      integer(int32), intent(inout) :: vector(:,:,:)
    end subroutine
    module subroutine bcast_i8_3(vector)
      integer(int64), intent(inout) :: vector(:,:,:)
    end subroutine
    module subroutine bcast_r4_3(vector)
      real(real32), intent(inout) :: vector(:,:,:)
    end subroutine
    module subroutine bcast_r8_3(vector)
      real(real64), intent(inout) :: vector(:,:,:)
    end subroutine
    module subroutine bcast_log_3(vector)
      logical, intent(inout) :: vector(:,:,:)
    end subroutine
    module subroutine bcast_char_3(vector)
      character(*), intent(inout) :: vector(:,:,:)
    end subroutine
  end interface

  interface distribute
    module subroutine dist_i1_0(scalar_out, scalarv_in)
      integer(int8), intent(out) :: scalar_out
      integer(int8), intent(in)  :: scalarv_in(:)
    end subroutine
    module subroutine dist_i4_0(scalar_out, scalarv_in)
      integer(int32), intent(out) :: scalar_out
      integer(int32), intent(in)  :: scalarv_in(:)
    end subroutine
    module subroutine dist_i8_0(scalar_out, scalarv_in)
      integer(int64), intent(out) :: scalar_out
      integer(int64), intent(in)  :: scalarv_in(:)
    end subroutine
    module subroutine dist_r4_0(scalar_out, scalarv_in)
      real(real32), intent(out) :: scalar_out
      real(real32), intent(in)  :: scalarv_in(:)
    end subroutine
    module subroutine dist_r8_0(scalar_out, scalarv_in)
      real(real64), intent(out) :: scalar_out
      real(real64), intent(in)  :: scalarv_in(:)
    end subroutine
    module subroutine dist_log_0(scalar_out, scalarv_in)
      logical, intent(out) :: scalar_out
      logical, intent(in)  :: scalarv_in(:)
    end subroutine
    module subroutine dist_i1_1(vector_out, vector_in)
      integer(int8), intent(out) :: vector_out(:)
      integer(int8), intent(in)  :: vector_in(:)
    end subroutine
    module subroutine dist_i4_1(vector_out, vector_in)
      integer(int32), intent(out) :: vector_out(:)
      integer(int32), intent(in)  :: vector_in(:)
    end subroutine
    module subroutine dist_i8_1(vector_out, vector_in)
      integer(int64), intent(out) :: vector_out(:)
      integer(int64), intent(in)  :: vector_in(:)
    end subroutine
    module subroutine dist_r4_1(vector_out, vector_in)
      real(real32), intent(out) :: vector_out(:)
      real(real32), intent(in)  :: vector_in(:)
    end subroutine
    module subroutine dist_r8_1(vector_out, vector_in)
      real(real64), intent(out) :: vector_out(:)
      real(real64), intent(in)  :: vector_in(:)
    end subroutine
    module subroutine dist_log_1(vector_out, vector_in)
      logical, intent(out) :: vector_out(:)
      logical, intent(in)  :: vector_in(:)
    end subroutine
    module subroutine dist_i1_2(vector_out, vector_in)
      integer(int8), intent(out) :: vector_out(:,:)
      integer(int8), intent(in)  :: vector_in(:,:)
    end subroutine
    module subroutine dist_i4_2(vector_out, vector_in)
      integer(int32), intent(out) :: vector_out(:,:)
      integer(int32), intent(in)  :: vector_in(:,:)
    end subroutine
    module subroutine dist_i8_2(vector_out, vector_in)
      integer(int64), intent(out) :: vector_out(:,:)
      integer(int64), intent(in)  :: vector_in(:,:)
    end subroutine
    module subroutine dist_r4_2(vector_out, vector_in)
      real(real32), intent(out) :: vector_out(:,:)
      real(real32), intent(in)  :: vector_in(:,:)
    end subroutine
    module subroutine dist_r8_2(vector_out, vector_in)
      real(real64), intent(out) :: vector_out(:,:)
      real(real64), intent(in)  :: vector_in(:,:)
    end subroutine
    module subroutine dist_log_2(vector_out, vector_in)
      logical, intent(out) :: vector_out(:,:)
      logical, intent(in)  :: vector_in(:,:)
    end subroutine
    module subroutine dist_i1_3(vector_out, vector_in)
      integer(int8), intent(out) :: vector_out(:,:,:)
      integer(int8), intent(in)  :: vector_in(:,:,:)
    end subroutine
    module subroutine dist_i4_3(vector_out, vector_in)
      integer(int32), intent(out) :: vector_out(:,:,:)
      integer(int32), intent(in)  :: vector_in(:,:,:)
    end subroutine
    module subroutine dist_i8_3(vector_out, vector_in)
      integer(int64), intent(out) :: vector_out(:,:,:)
      integer(int64), intent(in)  :: vector_in(:,:,:)
    end subroutine
    module subroutine dist_r4_3(vector_out, vector_in)
      real(real32), intent(out) :: vector_out(:,:,:)
      real(real32), intent(in)  :: vector_in(:,:,:)
    end subroutine
    module subroutine dist_r8_3(vector_out, vector_in)
      real(real64), intent(out) :: vector_out(:,:,:)
      real(real64), intent(in)  :: vector_in(:,:,:)
    end subroutine
    module subroutine dist_log_3(vector_out, vector_in)
      logical, intent(out) :: vector_out(:,:,:)
      logical, intent(in)  :: vector_in(:,:,:)
    end subroutine
  end interface

  interface collate
    module subroutine coll_i1_0(scalarv_out, scalar_in)
      integer(int8), intent(inout) :: scalarv_out(:)
      integer(int8), intent(in)  :: scalar_in
    end subroutine
    module subroutine coll_i4_0(scalarv_out, scalar_in)
      integer(int32), intent(inout) :: scalarv_out(:)
      integer(int32), intent(in)  :: scalar_in
    end subroutine
    module subroutine coll_i8_0(scalarv_out, scalar_in)
      integer(int64), intent(inout) :: scalarv_out(:)
      integer(int64), intent(in)  :: scalar_in
    end subroutine
    module subroutine coll_r4_0(scalarv_out, scalar_in)
      real(real32), intent(inout) :: scalarv_out(:)
      real(real32), intent(in)  :: scalar_in
    end subroutine
    module subroutine coll_r8_0(scalarv_out, scalar_in)
      real(real64), intent(inout) :: scalarv_out(:)
      real(real64), intent(in)  :: scalar_in
    end subroutine
    module subroutine coll_log_0(scalarv_out, scalar_in)
      logical, intent(inout) :: scalarv_out(:)
      logical, intent(in)  :: scalar_in
    end subroutine
    module subroutine coll_char_0(scalarv_out, scalar_in)
      character(*), intent(inout) :: scalarv_out(:)
      character(*), intent(in)  :: scalar_in
    end subroutine
    module subroutine coll_i1_1(vector_out, vector_in)
      integer(int8), intent(inout) :: vector_out(:)
      integer(int8), intent(in)  :: vector_in(:)
    end subroutine
    module subroutine coll_i4_1(vector_out, vector_in)
      integer(int32), intent(inout) :: vector_out(:)
      integer(int32), intent(in)  :: vector_in(:)
    end subroutine
    module subroutine coll_i8_1(vector_out, vector_in)
      integer(int64), intent(inout) :: vector_out(:)
      integer(int64), intent(in)  :: vector_in(:)
    end subroutine
    module subroutine coll_r4_1(vector_out, vector_in)
      real(real32), intent(inout) :: vector_out(:)
      real(real32), intent(in)  :: vector_in(:)
    end subroutine
    module subroutine coll_r8_1(vector_out, vector_in)
      real(real64), intent(inout) :: vector_out(:)
      real(real64), intent(in)  :: vector_in(:)
    end subroutine
    module subroutine coll_log_1(vector_out, vector_in)
      logical, intent(inout) :: vector_out(:)
      logical, intent(in)  :: vector_in(:)
    end subroutine
    module subroutine coll_char_1(vector_out, vector_in)
      character(*), intent(inout) :: vector_out(:)
      character(*), intent(in)  :: vector_in(:)
    end subroutine
    module subroutine coll_i1_2(vector_out, vector_in)
      integer(int8), intent(inout) :: vector_out(:,:)
      integer(int8), intent(in)  :: vector_in(:,:)
    end subroutine
    module subroutine coll_i4_2(vector_out, vector_in)
      integer(int32), intent(inout) :: vector_out(:,:)
      integer(int32), intent(in)  :: vector_in(:,:)
    end subroutine
    module subroutine coll_i8_2(vector_out, vector_in)
      integer(int64), intent(inout) :: vector_out(:,:)
      integer(int64), intent(in)  :: vector_in(:,:)
    end subroutine
    module subroutine coll_r4_2(vector_out, vector_in)
      real(real32), intent(inout) :: vector_out(:,:)
      real(real32), intent(in)  :: vector_in(:,:)
    end subroutine
    module subroutine coll_r8_2(vector_out, vector_in)
      real(real64), intent(inout) :: vector_out(:,:)
      real(real64), intent(in)  :: vector_in(:,:)
    end subroutine
    module subroutine coll_log_2(vector_out, vector_in)
      logical, intent(inout) :: vector_out(:,:)
      logical, intent(in)  :: vector_in(:,:)
    end subroutine
    module subroutine coll_char_2(vector_out, vector_in)
      character(*), intent(inout) :: vector_out(:,:)
      character(*), intent(in)  :: vector_in(:,:)
    end subroutine
  end interface

  interface global_any
    module function any_0(mask) result(global)
      logical, intent(in) :: mask
      logical :: global
    end function
    module function any_1(mask) result(global)
      logical, intent(in) :: mask(:)
      logical :: global
    end function
    module function any_2(mask) result(global)
      logical, intent(in) :: mask(:,:)
      logical :: global
    end function
  end interface

  interface global_all
    module function all_0(mask) result(global)
      logical, intent(in) :: mask
      logical :: global
    end function
    module function all_1(mask) result(global)
      logical, intent(in) :: mask(:)
      logical :: global
    end function
    module function all_2(mask) result(global)
      logical, intent(in) :: mask(:,:)
      logical :: global
    end function
  end interface

  interface global_count
    module function count_0(mask) result(n)
      logical, intent(in) :: mask
      integer :: n
    end function
    module function count_1(mask) result(n)
      logical, intent(in) :: mask(:)
      integer :: n
    end function
    module function count_2(mask) result(n)
      logical, intent(in) :: mask(:,:)
      integer :: n
    end function
  end interface

  interface global_sum
    module function sum_i4_0(a) result(s)
      integer(int32), intent(in) :: a
      integer(int32) :: s
    end function
    module function sum_i8_0(a) result(s)
      integer(int64), intent(in) :: a
      integer(int64) :: s
    end function
    module function sum_r4_0(a) result(s)
      real(real32), intent(in) :: a
      real(real32) :: s
    end function
    module function sum_r8_0(a) result(s)
      real(real64), intent(in) :: a
      real(real64) :: s
    end function
    module function sum_i4_1(a, mask) result(s)
      integer(int32), intent(in) :: a(:)
      logical, intent(in), optional :: mask(:)
      integer(int32) :: s
    end function
    module function sum_i8_1(a, mask) result(s)
      integer(int64), intent(in) :: a(:)
      logical, intent(in), optional :: mask(:)
      integer(int64) :: s
    end function
    module function sum_r4_1(a, mask) result(s)
      real(real32), intent(in) :: a(:)
      logical, intent(in), optional :: mask(:)
      real(real32) :: s
    end function
    module function sum_r8_1(a, mask) result(s)
      real(real64), intent(in) :: a(:)
      logical, intent(in), optional :: mask(:)
      real(real64) :: s
    end function
  end interface

  interface global_minval
    module function minval_i4_0(a) result(v)
      integer(int32), intent(in) :: a
      integer(int32) :: v
    end function
    module function minval_i8_0(a) result(v)
      integer(int64), intent(in) :: a
      integer(int64) :: v
    end function
    module function minval_r4_0(a) result(v)
      real(real32), intent(in) :: a
      real(real32) :: v
    end function
    module function minval_r8_0(a) result(v)
      real(real64), intent(in) :: a
      real(real64) :: v
    end function
    module function minval_i4_1(a, mask) result(v)
      integer(int32), intent(in) :: a(:)
      logical, intent(in), optional :: mask(:)
      integer(int32) :: v
    end function
    module function minval_i8_1(a, mask) result(v)
      integer(int64), intent(in) :: a(:)
      logical, intent(in), optional :: mask(:)
      integer(int64) :: v
    end function
    module function minval_r4_1(a, mask) result(v)
      real(real32), intent(in) :: a(:)
      logical, intent(in), optional :: mask(:)
      real(real32) :: v
    end function
    module function minval_r8_1(a, mask) result(v)
      real(real64), intent(in) :: a(:)
      logical, intent(in), optional :: mask(:)
      real(real64) :: v
    end function
  end interface

  interface global_maxval
    module function maxval_i4_0(a) result(v)
      integer(int32), intent(in) :: a
      integer(int32) :: v
    end function
    module function maxval_i8_0(a) result(v)
      integer(int64), intent(in) :: a
      integer(int64) :: v
    end function
    module function maxval_r4_0(a) result(v)
      real(real32), intent(in) :: a
      real(real32) :: v
    end function
    module function maxval_r8_0(a) result(v)
      real(real64), intent(in) :: a
      real(real64) :: v
    end function
    module function maxval_i4_1(a, mask) result(v)
      integer(int32), intent(in) :: a(:)
      logical, intent(in), optional :: mask(:)
      integer(int32) :: v
    end function
    module function maxval_i8_1(a, mask) result(v)
      integer(int64), intent(in) :: a(:)
      logical, intent(in), optional :: mask(:)
      integer(int64) :: v
    end function
    module function maxval_r4_1(a, mask) result(v)
      real(real32), intent(in) :: a(:)
      logical, intent(in), optional :: mask(:)
      real(real32) :: v
    end function
    module function maxval_r8_1(a, mask) result(v)
      real(real64), intent(in) :: a(:)
      logical, intent(in), optional :: mask(:)
      real(real64) :: v
    end function
  end interface

  interface global_dot_product
    module function dot_prod_r4(a, b) result(dp)
      real(real32), intent(in) :: a(:), b(:)
      real(real32) :: dp
    end function
    module function dot_prod_r8(a, b) result(dp)
      real(real64), intent(in) :: a(:), b(:)
      real(real64) :: dp
    end function
  end interface

  logical :: initialized = .false., flag = .false.

contains

  subroutine init_parallel_communication
    use pgslib_module
    integer :: ierr
    if (initialized) return ! should only be called once
    initialized = .true.
    call MPI_Initialized(flag, ierr)
    if (.not.flag) call MPI_Init(ierr)
    call PGSLib_INITIALIZE()
    call MPI_Comm_size(comm, npe, ierr)
    call MPI_Comm_rank(comm, this_pe, ierr)
    this_pe = this_pe + 1 ! start numbering at 1
    is_iop = (this_pe == io_pe)
  end subroutine

  subroutine halt_parallel_communication
    use pgslib_module
    integer :: ierr
    call PGSLib_finalize
    if (.not.flag) call MPI_Finalize(ierr)
  end subroutine

  subroutine abort_parallel_communication
    integer :: ierr
    call MPI_Abort(comm, 1, ierr)
  end subroutine

  !! Used in only two locations in edit_module.F90. The PGSLib C implementation
  !! used MPI_MAXLOC with the MPI_DOUBLE_INT type, which isn't available in
  !! the Fortran interface. We could fudge things by using MPI_MAXLOC with
  !! the MPI_2DOUBLE_PRECISION (and sticking the integer index into a double).
  !! Instead we do this manually.

  function global_maxloc(a) result(gid)
    real(real64), intent(in) :: a(:)
    integer :: gid(1), last_gid, local_gid, ierr
    real(real64) :: local_max, global_max
    local_max = maxval(a)
    call MPI_Allreduce(local_max, global_max, 1, MPI_REAL8, MPI_MAX, comm, ierr)
    call MPI_Scan(size(a), last_gid, 1, MPI_INTEGER, MPI_SUM, comm, ierr)
    if (size(a) > 0 .and. local_max == global_max) then
      local_gid = last_gid - size(a) + maxloc(a,dim=1)
    else
      local_gid = huge(1) ! invalid
    end if
    call MPI_Allreduce(local_gid, gid, 1, MPI_INTEGER, MPI_MIN, comm, ierr)
    if (gid(1) == huge(1)) gid(1) = 0 ! 0-sized arrays on all ranks
  end function

  function global_minloc(a) result(gid)
    real(real64), intent(in) :: a(:)
    integer :: gid(1), last_gid, local_gid, ierr
    real(real64) :: local_min, global_min
    local_min = minval(a)
    call MPI_Allreduce(local_min, global_min, 1, MPI_REAL8, MPI_MIN, comm, ierr)
    call MPI_Scan(size(a), last_gid, 1, MPI_INTEGER, MPI_SUM, comm, ierr)
    if (size(a) > 0 .and. local_min == global_min) then
      local_gid = last_gid - size(a) + minloc(a,dim=1)
    else
      local_gid = huge(1) ! invalid
    end if
    call MPI_Allreduce(local_gid, gid, 1, MPI_INTEGER, MPI_MIN, comm, ierr)
    if (gid(1) == huge(1)) gid(1) = 0 ! 0-sized arrays on all ranks
  end function

  !! Used only in FHT_norm_type.F90.

  subroutine global_maxloc_sub(array, pid, lindex, mask)

    real(real64), intent(in) :: array(:)
    integer, intent(out) :: pid, lindex
    logical, intent(in), optional :: mask(:)

    integer :: n, local_pair(2), global_pair(2), ierr
    real(real64) :: local_max, global_max

    if (present(mask)) then
      ASSERT(size(mask) == size(array))
      n = count(mask)
    else
      n = size(array)
    end if

#ifdef INTEL_BUG20220211
    if (present(mask)) then
      local_max = maxval(array, mask)
    else
      local_max = maxval(array)
    end if
#else
    local_max = maxval(array, mask)
#endif
    call MPI_Allreduce(local_max, global_max, 1, MPI_REAL8, MPI_MAX, comm, ierr)

    local_pair(1) = merge(this_pe, npe+1, (n > 0 .and. local_max == global_max))
#ifdef INTEL_BUG20220211
    if (present(mask)) then
      local_pair(2) = maxloc(array, dim=1, mask=mask) ! 0 if array is 0-sized
    else
      local_pair(2) = maxloc(array, dim=1 ) ! 0 if array is 0-sized
    end if
#else
    local_pair(2) = maxloc(array, dim=1, mask=mask) ! 0 if array is 0-sized
#endif
    call MPI_Allreduce(local_pair, global_pair, 1, MPI_2INTEGER, MPI_MINLOC, comm, ierr)

    pid = global_pair(1)
    if (pid > npe) pid = 0  ! 0-sized arrays on all ranks
    lindex = global_pair(2)

  end subroutine global_maxloc_sub

end module parallel_communication
