!!
!! SCL - Simple Communications Layer
!!

module scl

  implicit none
  private

  public :: scl_init, scl_finalize
  public :: scl_send, scl_recv
  public :: scl_bcast, scl_gather, scl_allgather, scl_scatter
  public :: scl_global_sum, scl_global_maxval
  public :: scl_size, scl_rank

  interface
    subroutine initialize ()
    end subroutine
    subroutine finalize ()
    end subroutine
    integer function f77_comm_rank()
    end function
    integer function f77_comm_size()
    end function
    subroutine MPI_Send_int (buf, count, dest, tag)
      integer, intent(in) :: buf(*), count, dest, tag
    end subroutine
    subroutine MPI_Send_float (buf, count, dest, tag)
      real, intent(in) :: buf(*)
      integer, intent(in) :: count, dest, tag
    end subroutine
    subroutine MPI_Send_double (buf, count, dest, tag)
      double precision, intent(in) :: buf(*)
      integer, intent(in) :: count, dest, tag
    end subroutine
    subroutine MPI_Recv_int (buf, count, source, tag)
      integer, intent(out) :: buf(*)
      integer, intent(in) :: count, source, tag
    end subroutine
    subroutine MPI_Recv_float (buf, count, source, tag)
      real, intent(out) :: buf(*)
      integer, intent(in) :: count, source, tag
    end subroutine
    subroutine MPI_Recv_double (buf, count, source, tag)
      double precision, intent(out) :: buf(*)
      integer, intent(in) :: count, source, tag
    end subroutine
    subroutine MPI_Bcast_int_scalar (buffer, root)
      integer, intent(inout) :: buffer
      integer, intent(in) :: root
    end subroutine
    subroutine MPI_Bcast_int_vector (buffer, len, root)
      integer, intent(inout) :: buffer(*)
      integer, intent(in) :: len, root
    end subroutine
    subroutine MPI_Bcast_logical_scalar (buffer, root)
      logical, intent(inout) :: buffer
      integer, intent(in) :: root
    end subroutine
    subroutine MPI_Bcast_logical_vector (buffer, len, root)
      logical, intent(inout) :: buffer(*)
      integer, intent(in) :: len, root
    end subroutine
    subroutine MPI_Bcast_float_scalar (buffer, root)
      real, intent(inout) :: buffer
      integer, intent(in) :: root
    end subroutine
    subroutine MPI_Bcast_float_vector (buffer, len, root)
      real, intent(inout) :: buffer(*)
      integer, intent(in) :: len, root
    end subroutine
    subroutine MPI_Bcast_double_scalar (buffer, root)
      double precision, intent(inout) :: buffer
      integer, intent(in) :: root
    end subroutine
    subroutine MPI_Bcast_double_vector (buffer, len, root)
      double precision, intent(inout) :: buffer(*)
      integer, intent(in) :: len, root
    end subroutine
    subroutine MPI_Bcast_char_vector (buffer, root)
      character(len=*), intent(inout) :: buffer
      integer, intent(in) :: root
    end subroutine
    subroutine MPI_Gather_int (sendbuf, recvbuf, root)
      integer, intent(in) :: sendbuf, root
      integer, intent(out) :: recvbuf(*)
    end subroutine
    subroutine MPI_Gatherv_int (sendbuf, sendcount, recvbuf, recvcounts, root)
      integer, intent(in) :: sendbuf(*), sendcount, recvcounts(*), root
      integer, intent(out) :: recvbuf(*)
    end subroutine
    subroutine MPI_Gather_float (sendbuf, recvbuf, root)
      real, intent(in) :: sendbuf
      real, intent(out) :: recvbuf(*)
      integer, intent(in) :: root
    end subroutine
    subroutine MPI_Gatherv_float (sendbuf, sendcount, recvbuf, recvcounts, root)
      real, intent(in) :: sendbuf(*)
      integer, intent(in) :: sendcount, recvcounts(*), root
      real, intent(out) :: recvbuf(*)
    end subroutine
    subroutine MPI_Gather_double (sendbuf, recvbuf, root)
      double precision, intent(in) :: sendbuf
      double precision, intent(out) :: recvbuf(*)
      integer, intent(in) :: root
    end subroutine
    subroutine MPI_Gatherv_double (sendbuf, sendcount, recvbuf, recvcounts, root)
      double precision, intent(in) :: sendbuf(*)
      integer, intent(in) :: sendcount, recvcounts(*), root
      double precision, intent(out) :: recvbuf(*)
    end subroutine
    subroutine MPI_Allgather_int (sendbuf, recvbuf)
      integer, intent(in) :: sendbuf
      integer, intent(out) :: recvbuf(*)
    end subroutine
    subroutine MPI_Allgatherv_int (sendbuf, sendcount, recvbuf, recvcounts)
      integer, intent(in) :: sendbuf(*), sendcount, recvcounts(*)
      integer, intent(out) :: recvbuf(*)
    end subroutine
    subroutine MPI_Allgather_float (sendbuf, recvbuf)
      real, intent(in) :: sendbuf
      real, intent(out) :: recvbuf(*)
    end subroutine
    subroutine MPI_Allgatherv_float (sendbuf, sendcount, recvbuf, recvcounts)
      real, intent(in) :: sendbuf(*)
      integer, intent(in) :: sendcount, recvcounts(*)
      real, intent(out) :: recvbuf(*)
    end subroutine
    subroutine MPI_Scatter_int (sendbuf, recvbuf, root)
      integer, intent(in) :: sendbuf(*), root
      integer, intent(out) :: recvbuf
    end subroutine
    subroutine MPI_Scatterv_int (sendbuf, sendcounts, recvbuf, recvcount, root)
      integer, intent(in) :: sendbuf(*), sendcounts(*), recvcount, root
      integer, intent(out) :: recvbuf(*)
    end subroutine
    subroutine MPI_Scatter_float (sendbuf, recvbuf, root)
      real, intent(in) :: sendbuf(*)
      real, intent(out) :: recvbuf
      integer, intent(in) :: root
    end subroutine
    subroutine MPI_Scatterv_float (sendbuf, sendcounts, recvbuf, recvcount, root)
      real, intent(in) :: sendbuf(*)
      integer, intent(in) :: sendcounts(*), recvcount, root
      real, intent(out) :: recvbuf(*)
    end subroutine
    subroutine MPI_Scatter_double (sendbuf, recvbuf, root)
      double precision, intent(in) :: sendbuf(*)
      double precision, intent(out) :: recvbuf
      integer, intent(in) :: root
    end subroutine
    subroutine MPI_Scatterv_double (sendbuf, sendcounts, recvbuf, recvcount, root)
      double precision, intent(in) :: sendbuf(*)
      integer, intent(in) :: sendcounts(*), recvcount, root
      double precision, intent(out) :: recvbuf(*)
    end subroutine
    subroutine MPI_Allreduce_sum_int (sendbuf, recvbuf)
      integer, intent(in) :: sendbuf
      integer, intent(out) :: recvbuf
    end subroutine
    subroutine MPI_Allreduce_sum_float (sendbuf, recvbuf)
      real, intent(in) :: sendbuf
      real, intent(out) :: recvbuf
    end subroutine
    subroutine MPI_Allreduce_sum_double (sendbuf, recvbuf)
      double precision, intent(in) :: sendbuf
      double precision, intent(out) :: recvbuf
    end subroutine
    subroutine MPI_Allreduce_max_int (sendbuf, recvbuf)
      integer, intent(in) :: sendbuf
      integer, intent(out) :: recvbuf
    end subroutine
    subroutine MPI_Allreduce_max_float (sendbuf, recvbuf)
      real, intent(in) :: sendbuf
      real, intent(out) :: recvbuf
    end subroutine
    subroutine MPI_Allreduce_max_double (sendbuf, recvbuf)
      double precision, intent(in) :: sendbuf
      double precision, intent(out) :: recvbuf
    end subroutine
  end interface

  interface scl_bcast
    module procedure bcast_i0, bcast_i1, bcast_i2
    module procedure bcast_l0, bcast_l1, bcast_l2
    module procedure bcast_r0, bcast_r1, bcast_r2
    module procedure bcast_d0, bcast_d1, bcast_d2
    module procedure bcast_c0, bcast_c1
  end interface

  interface scl_gather
    module procedure gather_i0, gather_i1
    module procedure gather_r0, gather_r1
    module procedure gather_d0, gather_d1
  end interface

  interface scl_allgather
    module procedure allgather_i0, allgather_i1
    module procedure allgather_r0, allgather_r1
  end interface

  interface scl_scatter
    module procedure scatter_i0, scatter_i1
    module procedure scatter_r0, scatter_r1
    module procedure scatter_d0, scatter_d1
  end interface

  interface scl_send
    module procedure send_i1
    module procedure send_r1
    module procedure send_d1
  end interface

  interface scl_recv
    module procedure recv_i1
    module procedure recv_r1
    module procedure recv_d1
  end interface

  interface scl_global_sum
    module procedure global_sum_i0, global_sum_i1
    module procedure global_sum_r0, global_sum_r1
    module procedure global_sum_d0, global_sum_d1
  end interface

  interface scl_global_maxval
    module procedure global_max_i0, global_max_i1
    module procedure global_max_r0, global_max_r1
    module procedure global_max_d0, global_max_d1
  end interface

  integer, save :: comm_size = -1, comm_rank = -1

contains

  subroutine scl_init ()
    call initialize ()
    comm_size = f77_comm_size()
    comm_rank = f77_comm_rank() + 1
  end subroutine scl_init

  subroutine scl_finalize ()
    call finalize ()
  end subroutine scl_finalize

  pure integer function scl_rank ()
    scl_rank = comm_rank
  end function

  pure integer function scl_size ()
    scl_size = comm_size
  end function

  subroutine bcast_i0 (x, root)
    integer, intent(inout) :: x
    integer, intent(in), optional :: root
    if (present(root)) then
      call MPI_Bcast_int_scalar (x, root-1)
    else
      call MPI_Bcast_int_scalar (x, 0)
    endif
  end subroutine

  subroutine bcast_i1 (x, root)
    integer, intent(inout) :: x(:)
    integer, intent(in), optional :: root
    if (present(root)) then
      call MPI_Bcast_int_vector (x, size(x), root-1)
    else
      call MPI_Bcast_int_vector (x, size(x), 0)
    endif
  end subroutine

  subroutine bcast_i2 (x, root)
    integer, intent(inout) :: x(:,:)
    integer, intent(in), optional :: root
    if (present(root)) then
      call MPI_Bcast_int_vector (x, size(x), root-1)
    else
      call MPI_Bcast_int_vector (x, size(x), 0)
    endif
  end subroutine

  subroutine bcast_r0 (x, root)
    real, intent(inout) :: x
    integer, intent(in), optional :: root
    if (present(root)) then
      call MPI_Bcast_float_scalar (x, root-1)
    else
      call MPI_Bcast_float_scalar (x, 0)
    endif
  end subroutine

  subroutine bcast_r1 (x, root)
    real, intent(inout) :: x(:)
    integer, intent(in), optional :: root
    if (present(root)) then
      call MPI_Bcast_float_vector (x, size(x), root-1)
    else
      call MPI_Bcast_float_vector (x, size(x), 0)
    endif
  end subroutine

  subroutine bcast_r2 (x, root)
    real, intent(inout) :: x(:,:)
    integer, intent(in), optional :: root
    if (present(root)) then
      call MPI_Bcast_float_vector (x, size(x), root-1)
    else
      call MPI_Bcast_float_vector (x, size(x), 0)
    endif
  end subroutine

  subroutine bcast_d0 (x, root)
    double precision, intent(inout) :: x
    integer, intent(in), optional :: root
    if (present(root)) then
      call MPI_Bcast_double_scalar (x, root-1)
    else
      call MPI_Bcast_double_scalar (x, 0)
    endif
  end subroutine

  subroutine bcast_d1 (x, root)
    double precision, intent(inout) :: x(:)
    integer, intent(in), optional :: root
    if (present(root)) then
      call MPI_Bcast_double_vector (x, size(x), root-1)
    else
      call MPI_Bcast_double_vector (x, size(x), 0)
    endif
  end subroutine

  subroutine bcast_d2 (x, root)
    double precision, intent(inout) :: x(:,:)
    integer, intent(in), optional :: root
    if (present(root)) then
      call MPI_Bcast_double_vector (x, size(x), root-1)
    else
      call MPI_Bcast_double_vector (x, size(x), 0)
    endif
  end subroutine

  subroutine bcast_l0 (x, root)
    logical, intent(inout) :: x
    integer, intent(in), optional :: root
    if (present(root)) then
      call MPI_Bcast_logical_scalar (x, root-1)
    else
      call MPI_Bcast_logical_scalar (x, 0)
    endif
  end subroutine

  subroutine bcast_l1 (x, root)
    logical, intent(inout) :: x(:)
    integer, intent(in), optional :: root
    if (present(root)) then
      call MPI_Bcast_logical_vector (x,size(x),root-1)
    else
      call MPI_Bcast_logical_vector (x,size(x),0)
    endif
  end subroutine

  subroutine bcast_l2 (x, root)
    logical, intent(inout) :: x(:,:)
    integer, intent(in), optional :: root
    if (present(root)) then
      call MPI_Bcast_logical_vector (x,size(x),root-1)
    else
      call MPI_Bcast_logical_vector (x,size(x),0)
    endif
  end subroutine

  subroutine bcast_c0 (x, root)
    character(len=*), intent(inout) :: x
    integer, intent(in), optional :: root
    if (present(root)) then
      call MPI_Bcast_char_vector (x, root-1)
    else
      call MPI_Bcast_char_vector (x, 0)
    endif
  end subroutine

  subroutine bcast_c1 (x, root)
    character(len=*), intent(inout) :: x(:)
    integer :: j
    integer, intent(in), optional :: root
    if (present(root)) then
      do j = 1, size(x)
        call MPI_Bcast_char_vector (x(j), root-1)
      end do
    else
      do j = 1, size(x)
        call MPI_Bcast_char_vector (x(j), 0)
      end do
    endif
  end subroutine

  subroutine gather_i0 (src, dest, root)
    integer, intent(in)  :: src
    integer, intent(out) :: dest(:)
    integer, intent(in), optional :: root
    if (present(root)) then
      call MPI_Gather_int (src, dest, root-1)
    else
      call MPI_Gather_int (src, dest, 0)
    endif
  end subroutine

  subroutine gather_i1 (src, dest, root)
    integer, intent(in)  :: src(:)
    integer, intent(out) :: dest(:)
    integer, intent(in), optional :: root
    integer :: sizes(comm_size)
    if (present(root)) then
      call MPI_Gather_int (size(src), sizes, root-1)
      call MPI_Gatherv_int (src, size(src), dest, sizes, root-1)
    else
      call MPI_Gather_int (size(src), sizes, 0)
      call MPI_Gatherv_int (src, size(src), dest, sizes, 0)
    endif
  end subroutine

  subroutine gather_r0 (src, dest, root)
    real, intent(in)  :: src
    real, intent(out) :: dest(:)
    integer, intent(in), optional :: root
    if (present(root)) then
      call MPI_Gather_float (src, dest, root-1)
    else
      call MPI_Gather_float (src, dest, 0)
    endif
  end subroutine

  subroutine gather_r1 (src, dest, root)
    real, intent(in)  :: src(:)
    real, intent(out) :: dest(:)
    integer, intent(in), optional :: root
    integer :: sizes(comm_size)
    if (present(root)) then
      call MPI_Gather_int (size(src), sizes, root-1)
      call MPI_Gatherv_float (src, size(src), dest, sizes, root-1)
    else
      call MPI_Gather_int (size(src), sizes, 0)
      call MPI_Gatherv_float (src, size(src), dest, sizes, 0)
    endif
  end subroutine

  subroutine gather_d0 (src, dest, root)
    double precision, intent(in)  :: src
    double precision, intent(out) :: dest(:)
    integer, intent(in), optional :: root
    if (present(root)) then
      call MPI_Gather_double (src, dest, root-1)
    else
      call MPI_Gather_double (src, dest, 0)
    endif
  end subroutine

  subroutine gather_d1 (src, dest, root)
    double precision, intent(in)  :: src(:)
    double precision, intent(out) :: dest(:)
    integer, intent(in), optional :: root
    integer :: sizes(comm_size)
    if (present(root)) then
      call MPI_Gather_int (size(src), sizes, root-1)
      call MPI_Gatherv_double (src, size(src), dest, sizes, root-1)
    else
      call MPI_Gather_int (size(src), sizes, 0)
      call MPI_Gatherv_double (src, size(src), dest, sizes, 0)
    endif
  end subroutine

  subroutine allgather_i0 (src, dest)
    integer, intent(in)  :: src
    integer, intent(out) :: dest(:)
    call MPI_Allgather_int (src, dest)
  end subroutine

  subroutine allgather_i1 (src, dest)
    integer, intent(in)  :: src(:)
    integer, intent(out) :: dest(:)
    integer :: sizes(comm_size)
    call MPI_Allgather_int (size(src), sizes)
    call MPI_Allgatherv_int (src, size(src), dest, sizes)
  end subroutine

  subroutine allgather_r0 (src, dest)
    real, intent(in)  :: src
    real, intent(out) :: dest(:)
    call MPI_Allgather_float (src, dest)
  end subroutine

  subroutine allgather_r1 (src, dest)
    real, intent(in)  :: src(:)
    real, intent(out) :: dest(:)
    integer :: sizes(comm_size)
    call MPI_Allgather_int (size(src), sizes)
    call MPI_Allgatherv_float (src, size(src), dest, sizes)
  end subroutine

  subroutine scatter_i0 (src, dest, root)
    integer, intent(in)  :: src(:)
    integer, intent(out) :: dest
    integer, intent(in), optional :: root
    if (present(root)) then
      call MPI_Scatter_int (src, dest, root-1)
    else
      call MPI_Scatter_int (src, dest, 0)
    endif
  end subroutine

  subroutine scatter_i1 (src, dest, root)
    integer, intent(in)  :: src(:)
    integer, intent(out) :: dest(:)
    integer, intent(in), optional :: root
    integer :: sizes(comm_size)
    if (present(root)) then
      call MPI_Gather_int (size(dest), sizes, root-1)
      call MPI_Scatterv_int (src, sizes, dest, size(dest), root-1)
    else
      call MPI_Gather_int (size(dest), sizes, 0)
      call MPI_Scatterv_int (src, sizes, dest, size(dest), 0)
    endif
  end subroutine

  subroutine scatter_r0 (src, dest, root)
    real, intent(in)  :: src(:)
    real, intent(out) :: dest
    integer, intent(in), optional :: root
    if (present(root)) then
      call MPI_Scatter_float (src, dest, root-1)
    else
      call MPI_Scatter_float (src, dest, 0)
    endif
  end subroutine

  subroutine scatter_r1 (src, dest, root)
    real, intent(in)  :: src(:)
    real, intent(out) :: dest(:)
    integer, intent(in), optional :: root
    integer :: sizes(comm_size)
    if (present(root)) then
      call MPI_Gather_int (size(dest), sizes, root-1)
      call MPI_Scatterv_float (src, sizes, dest, size(dest), root-1)
    else
      call MPI_Gather_int (size(dest), sizes, 0)
      call MPI_Scatterv_float (src, sizes, dest, size(dest), 0)
    endif
  end subroutine

  subroutine scatter_d0 (src, dest, root)
    double precision, intent(in)  :: src(:)
    double precision, intent(out) :: dest
    integer, intent(in), optional :: root
    if (present(root)) then
      call MPI_Scatter_double (src, dest, root-1)
    else
      call MPI_Scatter_double (src, dest, 0)
    endif
  end subroutine

  subroutine scatter_d1 (src, dest, root)
    double precision, intent(in)  :: src(:)
    double precision, intent(out) :: dest(:)
    integer, intent(in), optional :: root
    integer :: sizes(comm_size)
    if (present(root)) then
      call MPI_Gather_int (size(dest), sizes, root-1)
      call MPI_Scatterv_double (src, sizes, dest, size(dest), root-1)
    else
      call MPI_Gather_int (size(dest), sizes, 0)
      call MPI_Scatterv_double (src, sizes, dest, size(dest), 0)
    endif
  end subroutine

  subroutine send_i1 (buf, dest, tag)
    integer, intent(in) :: buf(:)
    integer, intent(in) :: dest, tag
    call MPI_send_int (buf, size(buf), dest-1, tag)
  end subroutine

  subroutine send_r1 (buf, dest, tag)
    real,    intent(in) :: buf(:)
    integer, intent(in) :: dest, tag
    call MPI_send_float (buf, size(buf), dest-1, tag)
  end subroutine

  subroutine send_d1 (buf, dest, tag)
    double precision, intent(in) :: buf(:)
    integer, intent(in) :: dest, tag
    call MPI_send_double (buf, size(buf), dest-1, tag)
  end subroutine

  subroutine recv_i1 (buf, source, tag)
    integer, intent(out) :: buf(:)
    integer, intent(in) :: source, tag
    call MPI_recv_int (buf, size(buf), source-1, tag)
  end subroutine

  subroutine recv_r1 (buf, source, tag)
    real, intent(out) :: buf(:)
    integer, intent(in) :: source, tag
    call MPI_recv_float (buf, size(buf), source-1, tag)
  end subroutine

  subroutine recv_d1 (buf, source, tag)
    double precision, intent(out) :: buf(:)
    integer, intent(in) :: source, tag
    call MPI_recv_double (buf, size(buf), source-1, tag)
  end subroutine

  function global_sum_i0 (a) result (s)
    integer, intent(in) :: a
    integer :: s
    call MPI_Allreduce_sum_int (a, s)
  end function

  function global_sum_i1 (a) result (s)
    integer, intent(in) :: a(:)
    integer :: s
    call MPI_Allreduce_sum_int (sum(a), s)
  end function

  function global_sum_r0 (a) result (s)
    real, intent(in) :: a
    real :: s
    call MPI_Allreduce_sum_float (a, s)
  end function

  function global_sum_r1 (a) result (s)
    real, intent(in) :: a(:)
    real :: s
    call MPI_Allreduce_sum_float (sum(a), s)
  end function

  function global_sum_d0 (a) result (s)
    double precision, intent(in) :: a
    double precision :: s
    call MPI_Allreduce_sum_double (a, s)
  end function

  function global_sum_d1 (a) result (s)
    double precision, intent(in) :: a(:)
    double precision :: s
    call MPI_Allreduce_sum_double (sum(a), s)
  end function

  function global_max_i0 (a) result (s)
    integer, intent(in) :: a
    integer :: s
    call MPI_Allreduce_max_int (a, s)
  end function

  function global_max_i1 (a) result (s)
    integer, intent(in) :: a(:)
    integer :: s
    call MPI_Allreduce_max_int (maxval(a), s)
  end function

  function global_max_r0 (a) result (s)
    real, intent(in) :: a
    real :: s
    call MPI_Allreduce_max_float (a, s)
  end function

  function global_max_r1 (a) result (s)
    real, intent(in) :: a(:)
    real :: s
    call MPI_Allreduce_max_float (maxval(a), s)
  end function

  function global_max_d0 (a) result (s)
    double precision, intent(in) :: a
    double precision :: s
    call MPI_Allreduce_max_double (a, s)
  end function

  function global_max_d1 (a) result (s)
    double precision, intent(in) :: a(:)
    double precision :: s
    call MPI_Allreduce_max_double (maxval(a), s)
  end function

end module scl
