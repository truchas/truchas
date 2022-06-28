!!
!! RE_DIST_VF_TYPE
!!
!! This module provides a derived type for describing the distributed view
!! factor data of a radiation enclosure and methods that operate on instances
!! of this type.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! 3 Apr 2008
!!
!! David Neill-Asanza <dhna@lanl.gov>
!! 6 Aug 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  CALL READ_DIST_VF (THIS, PATH) initializes THIS with the view factor
!!    data read from the radiation enclosure dataset PATH.  The dataset
!!    must contain view factor data.  This is a collective procedure and
!!    must be called from all process ranks.
!!
!!  CALL WRITE_DIST_VF (THIS, PATH) writes the view factor data contained in
!!    THIS to the radiation enclosure dataset PATH.  The data set must already
!!    contain compatible enclosure data but no view factor data.  This is a
!!    collective procedure and must be called from all process ranks.
!!


#include "f90_assert.fpp"

module re_dist_vf_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64, i8 => int64
  use scl
  implicit none
  private

  type, public :: dist_vf
    integer :: npatch = 0     ! number of patches (number of matrix rows) on this process
    integer :: offset = 0     ! difference between the local row index and global index
    integer :: npatch_tot = 0 ! total number of patches (number of matrix columns)
    !! View factor matrix in CSR format.
    real,     allocatable :: val(:)   ! values of the nonzero matrix elements stored by row
    integer,  allocatable :: ja(:)    ! column indices for the matrix elements in VAL
    integer,  allocatable :: ia(:)    ! start index in JA and VAL for the rows
    real(r8), allocatable :: area(:)  ! patch areas
    real(r8), allocatable :: w(:)     ! ratio of face area to patch area (face weights)
    real,     allocatable :: ambient(:)  ! ambient view factors
    logical :: has_ambient
    logical :: has_area
    logical :: has_weight
  contains
    procedure :: write
    procedure :: read
    procedure :: unpack_row
    procedure :: row_sum
    procedure :: get_ambient_vf
    procedure :: unpack_col
    procedure :: unpack_col_explicit, unpack_col_from_row
  end type

contains

  subroutine write (this, path)

    use rad_encl_file_type

    class(dist_vf), intent(in) :: this
    character(len=*), intent(in) :: path

    type(rad_encl_file) :: file
    integer :: n, nproc, my_rank, bsize(scl_size())
    integer, allocatable :: rowcount(:), ibuf(:)
    real, allocatable :: ambient(:), rbuf(:)
    integer(i8) :: start

    nproc   = scl_size()
    my_rank = scl_rank()

    if (nproc == 1) then

      call file%open_rw(path)
      call file%init_vf(size(this%val,kind=i8), this%has_ambient)
      call file%put_vf_rowcount(this%ia(2:)-this%ia(:size(this%ia)-1))
      call file%put_vf_rows(this%val, this%ja, start=1_i8)
      if (this%has_ambient) call file%put_ambient(this%ambient)
      if (this%has_area) call file%put_area(this%area)
      if (this%has_weight) call file%put_face_weight(this%w)
      call file%close

    else

      !! Gather VF matrix data
      n = merge(this%npatch_tot, 0, my_rank==1)
      allocate(rowcount(n))
      call scl_gather(this%ia(2:)-this%ia(:this%npatch), rowcount)
      call scl_allgather(size(this%val), bsize)

      !! Gather ambient view factors
      n = merge(this%npatch_tot, 0, my_rank==1 .and. this%has_ambient)
      allocate(ambient(n))
      if (this%has_ambient) call scl_gather(this%ambient, ambient)

      if (my_rank == 1) then
        call file%open_rw(path)
        call file%init_vf(sum(int(bsize,kind=i8)), this%has_ambient)
        call file%put_vf_rowcount(rowcount)
        if (this%has_ambient) call file%put_ambient(ambient)
        if (this%has_area) call file%put_area(this%area)
        if (this%has_weight) call file%put_face_weight(this%w)
      end if
      deallocate(rowcount, ambient)

      !! Write the VF matrix in process-sized blocks,
      !! receiving it from the owning process as we go.
      if (my_rank > 1) then
        call scl_send(this%val, dest=1, tag=my_rank)
        call scl_send(this%ja,  dest=1, tag=my_rank)
      else
        call file%put_vf_rows(this%val, this%ja, start=1_i8)
        n = maxval(bsize)
        allocate(ibuf(n), rbuf(n))
        start = 1 + bsize(1)
        do n = 2, nproc
          call scl_recv(rbuf(:bsize(n)), source=n, tag=n)
          call scl_recv(ibuf(:bsize(n)), source=n, tag=n)
          call file%put_vf_rows(rbuf(:bsize(n)), ibuf(:bsize(n)), start)
          start = start + bsize(n)
        end do
        deallocate(ibuf, rbuf)
        call file%close
      end if

    end if

  end subroutine write

  subroutine read (this, path)

    use rad_encl_file_type

    class(dist_vf), intent(out) :: this
    character(len=*), intent(in) :: path

    type(rad_encl_file) :: file
    integer :: j, n, nproc, my_rank, npatch, npatch_tot, nface_tot, bsize(scl_size())
    integer, allocatable :: rowcount(:), ibuf(:)
    real, allocatable :: ambient(:), rbuf(:)
    integer(i8) :: nnonz, start

    nproc   = scl_size()
    my_rank = scl_rank()

    if (nproc == 1) then

      call file%open_ro(path)
      call file%get_vf_dims(nface_tot, npatch_tot, nnonz)

      this%has_ambient = file%has_ambient()
      this%has_area = file%has_area()
      this%has_weight = file%has_face_weight()
      this%npatch = npatch_tot
      this%offset = 0
      this%npatch_tot = npatch_tot

      !! Read patch areas
      if (this%has_area) then
        allocate(this%area(npatch_tot))
        call file%get_area(this%area)
      end if

      !! Read face weights
      if (this%has_weight) then
        allocate(this%w(nface_tot))
        call file%get_face_weight(this%w)
      end if

      !! Read VF matrix
      allocate(this%ia(npatch_tot+1), this%val(nnonz), this%ja(nnonz))
      call file%get_vf_rowcount(this%ia(2:))
      call file%get_vf_rows(this%val, this%ja, start=1_i8)

      !! Read ambient view factors
      if (this%has_ambient) then
        allocate(this%ambient(npatch_tot))
        call file%get_ambient(this%ambient)
      end if

      call file%close

      this%ia(1) = 1
      do j = 2, ubound(this%ia,1)
        this%ia(j) = this%ia(j) + this%ia(j-1)
      end do

    else

      if (my_rank == 1) then
        call file%open_ro(path)
        call file%get_vf_dims(nface_tot, npatch_tot, nnonz)
        this%has_ambient = file%has_ambient()
        this%has_area = file%has_area()
        this%has_weight = file%has_face_weight()
      end if

      call scl_bcast(npatch_tot)
      call scl_bcast(this%has_ambient)
      call scl_bcast(this%has_area)
      call scl_bcast(this%has_weight)

      this%npatch_tot = npatch_tot

      !! Read and distribute the patch areas
      if (this%has_area) then
        allocate(this%area(npatch_tot))
        if (my_rank == 1) call file%get_area(this%area)
        call scl_bcast(this%area)
      end if

      !! Read the face weights. Only needed on rank 1.
      if (this%has_weight .and. my_rank==1) then
        allocate(this%w(nface_tot))
        call file%get_face_weight(this%w)
      end if

      !! Divvy up the rows.
      npatch = npatch_tot/nproc
      if (my_rank <= modulo(npatch_tot,nproc)) npatch = npatch + 1
      ASSERT(scl_global_sum(npatch) == npatch_tot)
      this%npatch = npatch

      call scl_allgather(npatch, bsize)
      this%offset = sum(bsize(1:my_rank-1))

      n = merge(npatch_tot, 0, my_rank == 1)

      !! Read and distribute the ambient viewfactors
      if (this%has_ambient) then
        allocate(ambient(n), this%ambient(npatch))
        if (my_rank == 1) call file%get_ambient(ambient)
        call scl_scatter(ambient, this%ambient)
      end if

      !! Read the VF matrix row counts
      allocate(rowcount(n))
      if (my_rank == 1) call file%get_vf_rowcount(rowcount)

      !! Distribute the row counts; generate the
      !! local IA indexing array from the counts.
      allocate(this%ia(npatch+1))
      call scl_scatter(rowcount, this%ia(2:))
      this%ia(1) = 1
      do j = 2, ubound(this%ia,1)
        this%ia(j) = this%ia(j) + this%ia(j-1)
      end do

      !! Determine the sizes of the distributed VF matrix.
      n = this%ia(this%npatch+1) - this%ia(1)
      allocate(this%val(n), this%ja(n))
      call scl_allgather(n, bsize)

      !! Read the VF matrix in process-sized blocks,
      !! sending it to the owning process as we go.
      if (my_rank > 1) then
        call scl_recv(this%val, source=1, tag=my_rank)
        call scl_recv(this%ja,  source=1, tag=my_rank)
      else
        call file%get_vf_rows(this%val, this%ja, start=1_i8)
        n = maxval(bsize)
        allocate(ibuf(n), rbuf(n))
        start = 1 + bsize(1)
        do n = 2, nproc
          call file%get_vf_rows(rbuf(:bsize(n)), ibuf(:bsize(n)), start)
          call scl_send(rbuf(:bsize(n)), dest=n, tag=n)
          call scl_send(ibuf(:bsize(n)), dest=n, tag=n)
          start = start + bsize(n)
        end do
        deallocate(ibuf, rbuf)
        call file%close
      end if

    end if

  end subroutine read

  function get_ambient_vf (dvf) result (u)
    class(dist_vf), intent(in) :: dvf
    real, allocatable :: u(:)
    integer :: n

    if (dvf%has_ambient) then
      n = merge(dvf%npatch, 0, scl_rank() == 1)
      allocate(u(n))

      call scl_gather (dvf%ambient, u)
    else
      allocate(u(0))
    end if

  end function get_ambient_vf

  function unpack_row (dvf, n) result (u)

    class(dist_vf), intent(in) :: dvf
    integer, intent(in) :: n
    real, allocatable :: u(:)

    integer :: l, vsize
    integer, allocatable :: sidx(:)
    real, allocatable :: svec(:)

    l = n - dvf%offset  ! local index

    if (l >= 1 .and. l <= dvf%npatch) then ! this process has the row
      sidx = dvf%ja(dvf%ia(l):dvf%ia(l+1)-1)
      svec = dvf%val(dvf%ia(l):dvf%ia(l+1)-1)
      vsize = dvf%npatch_tot
    else
      allocate(sidx(0), svec(0))
      vsize = 0
    end if
    u = unpack_dist_sparse_vector (vsize, svec, sidx)

  end function unpack_row

  function row_sum (dvf) result (u)

    class(dist_vf), intent(in) :: dvf
    real, allocatable :: u(:)

    integer :: j
    real :: u_l(dvf%npatch)

    do j = 1, dvf%npatch
      u_l(j) = sum(dvf%val(dvf%ia(j):dvf%ia(j+1)-1))
    end do

    if (scl_rank() == 1) then
      allocate(u(dvf%npatch_tot))
    else
      allocate(u(0))
    end if
    call scl_gather (u_l, u)

  end function row_sum

  !! Returns the unpacked n-th column of the VF matrix on process rank 1.
  function unpack_col (dvf, n) result (u)

    class(dist_vf), intent(in) :: dvf
    integer, intent(in) :: n
    real, allocatable :: u(:)

    if (dvf%has_area) then
      u = dvf%unpack_col_from_row(n)
    else
      u = dvf%unpack_col_explicit(n)
    end if

  end function unpack_col

  !! Returns the unpacked n-th column of the VF matrix on process rank 1.
  !! The columns is formed implicitly from the n-th row using reciprocity.
  !! Thus, the only the process owning row n communicates with rank 1.
  function unpack_col_from_row (dvf, n) result (u)

    class(dist_vf), intent(in) :: dvf
    integer, intent(in) :: n
    real, allocatable :: u(:)

    integer :: l, vsize
    integer, allocatable :: sidx(:)
    real, allocatable :: svec(:)

    l = n - dvf%offset  ! local index

    if (l >= 1 .and. l <= dvf%npatch) then ! this process has the column
      sidx = dvf%ja(dvf%ia(l):dvf%ia(l+1)-1)
      svec = dvf%val(dvf%ia(l):dvf%ia(l+1)-1) * dvf%area(n) / dvf%area(sidx)
      vsize = dvf%npatch_tot
    else
      allocate(sidx(0), svec(0))
      vsize = 0
    end if
    u = unpack_dist_sparse_vector (vsize, svec, sidx)

  end function unpack_col_from_row

  !! Returns the unpacked n-th column of the VF matrix on process rank 1.
  !! The column is formed explicitly by having each rank send any non-zeros
  !! in the n-th column of each local row.
  function unpack_col_explicit (dvf, n) result (u)

    class(dist_vf), intent(in) :: dvf
    integer, intent(in) :: n
    real, allocatable :: u(:)

    integer :: i, j, cnt
    real :: svec(dvf%npatch)
    integer :: sidx(dvf%npatch)

    !! Gather up the nonzeros in column N
    cnt = 0
    do j = 1, dvf%npatch
      i = loc_in_row (n, dvf%ja(dvf%ia(j):dvf%ia(j+1)-1))
      if (i > 0) then
        i = i + dvf%ia(j) - 1
        cnt = cnt + 1
        svec(cnt) = dvf%val(i)
        sidx(cnt) = j + dvf%offset
      end if
    end do

    u = unpack_dist_sparse_vector (dvf%npatch, svec(:cnt), sidx(:cnt))

  end function unpack_col_explicit

  function unpack_dist_sparse_vector (vsize, svec, sidx) result (u)

    integer, intent(in) :: vsize    ! local size of the full local vector
    real,    intent(in) :: svec(:)  ! sparse vector values ...
    integer, intent(in) :: sidx(:)  ! and the corresponding global indices
    real, allocatable :: u(:)

    integer :: n
    integer, allocatable :: bsize(:), sidx_g(:)
    real,    allocatable :: svec_g(:)

    ASSERT( size(svec) == size(sidx) )
    ASSERT( size(svec) <= vsize )

    allocate(bsize(scl_size()))
    call scl_gather (size(svec), bsize)
    if (scl_rank() == 1) then
      n = sum(bsize)
      allocate (svec_g(n), sidx_g(n))
    else
      allocate (svec_g(0), sidx_g(0))
    end if

    call scl_gather (svec, svec_g)
    call scl_gather (sidx, sidx_g)

    call scl_gather (vsize, bsize)
    if (scl_rank() == 1) then
      n = sum(bsize)
      allocate(u(n))
      u = 0.0
      do n = 1, size(sidx_g)
        u(sidx_g(n)) = svec_g(n)
      end do
    else
      allocate(u(0))
    end if

    deallocate(bsize, svec_g, sidx_g)

  end function unpack_dist_sparse_vector


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! LOC_IN_ROW
 !!
 !! This auxillary procedure returns the index of the value N in the array LIST
 !! or 0 if LIST does not contain N.  It is assumed that the values of LIST are
 !! in (strictly) increasing order.  A simple binary search algorithm is used.
 !! In its current use, LIST are the column indices of the nonzero elements in
 !! a row of a matrix, and N is a specific column index.
 !!
 !! Note: this improves upon the naive sequential search previously used, but
 !! I think additional significant performance gains could be obtained by
 !! recognizing that, in some uses, this routine is called with N running
 !! sequentially from 1.
 !!

  integer function loc_in_row (n, list) result (loc)

    integer, intent(in) :: n, list(:)

    integer :: i, i1, i2

    loc = 0
    if (size(list) == 0) return

    i1 = 1
    if (n <= list(i1)) then
      if (n == list(i1)) loc = i1
      return
    end if

    i2 = size(list)
    if (n >= list(i2)) then
      if (n == list(i2)) loc = i2
      return
    end if

    do while (i2-i1 > 1)
      i = (i1+i2)/2
      if (n < list(i)) then
        i2 = i
      else if (n > list(i)) then
        i1 = i
      else
        loc = i
        exit
      endif
    end do

  end function loc_in_row

end module re_dist_vf_type
