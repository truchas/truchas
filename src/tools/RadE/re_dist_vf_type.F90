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
!!  CALL VF_DIFF (A, B) replaces the distributed VF matrix B with the
!!    difference B - A.  The matrices need not have the same sparsity
!!    pattern, but A and B must be identically distributed across processes.
!!    This is a collective procedure.
!!
!!  CALL VF_DIFF_MAX_NORM (D, DNORM) computes the max norm of the
!!    distributed VF difference matrix D as an operator.  If e is the
!!    irradiance discrepancy resulting from D and a given radiosity q,
!!    then |e|_max <= |D|_max |q|_max, where |D|_max is maximum absolute
!!    row sum of D.  This is a collective procedure.
!!
!!  CALL VF_DIFF_ONE_NORM (D, DNORM) computes the L1 norm of the
!!    distributed VF difference matrix D as an operator.  If e is the
!!    irradiance discrepancy resulting from D and a given radiosity q,
!!    then |e|_1 <= |D|_1 |q|_1, where |D|_1 is the sum of the maximum
!!    absolute element value in each column of D.  This is a collective
!!    procedure.
!!
!!  CALL VF_DIFF_TWO_NORM (D, DNORM) computes the L2 norm of the
!!    distributed VF difference matrix D as an operator.  If e is the
!!    irradiance discrepancy resulting from D and a given radiosity q,
!!    then |e|_2 <= |D|_2 |q|_2, where |D|_2 is the square root of the sum
!!    over j of the dot product of the jth row with the jth column of D.
!!    This is a collective procedure.
!!
!!  CALL VF_DIFF_MAX (D, VAL, RLOC, CLOC) finds the element of the VF
!!    difference matrix D of maximum absolute value VAL and its row RLOC
!!    and column CLOC location.  This is a collective procedure.
!!


#include "f90_assert.fpp"

module re_dist_vf_type

  use kinds, only: i8
  use scl
  implicit none
  private

  public :: read_dist_vf, write_dist_vf
  public :: unpack_dvf_row, dvf_row_sum, unpack_dvf_col, get_ambient_vf
  public :: vf_diff, vf_diff_one_norm, vf_diff_two_norm, vf_diff_max_norm, vf_diff_max

  type, public :: dist_vf
    integer :: npatch = 0     ! number of patches (number of matrix rows) on this process
    integer :: offset = 0     ! difference between the local row index and global index
    integer :: npatch_tot = 0 ! total number of patches (number of matrix columns)
    !! View factor matrix in CSR format.
    real,    allocatable :: val(:)   ! values of the nonzero matrix elements stored by row
    integer, allocatable :: ja(:)    ! column indices for the matrix elements in VAL
    integer, allocatable :: ia(:)    ! start index in JA and VAL for the rows
    real,    allocatable :: ambient(:)  ! ambient view factors
    logical :: has_ambient
  end type

contains

  subroutine write_dist_vf (this, path)

    use rad_encl_file_type

    type(dist_vf), intent(in) :: this
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

  end subroutine write_dist_vf

  subroutine read_dist_vf (this, path)

    use rad_encl_file_type

    type(dist_vf), intent(out) :: this
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
      this%npatch = npatch_tot
      this%offset = 0
      this%npatch_tot = npatch_tot

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
      end if

      call scl_bcast(npatch_tot)
      call scl_bcast(this%has_ambient)

      this%npatch_tot = npatch_tot

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

  end subroutine read_dist_vf

  function get_ambient_vf (dvf) result (u)
    type(dist_vf), intent(in) :: dvf
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

  function unpack_dvf_row (dvf, n) result (u)

    type(dist_vf), intent(in) :: dvf
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

  end function

  function dvf_row_sum (dvf) result (u)

    type(dist_vf), intent(in) :: dvf
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

  end function dvf_row_sum

  function unpack_dvf_col (dvf, n) result (u)

    type(dist_vf), intent(in) :: dvf
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

  end function unpack_dvf_col

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

  end function

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! VF_DIFF
 !!
 !! Replaces the distributed VF matrix B with the difference B - A.  The
 !! matrices need not have the same sparsity pattern, which requires the
 !! redefinition of the internal structures of B.  However, we do require
 !! that A and B are identically distributed across processes.  We are
 !! careful to preserve the symmetric sparse structure of the difference,
 !! not compressing out any new zeros introduced by exact cancellation.
 !!

  subroutine vf_diff (a, b)

    type(dist_vf), intent(in)    :: a
    type(dist_vf), intent(inout) :: b

    integer :: i, j, n
    integer, allocatable :: ia(:), ja(:)
    real, allocatable :: val(:)
    real :: row(b%npatch_tot)
    logical :: tag(b%npatch_tot)

    ASSERT( a%npatch == b%npatch )

    !! First pass: generate the IA indexing array.
    allocate(ia(b%npatch+1))
    ia(1) = 1
    do j = 1, b%npatch
      !! Tag the nonzeros in the B row.
      tag = .false.
      do i = b%ia(j), b%ia(j+1)-1
        tag(b%ja(i)) = .true.
      end do
      !! Tag the nonzeros in the A row.
      do i = a%ia(j), a%ia(j+1)-1
        tag(a%ja(i)) = .true.
      end do
      ia(j+1) = ia(j) + count(tag)
    end do

    n = ia(b%npatch+1) - 1
    allocate(val(n), ja(n))

    !! Second pass: generate the sparse matrix rows.
    n = 1
    do j = 1, b%npatch
      !! Unpack the compressed B row.
      row = 0.0
      tag = .false.
      do i = b%ia(j), b%ia(j+1)-1
        tag(b%ja(i)) = .true.
        row(b%ja(i)) = b%val(i)
      end do
      !! Subtract the compressed A row.
      do i = a%ia(j), a%ia(j+1)-1
        tag(a%ja(i)) = .true.
        row(a%ja(i)) = row(a%ja(i)) - a%val(i)
      end do
      !! Pack the row.
      do i = 1, size(row)
        if (tag(i)) then
          val(n) = row(i)
          ja(n) = i
          n = n + 1
        end if
      end do
      INSIST( n == ia(j+1) )
    end do

    !! Redefine the internals of B.
    call move_alloc(ia, b%ia)
    call move_alloc(ja, b%ja)
    call move_alloc(val, b%val)

    ASSERT(a%has_ambient .eqv. b%has_ambient)
    if (a%has_ambient .and. b%has_ambient) then
      b%ambient = b%ambient - a%ambient
    end if

  end subroutine vf_diff

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! VF_DIFF_MAX_NORM
 !!
 !! Computes the max norm of the VF difference matrix D as an operator.  If e
 !! is the irradiance discrepancy resulting from D and a given radiosity q,
 !! then |e|_max <= |D|_max |q|_max, where |D|_max is maximum absolute row
 !! sum of D.
 !!
 !! NOTE: I've decide to not include the ambient VF contribution to the error
 !! since I don't know how to handle it in the other norms.
 !!

  subroutine vf_diff_max_norm (d, dnorm)

    type(dist_vf), intent(in) :: d
    real, intent(out) :: dnorm

    integer :: i, j
    real(kind(1.0d0)) :: s, drowsum(d%npatch)

    do j = 1, d%npatch
      s = 0.0 !abs(d%ambient(j))
      do i = d%ia(j), d%ia(j+1)-1
        s = s + abs(d%val(i))
      end do
      drowsum(j) = s
    end do

    dnorm = scl_global_maxval(drowsum)

  end subroutine vf_diff_max_norm

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! VF_DIFF_ONE_NORM
 !!
 !! Computes the L_1 norm of the VF difference matrix D as an operator.  If e
 !! is the irradiance discrepancy resulting from D and a given radiosity q,
 !! then |e|_1 <= |D|_1 |q|_1, where |D|_1 is the sum of the maximum abs value
 !! in each column of D.  Note that the vector norms are properly weighted by
 !! the area of the faces.  This exploits VF reciprocity.
 !!
 !! NOTE: I'm not sure how to handle the ambient VF part because we don't have
 !! that virtual row of the matrix.  So the norm reflects no ambient radiation.
 !!

  subroutine vf_diff_one_norm (d, dnorm)

    type(dist_vf), intent(in) :: d
    real, intent(out) :: dnorm

    integer :: i, j, n
    real :: s2
    real(kind(1.0d0)) :: s1

    s1 = 0.0d0
    do n = 1, d%npatch_tot
      !! Everybody find their max absolute value in column N.
      s2 = 0.0
      do j = 1, d%npatch
        i = loc_in_row(n, d%ja(d%ia(j):d%ia(j+1)-1))
	if (i > 0) then
	  i = i + d%ia(j) - 1
	  s2 = max(s2, abs(d%val(i)))
	end if
      end do
      s1 = s1 + scl_global_maxval(s2)
    end do
    dnorm = s1

  end subroutine vf_diff_one_norm

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! VF_DIFF_TWO_NORM
 !!
 !! Computes the L_2 norm of the VF difference matrix D as an operator.  If e
 !! is the irradiance discrepancy resulting from D and a given radiosity q,
 !! then |e|_2 <= |D|_2 |q|_2, where |D|_2 is the square root of the sum over
 !! j of the dot product of the jth row with the jth column of D.  Note that
 !! the vector norms are properly weighted by the area of the faces.  This
 !! exploits VF reciprocity.
 !!
 !! NOTE: I'm not sure how to handle the ambient VF part because we don't have
 !! that virtual row of the matrix.  So the norm reflects no ambient radiation.
 !!

  subroutine vf_diff_two_norm (d, dnorm)

    type(dist_vf), intent(in) :: d
    real, intent(out) :: dnorm

    integer :: i, j, n, cnt, recvcnt, root, last(scl_size())
    integer :: sidx(d%npatch), sidx_g(d%npatch_tot)
    real    :: svec(d%npatch), svec_g(d%npatch_tot)
    real(kind(1.0d0)) :: s(d%npatch)

    !! The last global row index owned by each process.
    call scl_allgather (d%offset+d%npatch, last)

    ASSERT( last(scl_size()) == d%npatch_tot )

    root = 1
    do n = 1, d%npatch_tot

      !! Update the root process rank: owner of row N.
      do while (n > last(root))
        root = root + 1
      end do

      !! Everybody gather up their nonzeros in column N.
      cnt = 0
      do j = 1, d%npatch
      	i = loc_in_row (n, d%ja(d%ia(j):d%ia(j+1)-1))
	if (i > 0) then
	  i = i + d%ia(j) - 1
          cnt = cnt + 1
          svec(cnt) = d%val(i)
          sidx(cnt) = j + d%offset
	end if
      end do

      !! Gather the compressed sparse column.  The row indices SIDX_G should
      !! be identical to the column indices of row N because of reciprocity.
      recvcnt = scl_global_sum(cnt)
      call scl_gather (svec(:cnt), svec_g(:recvcnt), root)
      call scl_gather (sidx(:cnt), sidx_g(:recvcnt), root)

      !! Form the dot product of row N and column N.
      if (scl_rank() == root) then
        j = n - d%offset  ! local row index
        INSIST( recvcnt == d%ia(j+1)-d%ia(j) )
        INSIST( all(sidx_g(:recvcnt) == d%ja(d%ia(j):d%ia(j+1)-1)) )
        s(j) = dot_product(d%val(d%ia(j):d%ia(j+1)-1), svec_g(:recvcnt))
      end if

    end do

    dnorm = sqrt(scl_global_sum(s))

  end subroutine vf_diff_two_norm

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! VF_DIFF_MAX
 !!
 !! Finds the element of the VF difference matrix D of maximum absolute value
 !! and its row and column location.  The ambient VF difference is ignored
 !! here for consistency with the other norms.
 !!

  subroutine vf_diff_max (d, val, rloc, cloc)

    type(dist_vf), intent(in) :: d
    real, intent(out) :: val
    integer, intent(out) :: rloc, cloc

    integer :: i, i1, i2, n
    integer :: rlocs(scl_size()), clocs(scl_size())
    real :: vals(scl_size())

    val = 0.0
    rloc = 0
    cloc = 0

    if (size(d%val) > 0) then
      n = maxloc(abs(d%val), dim=1)
      val = abs(d%val(n))
      cloc = d%ja(n)

      !! Find the row number using a binary search on the ordered array IA:
      !! I1, I2 satisfy IA(I1)<=N<IA(I2); seach terminates when I2-I1==1
      !! with I1 the row number.
      i1 = 1; i2 = size(d%ia)
      do while (i2-i1 > 1)
        i = (i1 + i2)/2
        if (n >= d%ia(i)) then
          i1 = i
        else
          i2 = i
        end if
      end do
      rloc = i1 + d%offset
    end if

    call scl_allgather (val,  vals)
    call scl_allgather (rloc, rlocs)
    call scl_allgather (cloc, clocs)

    n = maxloc(vals, dim=1)

    val  = vals(n)
    rloc = rlocs(n)
    cloc = clocs(n)

  end subroutine vf_diff_max

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
