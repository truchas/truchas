!!
!! VF_MATRIX_CLASS
!!
!! A common interface to distributed face-based view factor matrices used by the
!! enclosure radiosity solver.
!!
!! David Neill-Asanza <dhna@lanl.gov>
!! 18 July 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!   This module defines the abstract base class VF_MATRIX which defines a
!!   generic face-base distributed view factor matrix and provides a common
!!   interface for operating on concrete instances of this class.  Extensions
!!   to this class must accept and return face-based arguments, regardless of
!!   how the enclosure patches were defined for the view factor calculation.
!!   In particular, if the VF_MATRIX includes ambient view factor data, then
!!   THIS%AMB_VF must be a face-based array.
!!
!!   Application code is expected to use polymorphic variables of this type and
!!   not work directly with extensions of the type.  The base type defines the
!!   following type bound procedures:
!!
!!   INIT(THIS, FILE) initializes the matrix data structures from FILE, a
!!   radiation enclosure file, and defines a distribution of matrix rows among
!!   the ranks.  Since each enclosure face corresponds to a patch (a row), the
!!   distribution of rows imposes a distribution on the faces.  The enclosure
!!   face distribution can be retrieved with the PARTITION_ER_FACES subroutine.
!!   The actual view factor data is not read in INIT, and must be loaded with
!!   the LOAD_VIEW_FACTORS subroutine.
!!
!!   PARTITION_ER_FACES(THIS, COLOR) returns the distribution of faces among
!!   ranks that corresponds to the expected distribution of matrix rows defined
!!   by THIS%INIT.  COLOR(:) is a rank-1 integer array that maps faces to the
!!   rank that owns them (i.e. their color).
!!
!!   LOAD_VIEW_FACTORS(THIS, FILE, ENCL) reads the view factor data from FILE
!!   and distributes the matrix rows according to the distribution defined in
!!   THIS%INIT.  ENCL is a RAD_ENCL type that stores only the subset of the mesh
!!   corresponding to the faces owned by this rank, as specified by the output
!!   of PARTITION_ER_FACES.  Therefore, the enclosure mesh must be distributed
!!   before calling LOAD_VIEW_FACTORS.
!!
!!   PHI_X(THIS, LHS, X) computes the global product LHS=PHI*X, where PHI is
!!   this view factor matrix, and LHS and X are rank-1 real arrays local to this
!!   process.  The elements of X must correspond to the faces owned by this
!!   rank, and must be in the same relative ordering defined by the output of
!!   PARTITION_ER_FACES.  Likewise, the elements of LHS correspond to the faces
!!   owned by this rank, and preserve the ordering of X.
!!


#include "f90_assert.fpp"

module vf_matrix_class

  use kinds, only: r8
  use rad_encl_type
  use rad_encl_file_type
  use parallel_communication
  implicit none
  private

  type, abstract, public :: vf_matrix
    integer :: nface      ! number of faces on this process
    integer :: nface_tot  ! total number of faces (number of columns)
    integer, allocatable :: ia(:), ja(:)
    real, allocatable :: vf(:), amb_vf(:)
    logical :: has_ambient

  contains
    procedure(init), deferred                :: init
    procedure(partition_ER_faces), deferred  :: partition_ER_faces
    procedure(load_view_factors), deferred   :: load_view_factors
    procedure(phi_x), deferred               :: phi_x
    procedure :: distribute_vf_rows
  end type

  abstract interface

    subroutine init (this, file)
      import :: vf_matrix, rad_encl_file
      class(vf_matrix), intent(out) :: this
      type(rad_encl_file), intent(in) :: file
    end subroutine

    subroutine partition_ER_faces (this, color)
      import :: vf_matrix
      class(vf_matrix), intent(inout) :: this
      integer, allocatable, intent(out) :: color(:)
    end subroutine partition_ER_faces

    subroutine load_view_factors (this, file, encl)
      import :: vf_matrix, rad_encl_file, rad_encl
      class(vf_matrix), intent(inout) :: this
      type(rad_encl_file), intent(in) :: file
      type(rad_encl), intent(in) :: encl
    end subroutine load_view_factors

    subroutine phi_x (this, lhs, x)
      import :: vf_matrix, r8
      class(vf_matrix), intent(in) :: this
      real(r8), intent(out) :: lhs(:)
      real(r8), intent(in) :: x(:)
    end subroutine phi_x

  end interface

contains

  !! An auxiliary routine for LOAD_VIEW_FACTORS that reads view factor data from
  !! FILE and distributes it in a block pattern.  NROWS is the number of rows
  !! owned by this process, and NROWS_TOT is the total number of rows (and
  !! columns) of the matrix.
  subroutine distribute_vf_rows (this, file, nrows, nrows_tot)

    use,intrinsic :: iso_fortran_env, only: i8 => int64

    class(vf_matrix), intent(inout) :: this
    type(rad_encl_file), intent(in) :: file
    integer, intent(in) :: nrows      ! number of rows on this process
    integer, intent(in) :: nrows_tot  ! total number of rows (and columns)

    integer(i8) :: start
    integer :: j, n
    integer :: vf_bsize(nPE), lengths(nPE), idum0(0)
    real :: rdum0(0)
    integer, allocatable :: ibuf(:)
    real, allocatable :: rbuf(:)
    logical :: has_ambient

    !! Read and distribute the ambient view factors.
    if (is_IOP) has_ambient = file%has_ambient()
    call broadcast(has_ambient)
    this%has_ambient = has_ambient

    if (has_ambient) then
      allocate(rbuf(merge(nrows_tot, 0, is_IOP)))
      if (is_IOP) call file%get_ambient(rbuf)
      allocate(this%amb_vf(nrows))
      call distribute (this%amb_vf, rbuf)
      deallocate(rbuf)
    end if

    !! Read and distribute the VF matrix nonzero row counts.
    allocate(ibuf(merge(nrows_tot, 0, is_IOP)))
    if (is_IOP) call file%get_vf_rowcount(ibuf)
    allocate(this%ia(nrows+1))
    call distribute (this%ia(2:), ibuf)
    deallocate(ibuf)

    !! Convert the row counts into the local IA indexing array.
    this%ia(1) = 1
    do j = 2, size(this%ia)
      this%ia(j) = this%ia(j) + this%ia(j-1)
    end do

    !! Determine the sizes of the distributed VF matrix.
    n = this%ia(nrows+1) - this%ia(1)
    allocate(this%vf(n), this%ja(n))
    call collate (vf_bsize, n)

    !! Read the VF matrix in process-sized blocks, sending them to the owning
    !! processes as we go.  PGSLib only provides the distribute collective to
    !! do this, and so we distribute with 0-sized destination arrays for all
    !! processes but the receiving one.  We also use the optional LENGTHS
    !! argument to distribute (only referenced on the IO process) that gives
    !! the number of items to distribute to each process, resulting in the
    !! actual sizes of the array arguments being ignored except to check that
    !! they are sufficiently large.  This allows us to simplify the calls.
    !! The code is structured for future use of MPI_Isend/Irecv, abandoning
    !! the usual 'single code path' pattern.

    if (is_IOP) then
      call file%get_vf_rows(this%vf, this%ja, start=1_i8)
      if (nPE > 1) then
        n = maxval(vf_bsize(2:))
        allocate(ibuf(n), rbuf(n))
        lengths = 0
        start = 1 + vf_bsize(1)
        do n = 2, nPE
          call file%get_vf_rows(rbuf(:vf_bsize(n)), ibuf(:vf_bsize(n)), start)
          lengths(n) = vf_bsize(n)
          call distribute(rdum0, rbuf, lengths)
          call distribute(idum0, ibuf, lengths)
          lengths(n) = 0
          start = start + vf_bsize(n)
        end do
        deallocate(ibuf, rbuf)
      end if
    else  ! everybody else just participates in the distribute calls.
      do n = 2, nPE
        if (n == this_PE) then
          call distribute(this%vf, rdum0, lengths)
          call distribute(this%ja, idum0, lengths)
        else
          call distribute(rdum0, rdum0, lengths)
          call distribute(idum0, idum0, lengths)
        end if
      end do
    end if

  end subroutine distribute_vf_rows

end module vf_matrix_class
