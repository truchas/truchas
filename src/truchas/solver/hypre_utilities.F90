!!
!! HYPRE_UTILITIES
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module hypre_utilities

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use fhypre
  implicit none
  private

  public :: copy_to_ijmatrix

  interface copy_to_ijmatrix
    module procedure csr_to_ijmatrix
  end interface

contains

  !! Copy a CSR_MATRIX object SRC to an equivalent HYPRE_IJMatrix object.
  !! The HYPRE matrix is created if it does not exist. Otherwise the elements
  !! of the existing HYPRE matrix are overwritten with the values from SRC.
  !! In the latter case the sparsity pattern of the two matrices must be
  !! identical. This is a serial procedure and the resulting HYPRE matrix is
  !! associated with the MPI_COMM_SELF communicator. On entry MATRIX must
  !! either be a valid HYPRE object or have the value HYPRE_NULL_OBJ.

  subroutine csr_to_ijmatrix(src, matrix)

    use csr_matrix_type
    use mpi, only: MPI_COMM_SELF

    class(csr_matrix), intent(in) :: src
    type(hypre_obj), intent(inout) :: matrix

    integer :: j, ierr
    integer, allocatable :: ncols_offP(:), ncols(:), rows(:)

    call fHYPRE_ClearAllErrors

    rows = [(j, j = 1, src%nrow)]
    ncols = src%graph%xadj(2:) - src%graph%xadj(:src%nrow)

    if (.not.hypre_associated(matrix)) then
      allocate(ncols_offP(src%nrow), source=0)
      call fHYPRE_IJMatrixCreate(MPI_COMM_SELF, 1, src%nrow, 1, src%ncol, matrix, ierr)
      call fHYPRE_IJMatrixSetDiagOffdSizes(matrix, ncols, ncols_offP, ierr)
      call fHYPRE_IJMatrixSetMaxOffProcElmts(matrix, 0, ierr) ! probably not necessary
      INSIST(ierr == 0)
    end if

    !! After initialization the HYPRE matrix elements can be set.
    call fHYPRE_IJMatrixInitialize(matrix, ierr)
    INSIST(ierr == 0)

    !! Copy the matrix elements into the HYPRE matrix.  This defines both the
    !! nonzero structure of the matrix and the values of those elements.
    call fHYPRE_IJMatrixSetValues(matrix, src%nrow, ncols, rows, src%graph%adjncy, src%values, ierr)
    INSIST(ierr == 0)

    !! After assembly the HYPRE matrix is ready to use.
    call fHYPRE_IJMatrixAssemble(matrix, ierr)
    INSIST(ierr == 0)

  end subroutine csr_to_ijmatrix

end module hypre_utilities
