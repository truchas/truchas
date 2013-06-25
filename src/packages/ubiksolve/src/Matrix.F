! -*- Mode: f90 -*-

module Ubik_Matrix_module
  use ubik_types_module
  implicit none

  interface assignment(=)
     module procedure CSC2ELL
  end interface

  private
  public :: assignment(=)

contains

  !*****************************************************************
  !
  ! <para>
  ! Set ELL matrix from CSC matrix.
  ! </para>
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   csc - Compressed Sparse Column matrix
  ! 
  !  Output:
  !   ell - ELLPACK-ITPACK matrix
  !
  !*****************************************************************
  subroutine CSC2ELL (ell, csc)
    use Ubik_CSC_module
    use Ubik_ELL_module

    ! <arguments>
    type(CSC_matrix), intent(in) :: csc
    type(ELL_matrix), intent(inout) :: ell
    ! </arguments>

    ! Local:
    integer(int_type) :: i, ii, j, jj, status
    integer(int_type) :: nnzero
    integer(int_type) :: maxnz, nrows
    integer(int_type), allocatable, dimension(:) :: colind, jpos, nnzrow

    nnzero = SIZE(csc%values)
    nrows = SIZE(csc%rowind)

    allocate (nnzrow(nrows), colind(nnzero), stat=status)
    if (status /= 0) then
       call set_status_alloc_failure(ell)
       return
    end if

    ! need to find max number of nonzero elements in any row for allocation of ELL
    nnzrow = 0
    ii = 1
    do i = 1, nnzero
       if (i >= csc%colptr(ii+1)) then
          ii = ii + 1
       end if
       colind(i) = ii
       nnzrow(ii) = nnzrow(ii) + 1
    end do
    maxnz = MaxVal(nnzrow)

    deallocate(nnzrow, stat=status)
    if (status /= 0) then
       call set_status_alloc_failure(ell)
       return
    end if

    allocate(jpos(nrows), stat=status)
    if (status /= 0) then
       call set_status_alloc_failure(ell)
       return
    end if

    ! create will resize if necessary
    call create (ell, nrows, maxnz)
    if (alloc_failure(ell)) then
       return
    end if

    jpos = 0
    ell%values = 0.0
    ell%map = 0
    do i = 1, nnzero
       ii = csc%rowind(i)
       jj = colind(i)
       jpos(ii) = jpos(ii) + 1
       j = jpos(ii)
       ell%values(ii,j) = csc%values(i)
       ell%map(ii,j) = jj
    end do

    deallocate (jpos, colind, stat=status)
    if (status /= 0) then
       call set_status_alloc_failure(ell)
    end if

    return
  end subroutine CSC2ELL

end module Ubik_Matrix_module
