! -*- Mode: f90 -*-

!*********************************************************************
!
! $Id: CSC.F,v 1.11 2002/05/12 18:58:48 turner Exp $
! 
!*********************************************************************
module Ubik_CSC_module
  use ubik_types_module
  implicit none

  ! Type definition.

  type CSC_matrix
     real(real_type), dimension(:), pointer :: values => NULL()
     integer(int_type), dimension(:), pointer :: rowind => NULL()
     integer(int_type), dimension(:), pointer :: colptr => NULL()
     integer(int_type) :: status
  end type CSC_matrix

  ! Operators.

  interface assignment(=)
     module procedure CSC_eq_int_type
     module procedure CSC_eq_real_type
     module procedure CSC_eq_CSC
  end interface

  ! Statements.

  interface create
     module procedure create_and_size
     module procedure create_and_size_and_init_int
     module procedure create_and_size_and_init_float
  end interface

  interface destroy
     module procedure destroy_CSC
  end interface

  ! Subroutines.

  interface set_status_alloc_failure
     module procedure set_status_alloc_failure_CSC
  end interface

  interface alloc_failure
     module procedure alloc_failure_CSC
  end interface

  interface Load
     module procedure LoadMatrix_CSC
     module procedure LoadSystem_CSC
  end interface

  private
  public :: &
       CSC_matrix, &
       assignment(=), &
       create, &
       destroy, &
       set_status_alloc_failure, &
       alloc_failure, &
       Load

contains

  !*****************************************************************
  !
  ! status setting and accessing
  !
  !*****************************************************************
  subroutine set_status_alloc_failure_CSC(a)
    use ubik_constants_module

    type(CSC_matrix), intent(inout) :: a

    a%status = alloc_failure__
    return
  end subroutine set_status_alloc_failure_CSC

  function alloc_failure_CSC(a)
    use ubik_constants_module

    logical(logical_type) :: alloc_failure_CSC
    type(CSC_matrix), intent(in) :: a

    alloc_failure_CSC = (a%status == alloc_failure__)
    return
  end function alloc_failure_CSC

  !*****************************************************************
  !
  ! Set CSC matrix to scalar.
  ! 
  !*****************************************************************
  subroutine CSC_eq_int_type (a, i)
    integer(int_type), intent(in) :: i
    type(CSC_matrix), intent(inout) :: a

    a%values = DBLE(i)

    return
  end subroutine CSC_eq_int_type
  subroutine CSC_eq_real_type (a, s)
    real(real_type), intent(in) :: s
    type(CSC_matrix), intent(inout) :: a

    a%values = s

    return
  end subroutine CSC_eq_real_type
  subroutine CSC_eq_CSC (a, b)
    type(CSC_matrix), intent(in) :: b
    type(CSC_matrix), intent(inout) :: a

    a%values = b%values
    a%rowind = b%rowind
    a%colptr = b%colptr

    return
  end subroutine CSC_eq_CSC

  !*******************************************************************
  !
  ! constructors and destructors
  !
  !**********************************************************************

  subroutine destroy_CSC (a)

    ! Args:
    type(CSC_matrix), intent(inout) :: a

    ! Local:
    integer(int_type) :: status

    if (associated(a%values)) then
       deallocate(a%values, stat=status)
       if (alloc_failure(a)) then
          return
       end if
    end if
    if (associated(a%rowind)) then
       deallocate(a%rowind, stat=status)
       if (alloc_failure(a)) then
          return
       end if
    end if
    if (associated(a%colptr)) then
       deallocate(a%colptr, stat=status)
       if (alloc_failure(a)) then
          return
       end if
    end if

    return
  end subroutine destroy_CSC

  subroutine create_and_size (a, nnzero, ncol)

    ! Args:
    integer(int_type), intent(in) :: nnzero, ncol
    type(CSC_matrix), intent(inout) :: a

    ! Local:
    integer(int_type) :: status

    ! if resizing, need to deallocate first - note that must check
    ! association status first, since SIZE requires a status of
    ! associated
    if (associated(a%values) .and. &
        associated(a%rowind) .and. &
        associated(a%colptr)) then
       if ((size(a%values) /= nnzero) .or. &
           (size(a%rowind) /= nnzero) .or. &
           (size(a%colptr) /= ncol+1)) then
          call destroy (a)
          if (alloc_failure(a)) then
             return
          end if
       end if
    end if

    ! allocate
    if (.not.associated(a%values) .and. &
        .not.associated(a%rowind) .and. &
        .not.associated(a%colptr)) then
       allocate (a%values(nnzero), a%rowind(nnzero), a%colptr(ncol+1), stat=status)
       if (status /= 0) then
          call set_status_alloc_failure(a)
       end if
    end if

    return
  end subroutine create_and_size

  subroutine create_and_size_and_init_int (a, nnzero, ncol, i)

    integer(int_type), intent(in) :: nnzero, ncol, i
    type(CSC_matrix), intent(inout) :: a

    call create (a, nnzero, ncol)
    if (alloc_failure(a)) then
       return
    end if

    a = i

    return
  end subroutine create_and_size_and_init_int

  subroutine create_and_size_and_init_float (a, nnzero, ncol, x)

    integer(int_type), intent(in) :: nnzero, ncol
    real(real_type) :: x
    type(CSC_matrix), intent(inout) :: a

    call create (a, nnzero, ncol)
    if (alloc_failure(a)) then
       return
    end if

    a = x

    return
  end subroutine create_and_size_and_init_float

  !*********************************************************************
  !
  ! Subroutine to load the nonzero elements of a coefficient matrix
  ! written in Harwell-Boeing format into a CSC matrix.
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   lu - logical unit from which to read
  ! 
  !  Output:
  !   a - coefficient matrix
  !
  !*********************************************************************
  subroutine LoadMatrix_CSC (lu, a)

    ! Args:
    integer(int_type), intent(in) :: lu
    type(CSC_matrix), intent(inout) :: a

    ! Local:
    character :: title*72, key*8, mxtype*3, ptrfmt*16, indfmt*16, valfmt*20, rhsfmt*20, rhstype*3
    integer(int_type) :: i
    integer(int_type) :: totcrd, ptrcrd, indcrd, valcrd, rhscrd, nrow, ncol, &
         nnzero, neltvl, nrhs, nrhsix

    ! header
    read (lu, 1) title, key
    read (lu, 2) totcrd, ptrcrd, indcrd, valcrd, rhscrd
    read (lu, 3) mxtype, nrow, ncol, nnzero, neltvl
    read (lu, 4) ptrfmt, indfmt, valfmt, rhsfmt
    read (lu, 5) rhstype, nrhs, nrhsix
1   format (a72, a8)
2   format (5i14)
3   format (a3, 11x, 4i14)
4   format (2a16, 2a20)
5   format (a3, 11x, 2i14)

    call create (a, nnzero, ncol)

    ! coeff structure
    read (lu, ptrfmt) (a%colptr(i), i = 1, ncol+1)
    read (lu, indfmt) (a%rowind(i), i = 1, nnzero)

    ! coeff values
    if (valcrd > 0) then
       read (lu, valfmt) (a%values(i), i = 1, nnzero)
    end if

    return
  end subroutine LoadMatrix_CSC

  !*********************************************************************
  !
  ! Subroutine to load an entire linear system.
  ! 
  ! <PARAMETER LIST>
  ! 
  !  Input:
  !   lu - logical unit from which to read
  !   a - coefficient matrix
  !   x - vector of unknowns
  !   b - vector containing right-hand side
  !
  !*********************************************************************
  subroutine LoadSystem_CSC (lu, a, x, b)

    ! Args:
    integer(int_type), intent(in) :: lu
    type(CSC_matrix), intent(inout) :: a
    real(real_type), intent(out), dimension(:) :: x, b

    ! Local:
    character :: title*72, key*8, mxtype*3, ptrfmt*16, indfmt*16, valfmt*20, rhsfmt*20, rhstype*3
    integer(int_type) :: i
    integer(int_type) :: totcrd, ptrcrd, indcrd, valcrd, rhscrd, nrow, ncol, &
         nnzero, neltvl, nrhs, nrhsix

    ! header
    read (lu, 1) title, key
    read (lu, 2) totcrd, ptrcrd, indcrd, valcrd, rhscrd
    read (lu, 3) mxtype, nrow, ncol, nnzero, neltvl
    read (lu, 4) ptrfmt, indfmt, valfmt, rhsfmt
    read (lu, 5) rhstype, nrhs, nrhsix
1   format (a72, a8)
2   format (5i14)
3   format (a3, 11x, 4i14)
4   format (2a16, 2a20)
5   format (a3, 11x, 2i14)

    call create (a, nnzero, ncol)

    ! coeff structure
    read (lu, ptrfmt) (a%colptr(i), i = 1, ncol+1)
    read (lu, indfmt) (a%rowind(i), i = 1, nnzero)

    ! coeff values
    if (valcrd > 0) then
       read (lu, valfmt) (a%values(i), i = 1, nnzero)
    end if

    ! unknown and source vectors
    if (SIZE(x) /= nrow) then
       stop
!       deallocate(x)
!       allocate(x(nrow))
    end if
    x = 1.0d0

    if (SIZE(b) /= nrow) then
       stop
!       deallocate(b)
!       allocate(b(nrow))
    end if
    read (lu, rhsfmt) (b(i), i = 1, nrow)

    return
  end subroutine LoadSystem_CSC

end module Ubik_CSC_module
