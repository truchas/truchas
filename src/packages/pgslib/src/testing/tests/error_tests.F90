module ERROR_TESTS
  use pgslib_module
  implicit none 
  ! Routines to check for significant differences between real numbers
  PRIVATE
  public :: WITHIN_TOLERANCE

  real(pgslib_int_type),    save  :: integer_tol
  real(pgslib_single_type), save  :: single_tol
  real(pgslib_double_type), save  :: double_tol

  logical, save :: tol_initialized = .false.
  
  INTERFACE WITHIN_TOLERANCE
     MODULE PROCEDURE single_within_tol
     MODULE PROCEDURE double_within_tol
     MODULE PROCEDURE integer_within_tol
  END INTERFACE

CONTAINS
  subroutine Initialize_Tolerances()
    ! Initialize the tolerances if necessary
    integer(pgslib_int_type) :: npe

    if (.not. tol_initialized) then
       npe = PGSLib_Inquire_nPE()
       single_tol = npe*EPSILON(single_tol)
       double_tol = npe*EPSILON(double_tol)
       integer_tol = 0
       tol_initialized = .TRUE.
    end if

    return
  end subroutine Initialize_Tolerances

  function GET_SCOPE(SCOPE)
    ! If SCOPE is present, then return scope,
    ! else return the default scope, which is GLOBAL
    type(PGSLib_SCOPE), intent(IN), OPTIONAL :: SCOPE
    type(PGSLib_SCOPE) :: GET_SCOPE

    if (PRESENT(SCOPE)) then
       GET_SCOPE = SCOPE
    else
       GET_SCOPE = PGSLib_Global
    end if
    return
  end function GET_SCOPE

  function single_within_tol(DESIRED, FOUND, SCOPE)
    logical(PGSLib_Log_Type) :: single_within_tol
    real(pgslib_single_type), intent(IN) :: DESIRED
    real(pgslib_single_type), intent(IN) :: FOUND
    type(PGSLib_SCOPE), intent(IN), OPTIONAL :: SCOPE
    
    single_within_tol = pgslib_global_any((ABS(DESIRED - FOUND) < single_tol), GET_SCOPE(SCOPE))
    return
  end function single_within_tol

  function double_within_tol(DESIRED, FOUND, SCOPE)
    logical(PGSLib_Log_Type) :: double_within_tol
    real(pgslib_double_type), intent(IN) :: DESIRED
    real(pgslib_double_type), intent(IN) :: FOUND
    type(PGSLib_SCOPE), intent(IN), OPTIONAL :: SCOPE

    double_within_tol = pgslib_global_any((ABS(DESIRED - FOUND) < double_tol), GET_SCOPE(SCOPE))
    return
  end function double_within_tol

  function integer_within_tol(DESIRED, FOUND, SCOPE)
    logical(PGSLib_Log_Type) :: integer_within_tol
    integer(pgslib_int_type), intent(IN) :: DESIRED
    integer(pgslib_int_type), intent(IN) :: FOUND
    type(PGSLib_SCOPE), intent(IN), OPTIONAL :: SCOPE

    integer_within_tol = pgslib_global_any((ABS(DESIRED - FOUND) < integer_tol), GET_SCOPE(SCOPE))
    return
  end function integer_within_tol

end module ERROR_TESTS
