MODULE PARALLEL_SCOPE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define data types and parameters to support parallel scope
  !   distinctions, such as global vs. local operations.
  !
  ! Contains: None
  !
  ! Author(s): Robert Ferrell (ferrell@cpca.com)
  !
  !=======================================================================
  use kind_module
  implicit none
  save

  ! Public Module
  public

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! SCOPE definitions
  type PL_SCOPE
     integer(int_kind) :: SCOPE
  end type PL_SCOPE
  ! SCOPE Parameters
  type (PL_SCOPE), parameter :: GLOBAL_SCOPE = PL_SCOPE(1)
  type (PL_SCOPE), parameter :: LOCAL_SCOPE  = PL_SCOPE(2)

  ! Interface definitions
  INTERFACE OPERATOR (.EQ.)
     MODULE PROCEDURE EQ_SCOPE
  END INTERFACE

  INTERFACE OPERATOR (.NE.)
     MODULE PROCEDURE NEQ_SCOPE
  END INTERFACE
CONTAINS

  function EQ_SCOPE(A, B)
    implicit none
    type (PL_SCOPE), intent(IN) :: A, B
    logical (Log_kind)          :: EQ_SCOPE

    EQ_SCOPE = (A%Scope == B%Scope)

    RETURN
  end function EQ_SCOPE
  
  function NEQ_SCOPE(A, B)
    implicit none
    type (PL_SCOPE), intent(IN) :: A, B
    logical (Log_kind)          :: NEQ_SCOPE

    NEQ_SCOPE = (A%Scope /= B%Scope)

    RETURN
  end function NEQ_SCOPE
END MODULE PARALLEL_SCOPE

  
  
