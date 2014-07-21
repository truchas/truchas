MODULE GATHER_MODULE
  !=======================================================================
  ! PURPOSE - 
  !   Routines for gather.  These are used by both 
  !   Element<->Element and Element<->Node.
  ! 
  !=======================================================================
  use truchas_logging_services
  use gs_info_module
  use kinds, only: r8
  use mesh_module,  only: MESH_CONNECTIVITY
  use pgslib_module,only: PGSLib_GS_Trace,      &
                          PGSLib_Size_Of_Dup,   &
                          PGSLib_Size_Of_Sup,   &
                          PGSLib_Dup_Index,     &
                          PGSLib_Gather_Buffer

  implicit none
  private
  save

  ! Public
  public :: GATHER
  public :: Gather_BoundaryData

  ! Interface blocks
  INTERFACE GATHER
    !=======================================================================
    ! PURPOSE - 
    !   Gather either cell-centered data or vertex centered data
    !   into a cell-centered vector of length ShortDim
    !   
    !    Input: 
    !          Mesh - Mesh connectivity structure
    !          Src  - Cell or vertex -centered quantity 
    !   Output: 
    !          Dest - Cell-centered vector of length ShortDim containing
    !                 the value of Src for all neighbors
    !          
    !=======================================================================
     MODULE PROCEDURE Gather_INT
     MODULE PROCEDURE Gather_SINGLE
     MODULE PROCEDURE Gather_DOUBLE
     MODULE PROCEDURE Gather_LOG
     MODULE PROCEDURE Gather_V_V_INT
     MODULE PROCEDURE Gather_V_V_SINGLE
     MODULE PROCEDURE Gather_V_V_DOUBLE
     MODULE PROCEDURE Gather_V_V_LOG
  END INTERFACE
  
  INTERFACE Gather_BoundaryData
    MODULE PROCEDURE Gather_BoundaryData_Int
    MODULE PROCEDURE Gather_BoundaryData_Double
  END INTERFACE

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

    !=======================================================================
    ! PURPOSE - 
    !   Gather either cell-centered data or vertex centered integer data
    !   into an integer cell-centered vector of length ShortDim
    !=======================================================================

#define _ROUTINE_NAME_ GATHER_INT
#define _DATA_TYPE_    integer
#define _OP_ID_        0

#include "gather_parallel_v_s_include.fpp"

    !=======================================================================
    ! PURPOSE - 
    !   Gather either cell-centered data or vertex centered single data
    !   into a single cell-centered vector of length ShortDim
    !=======================================================================

#define _ROUTINE_NAME_ GATHER_SINGLE
#define _DATA_TYPE_    real
#define _OP_ID_        0.0

#include "gather_parallel_v_s_include.fpp"
  
    !=======================================================================
    ! PURPOSE - 
    !   Gather either cell-centered data or vertex centered double data
    !   into a double cell-centered vector of length ShortDim
    !=======================================================================

#define _ROUTINE_NAME_ GATHER_DOUBLE
#define _DATA_TYPE_    real(r8)
#define _OP_ID_        0.0_r8

#include "gather_parallel_v_s_include.fpp"
  
    !=======================================================================
    ! PURPOSE - 
    !   Gather either cell-centered data or vertex centered logical data
    !   into a logical cell-centered vector of length ShortDim
    !=======================================================================

#define _ROUTINE_NAME_ GATHER_LOG
#define _DATA_TYPE_    logical
#define _OP_ID_        .FALSE.

#include "gather_parallel_v_s_include.fpp"
  
    !=======================================================================
    ! PURPOSE - 
    !   Gather either cell-centered data or vertex centered integer data
    !   into an integer cell-centered vector of length ShortDim
    !=======================================================================

#define _ROUTINE_NAME_ GATHER_V_V_INT
#define _DATA_TYPE_    integer
#define _OP_ID_        0

#include "gather_parallel_v_v_include.fpp"

    !=======================================================================
    ! PURPOSE - 
    !   Gather either cell-centered data or vertex centered single data
    !   into a single cell-centered vector of length ShortDim
    !=======================================================================

#define _ROUTINE_NAME_ GATHER_V_V_SINGLE
#define _DATA_TYPE_    real
#define _OP_ID_        0.0

#include "gather_parallel_v_v_include.fpp"

    !=======================================================================
    ! PURPOSE - 
    !   Gather either cell-centered data or vertex centered double data
    !   into a double cell-centered vector of length ShortDim
    !=======================================================================

#define _ROUTINE_NAME_ GATHER_V_V_DOUBLE
#define _DATA_TYPE_    real(r8)
#define _OP_ID_        0.0_r8

#include "gather_parallel_v_v_include.fpp"
  
    !=======================================================================
    ! PURPOSE - 
    !   Gather either cell-centered data or vertex centered logical data
    !   into a logical cell-centered vector of length ShortDim
    !=======================================================================

#define _ROUTINE_NAME_ GATHER_V_V_LOG
#define _DATA_TYPE_    logical
#define _OP_ID_        .FALSE.

#include "gather_parallel_v_v_include.fpp"
  
  
  !-----------------------------------------------------------------------------

  subroutine Gather_BoundaryData_Int (BOUNDARY, SOURCE)
    !---------------------------------------------------------------------------
    ! RCF June 24, 1998
    ! Moved this from preconditioner_module.F90.
    ! May/Will modify/generalize/streamline as we clean up solver/precond code
    !---------------------------------------------------------------------------    

    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !   update the X values from neighboring domains
    !
    !   The functionality of this bit of code is included in
    !   EE_GATHER, but EE_GATHER does more.  This bit gets the
    !   neighboring domain x values and puts them in x%aux1.
    !   EE_GATHER continues by filling them into an ncells x nfc
    !   temporary.  Use EE_GATHER if you want to do vector operations
    !   with the temporary.  Use this if you have to do indirect
    !   addressing yourself.
    !
    !   This code should probably reside in PGSLIB somewhere, but for
    !   now it stays here.  It is only called locally from other
    !   routines in preconditioner_module, and is not public.
    !---------------------------------------------------------------------------

    ! arguments
    integer, INTENT(IN) :: SOURCE(:)
    integer, POINTER :: BOUNDARY(:)
    

    ! local variables
    integer :: status
    integer, POINTER, Dimension(:) :: Duplicate_Data

    !---------------------------------------------------------------------------

    ! if BOUNDARY doesn't point at anything, get some space
    ! if BOUNDARY is already allocated, don't do anything.
    If (.Not. Associated(BOUNDARY)) Then
       ALLOCATE(Duplicate_Data(PGSLib_Size_Of_Dup(EE_Trace)))
       Duplicate_Data = SOURCE(PGSLib_Dup_Index(EE_Trace))
       Allocate(BOUNDARY(PGSLib_Size_OF_Sup(EE_Trace)), STAT=status)
       Call TLS_fatal_if_any ((status /= 0), 'Gather_BoundaryData: error allocating BOUNDARY')

       ! the communication, takes data from Duplicate buffer, puts data into BOUNDARY
       BOUNDARY = PGSLib_gather_Buffer(Duplicate_Data, EE_Trace)
       ! Clean up
       DEALLOCATE(Duplicate_Data)
    End If

  End subroutine Gather_BoundaryData_Int


  subroutine Gather_BoundaryData_Double (BOUNDARY, SOURCE)
    !---------------------------------------------------------------------------
    ! RCF June 24, 1998
    ! Moved this from preconditioner_module.F90.
    ! May/Will modify/generalize/streamline as we clean up solver/precond code
    !---------------------------------------------------------------------------    

    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !   update the X values from neighboring domains
    !
    !   The functionality of this bit of code is included in
    !   EE_GATHER, but EE_GATHER does more.  This bit gets the
    !   neighboring domain x values and puts them in x%aux1.
    !   EE_GATHER continues by filling them into an ncells x nfc
    !   temporary.  Use EE_GATHER if you want to do vector operations
    !   with the temporary.  Use this if you have to do indirect
    !   addressing yourself.
    !
    !   This code should probably reside in PGSLIB somewhere, but for
    !   now it stays here.  It is only called locally from other
    !   routines in preconditioner_module, and is not public.
    !---------------------------------------------------------------------------

    ! arguments
    real(r8), INTENT(IN) :: SOURCE(:)
    real(r8), POINTER :: BOUNDARY(:)
    

    ! local variables
    integer :: status
    real(r8), POINTER, Dimension(:) :: Duplicate_Data

    !---------------------------------------------------------------------------

    ! if BOUNDARY doesn't point at anything, get some space
    ! if BOUNDARY is already allocated, don't do anything.
    If (.Not. Associated(BOUNDARY)) Then
       ALLOCATE(Duplicate_Data(PGSLib_Size_Of_Dup(EE_Trace)))
       Duplicate_Data = SOURCE(PGSLib_Dup_Index(EE_Trace))
       Allocate(BOUNDARY(PGSLib_Size_OF_Sup(EE_Trace)), STAT=status)
       Call TLS_fatal_if_any ((status /= 0), 'Gather_BoundaryData: error allocating BOUNDARY')

       ! the communication, takes data from Duplicate buffer, puts data into BOUNDARY
       BOUNDARY = PGSLib_gather_Buffer(Duplicate_Data, EE_Trace)
       ! Clean up
       DEALLOCATE(Duplicate_Data)
    End If

  End subroutine Gather_BoundaryData_Double

END MODULE GATHER_MODULE
