MODULE PGSLib_Type_MODULE
  IMPLICIT NONE
  SAVE
  PRIVATE
  PUBLIC :: PGSLib_Int_Type, PGSLib_Real_Type, PGSLib_Single_Type, PGSLib_Double_Type
  PUBLIC :: PGSLib_Log_Type, PGSLib_Char_Type
  PUBLIC :: PGSLib_TRUE, PGSLib_FALSE
  PUBLIC :: PGSLib_Timer, PGSLib_Timer_Slot
  PUBLIC :: PGSLib_Instrument_T, PGSLib_Instrument_Slot
  PUBLIC :: PGSLib_PEInfo_Struct, PGSLib_GID, PGSLib_GS_Trace
  PUBLIC :: PGSLib_Access_Table
  PUBLIC :: PGSLib_SCOPE, OPERATOR(==)
  PUBLIC :: PGSLib_Size_Of_Dup, PGSLib_Size_Of_Sup
  PUBLIC :: PGSLib_Dup_Index
  PUBLIC :: PGSLib_Sup_Global_PEs, PGSLib_Sup_Global_Indices

  PUBLIC :: PGSLib_Error_Type, PGSLib_Error_Test, PGSLib_Error_Set
  ! This module sets the types used in the PGSLib library.
  ! It also defines the derived types.
  ! It also has operators for a few of the derived types
  ! It also provides accessor functions to the components

  INTERFACE PGSLib_Dup_Index
     MODULE PROCEDURE PGSLib_Dup_Index_S
     MODULE PROCEDURE PGSLib_Dup_Index_V
  END INTERFACE
  

  ! $Id: pgslib_type_module.F,v 1.1.1.1 2000/10/11 22:44:32 ferrell Exp $

!!!!!!!!!! INTRINSIC DATA TYPES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Most routines support three data types.
  !          PGSLib_INT_TYPE    (typically the default integer on each system)
  !          PGSLib_REAL_TYPE   (typically the default single precision) (Archaic)
  !          PGSLib_Single_TYPE (typically the default single precision) (Modern)
  !          PGSLib_DOUBLE_TYPE (typically the default double precision)
  ! The restriction on the types is that they must correspond to types in C and MPI.
  !
  !          F90 type              C Type        MPI Type
  !
  !          PGSLib_INT_TYPE       int           MPI_INT
  !          PGSLib_REAL_TYPE      float         MPI_FLOAT 
  !          PGSLib_Single_TYPE      float         MPI_FLOAT
  !          PGSLib_DOUBLE_TYPE    double        MPI_DOUBLE
  !
  ! These restrictions are not automatically satisfied.  It is required that
  !          these restrictions be considered on each new system.

  INTEGER, PARAMETER:: PGSLib_INT_TYPE    = KIND(1)  
  INTEGER, PARAMETER:: PGSLib_Single_TYPE = KIND(1.0)
  INTEGER, PARAMETER:: PGSLib_REAL_TYPE   = PGSLib_Single_Type
  INTEGER, PARAMETER:: PGSLib_DOUBLE_TYPE = KIND(1.0d0)
  INTEGER, PARAMETER:: PGSLib_Log_Type    = KIND(.TRUE.)
  INTEGER, PARAMETER:: PGSLib_Char_Type   = KIND('A')


  ! Passing logicals to C is tricky.  One sure way is to convert F90 logicals
  ! to C integers.  We do that, accordinge to this convention
  INTEGER, PARAMETER :: PGSLib_TRUE  = 1
  INTEGER, PARAMETER :: PGSLib_FALSE = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
!!!!!!!!!! DERIVED TYPES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! Type for returning Error Information (not much there yet)
  TYPE PGSLib_Error_Type
     integer (PGSLib_Int_Type) :: Error_Number
  END TYPE PGSLib_Error_Type


!!!!! Type for holding general parallel info
  TYPE PGSLib_PEInfo_Struct
     integer ( PGSLib_INT_TYPE):: nPE
     integer ( PGSLib_INT_TYPE):: thisPE
     integer ( PGSLib_INT_TYPE):: IO_ROOT_PE
  END TYPE PGSLib_PEInfo_Struct




!!!!! Type for Timing
  integer, parameter:: T_STRIN_LEN=1024
  type PGSLib_Timer
     integer                 :: Timer_N
     integer                 :: T_Ref_Counter
     integer, dimension(1:8) :: Start_Clock, Stop_Clock
     real                    :: Start_Time, Stop_Time
     real                    :: Elapsed_Clock
     real                    :: Elapsed_Time
     real                    :: Internal_Time
     logical                 :: Running
     character (LEN=T_STRIN_LEN)    :: Timer_String
  end type PGSLib_Timer

  type PGSLib_Timer_slot
     integer                 :: Slot
  end type PGSLib_Timer_slot

!!!!! Instrumentation Information
  type PGSLib_Instrument_Slot
     integer               :: Slot
  end type PGSLib_Instrument_Slot

  type PGSLib_Instrument_T
     type (PGSLib_Timer)           :: Timer
     type (PGSLib_Instrument_Slot) :: Slot
  end type PGSLib_Instrument_T



!!!!! PGSLib_GID holds the Global Index Descriptor (GID)
!!!!! for each data structure.  This assumes data structures are all 1D.

  TYPE PGSLib_GID
     logical (PGSLib_Log_Type)                       ::       SetupP
     integer (PGSLib_Int_TYPE)                       ::  TotalExtent ! Global Size
     integer (PGSLib_Int_TYPE)                       :: ExtentThisPE ! Local Size
     integer (PGSLib_Int_TYPE), POINTER, DIMENSION(:):: LocalExtents ! Local Sizes
     integer (PGSLib_Int_TYPE), POINTER, DIMENSION(:)::     LowIndex ! Low Global Index
     integer (PGSLib_Int_TYPE), POINTER, DIMENSION(:)::    HighIndex ! High Global Index
  END TYPE  PGSLib_GID

  TYPE PGSLib_Access_Table
     ! The table is actually in C, so this is just a place to hold a pointer
     integer (PGSLib_INT_TYPE),    POINTER, DIMENSION(:) :: C_Table
  END TYPE PGSLib_Access_Table
     

  TYPE PGSLib_GS_TRACE
     logical (PGSLib_Log_Type) :: SetupP

     logical (PGSLib_Log_Type) :: Present_PMask
     integer (PGSLib_INT_TYPE) :: N_Supplement, N_Duplicate

     integer (PGSLib_INT_TYPE),    POINTER, DIMENSION(:):: Supplement_Global_Indices
     integer (PGSLib_INT_TYPE),    POINTER, DIMENSION(:):: Supplement_Global_PEs
     integer (PGSLib_INT_TYPE),    POINTER, DIMENSION(:):: Supplement_I_Data
     real    (PGSLib_REAL_TYPE),   POINTER, DIMENSION(:):: Supplement_R_Data
     real    (PGSLib_DOUBLE_TYPE), POINTER, DIMENSION(:):: Supplement_D_Data
     logical (PGSLib_Log_TYPE),    POINTER, DIMENSION(:):: Supplement_L_Data
     integer (PGSLib_INT_TYPE),    POINTER, DIMENSION(:,:):: Supplement_2_I_Data
     real    (PGSLib_REAL_TYPE),   POINTER, DIMENSION(:,:):: Supplement_2_R_Data
     real    (PGSLib_DOUBLE_TYPE), POINTER, DIMENSION(:,:):: Supplement_2_D_Data
     logical (PGSLib_Log_TYPE),    POINTER, DIMENSION(:,:):: Supplement_2_L_Data

     integer (PGSLib_INT_TYPE),    POINTER, DIMENSION(:):: Duplicate_Indices
     integer (PGSLib_INT_TYPE),    POINTER, DIMENSION(:):: Duplicate_I_Data
     real    (PGSLib_REAL_TYPE),   POINTER, DIMENSION(:):: Duplicate_R_Data
     real    (PGSLib_DOUBLE_TYPE), POINTER, DIMENSION(:):: Duplicate_D_Data
     logical (PGSLib_Log_TYPE),    POINTER, DIMENSION(:):: Duplicate_L_Data
     integer (PGSLib_INT_TYPE),    POINTER, DIMENSION(:,:):: Duplicate_2_I_Data
     real    (PGSLib_REAL_TYPE),   POINTER, DIMENSION(:,:):: Duplicate_2_R_Data
     real    (PGSLib_DOUBLE_TYPE), POINTER, DIMENSION(:,:):: Duplicate_2_D_Data
     logical (PGSLib_Log_TYPE),    POINTER, DIMENSION(:,:):: Duplicate_2_L_Data


     integer (PGSLib_Int_Type),    pointer, dimension(:)  :: Global_Index_1

     ! Timers for tracking times and use of a trace
     type (PGSLib_Timer) :: GatherTotal_Timer, ScatterTotal_Timer
     type (PGSLib_Timer) :: GatherBuffer_Timer, ScatterBuffer_Timer
     type (PGSLib_Timer) :: Setup_Timer


     integer (PGSLib_INT_TYPE),    POINTER, DIMENSION(:) :: GS_Trace
     ! Here is the key for the array element''s meaning in the serial code.
     ! 1 = GS_Trace allocated
     ! 2 = Hash table initialized
     ! 3 = Hash table ready for redirection.

  END TYPE PGSLib_GS_TRACE

  ! This type defines the scope of some of the operations
  TYPE PGSLib_SCOPE
     integer (PGSLib_INT_Type) :: SCOPE
  END TYPE PGSLib_SCOPE
  
  ! This operator tests for equality of scope
  INTERFACE OPERATOR (.EQ.)
     MODULE PROCEDURE EQ_SCOPE
  END INTERFACE
  
CONTAINS
  function EQ_SCOPE(A, B)
    implicit none
    type (PGSLib_Scope), intent(IN) :: A, B
    logical (PGSLib_Log_Type)       :: EQ_SCOPE

    EQ_SCOPE = (A%Scope == B%Scope)

    RETURN
      end function EQ_SCOPE
  
  !=====================================================================
  !          Routines to Access fields in PGSLib_GS_Trace
  !=====================================================================

  function PGSLib_Size_Of_Dup(Trace)
    implicit none
    type (PGSLib_GS_Trace), INTENT(IN) :: Trace
    integer (PGSLib_Int_Type)           :: PGSLib_Size_Of_Dup
    PGSLib_Size_Of_Dup = Trace%N_Duplicate
    RETURN
  end function PGSLib_Size_Of_Dup

  function PGSLib_Size_Of_Sup(Trace)
    implicit none
    type (PGSLib_GS_Trace), INTENT(IN) :: Trace
    integer (PGSLib_Int_Type)           :: PGSLib_Size_Of_Sup
    PGSLib_Size_Of_Sup = Trace%N_Supplement
    RETURN
  end function PGSLib_Size_Of_Sup

  function PGSLib_Dup_Index_S(Trace, I)
    implicit none
    type (PGSLib_GS_Trace),   INTENT(IN) :: Trace
    integer (PGSLib_Int_Type), INTENT(IN) :: I
    integer (PGSLib_Int_Type)             :: PGSLib_Dup_Index_S
    PGSLib_Dup_Index_S = Trace%Duplicate_Indices(I)
    RETURN
  end function PGSLib_Dup_Index_S

  function PGSLib_Dup_Index_V(Trace)
    implicit none
    type (PGSLib_GS_Trace),   INTENT(IN)   :: Trace
    integer (PGSLib_Int_Type),               &
               DIMENSION(Trace%N_Duplicate) :: PGSLib_Dup_Index_V
    PGSLib_Dup_Index_V = Trace%Duplicate_Indices
    RETURN
  end function PGSLib_Dup_Index_V
  
  function PGSLib_Sup_Global_PEs(Trace)
    implicit none
    type (PGSLib_GS_Trace),   INTENT(IN)   :: Trace
    integer (PGSLib_Int_Type),               &
               DIMENSION(Trace%N_Supplement) :: PGSLib_Sup_Global_PEs
    PGSLib_Sup_Global_PEs = Trace%Supplement_Global_PEs
    RETURN
  end function PGSLib_Sup_Global_PEs
  
  function PGSLib_Sup_Global_Indices(Trace)
    implicit none
    type (PGSLib_GS_Trace),   INTENT(IN)   :: Trace
    integer (PGSLib_Int_Type),               &
               DIMENSION(Trace%N_Supplement) :: PGSLib_Sup_Global_Indices
    PGSLib_Sup_Global_Indices = Trace%Supplement_Global_Indices
    RETURN
  end function PGSLib_Sup_Global_Indices
  
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routines for manipulating error structures
  function PGSLib_Error_Test(Error)
    ! Any non-zero error number is an error
    implicit none
    type (PGSLib_Error_Type), intent(IN)  :: Error
    logical (PGSLib_Log_Type) :: PGSLib_Error_Test
    PGSLib_Error_Test = Error%Error_Number /= 0
    return
  end function PGSLib_Error_Test

  subroutine PGSLib_Error_Set(Error, Number)
    implicit none
    type (PGSLib_Error_Type), intent(INOUT) :: Error
    integer (PGSLib_Int_Type), intent(IN   ) :: Number
    Error%Error_Number = Number
    return
  end subroutine PGSLib_Error_Set
  
END MODULE PGSLib_Type_MODULE
