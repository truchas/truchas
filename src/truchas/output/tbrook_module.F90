!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! 
!!!  WARNING: Do not modify this file
!!!           Do not put use associations in this file
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE tbrook_module
  ! This module is a parallel wrapper around the brooks module.  We
  ! will do all the parallel calls here so that users do not have to
  ! worry about parallelization.  Note that no permutation is done 
  ! to indices.  The assumption is that the output will have the 
  ! information for doing the permutation.
  !
  ! A new data structure will be defined here for our wrapper called
  ! the TBROOK structure.  This is a private structure, not visible
  ! to the outside world.  We will define a set of integer parameters
  ! as IDs for specific output files.
  !
  ! Note that since we use only one output file, all of the different
  ! TBrooks defined here write to the same output file.  We just
  ! maintain the different brooks to have different tags for post
  ! processing.  
  !
  ! The only parallel call is in function I_SHALL_WRITE
  ! If you use PGSLib, you can directly use this module.
  ! If you use other parallelization, change the selection
  ! of the IOP in function I_SHALL_WRITE


  use brook_module

  ! The following are dependencies in the include files
  use pgslib_module,        only: PGSLib_GLOBAL_SUM
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: TBROOK
  PUBLIC :: BROOK             ! Provide the brook from brook_module
  PUBLIC :: BROOK_POINTER     ! Provide the brook from brook_module
  PUBLIC :: B_STDOUT          ! Provide the brook from brook_module
  PUBLIC :: NULLBROOK
  PUBLIC :: B_IFORM_BINARY
  PUBLIC :: B_IFORM_ASCII
  PUBLIC :: B_IFORM_XML
  PUBLIC :: BROOK_SET         ! Provide the brook from brook_module
  PUBLIC :: TBROOK_SET 
  PUBLIC :: ASSIGNMENT(=)
  PUBLIC :: OPERATOR(+)
  PUBLIC :: TBROOK_FLUSH

  PUBLIC :: OUT_TBROOK
  PUBLIC :: TTY_TBROOK
  PUBLIC :: ERR_TBROOK
  PUBLIC :: AUX_TBROOK
  PUBLIC :: INT_TBROOK
  PUBLIC :: OUT_TTY_TBROOK
  PUBLIC :: OUT_ERR_TBROOK
  PUBLIC :: OUT_ERR_TTY_TBROOK
  PUBLIC :: ERR_TTY_TBROOK
  PUBLIC :: XML_TBROOK
  PUBLIC :: BaseBrook
  PUBLIC :: BLastStep
  PUBLIC :: BaseBrookAscii

  PUBLIC :: TB_SCOPE_LOCAL
  PUBLIC :: TB_SCOPE_GLOBAL

  PUBLIC :: TBrook_initialize
  PUBLIC :: TBrook_SetBaseFile
  PUBLIC :: TBrook_CleanUp
  PUBLIC :: TBrook_Remove_Duplicates
  PUBLIC :: TBrook_GET_BROOK
  PUBLIC :: TBrook_NEXT
  PUBLIC :: TBrook_GET_SCOPE
  PUBLIC :: TBrook_GET_Tag
  PUBLIC :: TBrook_GET_LastTag
  PUBLIC :: TBrook_SET_LastTag
  PUBLIC :: TBROOK_WRITE
  PUBLIC :: TBROOK_WRITEXMLTAG
  PUBLIC :: TBROOK_OPENXMLTAG
  PUBLIC :: TBROOK_CLOSEXMLTAG
  PUBLIC :: TBROOK_WRITEXMLComment
  PUBLIC :: TBROOK_CLOSE
  PUBLIC :: TBROOK_DESTROY
  PUBLIC :: B_STRLEN
  PUBLIC :: TBROOK_FILE
  PUBLIC :: TBROOK_UNIT
  PUBLIC :: TBROOK_DEDUPE_IDS
  PUBLIC :: TBU_MAKE_FILE_ENTRY

  INTERFACE TBROOK_FLUSH
     module procedure BROOK_FLUSH
     module procedure TB_FLUSH
  END INTERFACE

  INTERFACE TBROOK_UNIT
     module procedure BROOK_UNIT
     module procedure TB_UNIT
  END INTERFACE

  INTERFACE TBROOK_FILE
     module procedure BROOK_FILE
     module procedure TB_FILE
  END INTERFACE

  INTERFACE TBROOK_REMOVE_DUPLICATES
     module procedure TB_ID_DEDUPE
     module procedure Brook_Remove_Duplicates
  END INTERFACE

  INTERFACE TBROOK_WRITEXMLCOMMENT
     MODULE PROCEDURE TB_ID_WXMLC_STRING
     MODULE PROCEDURE TB_ID_WXMLC_STRINGARRAY
     MODULE PROCEDURE B_WXMLC_STRING
     MODULE PROCEDURE B_WXMLC_STRINGARRAY
  END INTERFACE

  INTERFACE TBROOK_CLOSE
     MODULE PROCEDURE  TB_CLOSE
     MODULE PROCEDURE  TB_ID_CLOSE
     MODULE PROCEDURE  TB_BROOK_CLOSE
  END INTERFACE

  INTERFACE TBROOK_DESTROY
     MODULE PROCEDURE TB_DESTROY
     MODULE PROCEDURE TB_ID_DESTROY
     MODULE PROCEDURE TB_BROOK_DESTROY
  END INTERFACE
  ! Increase this number if you add any more brooks
  INTEGER, PARAMETER :: MAX_TBROOKS          = 10
  INTEGER, PARAMETER :: OUT_TBROOK           = 1
  INTEGER, PARAMETER :: ERR_TBROOK           = 2
  INTEGER, PARAMETER :: TTY_TBROOK           = 3
  INTEGER, PARAMETER :: AUX_TBROOK           = 4
  INTEGER, PARAMETER :: INT_TBROOK           = 5
  INTEGER, PARAMETER :: OUT_TTY_TBROOK       = 6
  INTEGER, PARAMETER :: ERR_TTY_TBROOK       = 7
  INTEGER, PARAMETER :: OUT_ERR_TBROOK       = 8
  INTEGER, PARAMETER :: XML_TBROOK           = 9
  INTEGER, PARAMETER :: OUT_ERR_TTY_TBROOK   = 10
  
  INTEGER, PARAMETER :: TB_SCOPE_LOCAL       = 1
  INTEGER, PARAMETER :: TB_SCOPE_GLOBAL      = 2
  
  INTEGER, PARAMETER :: TB_TAGLEN = 32
  INTEGER, PARAMETER :: B_STRLEN  = 1024

  CHARACTER(LEN=32), PARAMETER :: DO_NOT_INITIALIZE = 'uninitialized_file'
  LOGICAL, SAVE :: no_file_set = .true.
  TYPE TBROOK
     CHARACTER(LEN=TB_TAGLEN) :: LASTTAG = 'UNKNOWN'
     CHARACTER(LEN=TB_TAGLEN) :: TAG     = 'UNKNOWN'
     TYPE(BROOK), pointer     :: B => NULL()
  END TYPE TBROOK

  public :: TB
  public :: MAX_TBROOKS
  TYPE(TBROOK), TARGET, ALLOCATABLE, DIMENSION(:), SAVE :: TB
  TYPE(Brook),  TARGET,                            SAVE :: NullBrook
  TYPE(Brook),  TARGET,                            SAVE :: BaseBrook
  TYPE(Brook),  TARGET,                            SAVE :: BLastStep
  TYPE(Brook),  TARGET,                            SAVE :: BaseBrookAscii
  LOGICAL,      TARGET, ALLOCATABLE, DIMENSION(:), SAVE :: TB_OPENTAG
  LOGICAL, SAVE :: initialized = .FALSE.
  CHARACTER(LEN=32),        save :: TBrookVersion='1.001'
  CHARACTER(LEN=TB_TAGLEN), save :: TB_LAST_TAG  = 'UNKNOWN_TAG'
  
  ! The following are needed for the TBrook_Write interface, but the
  ! average Joe need never know that they exist.  Hence they are
  ! tucked in here.
  PUBLIC :: BROOK_CLOSE
  PUBLIC :: BROOK_DESTROY
  PUBLIC :: TB_CLOSE
  PUBLIC :: TB_ID_CLOSE
  PUBLIC :: TB_BROOK_CLOSE
  PUBLIC :: TB_DESTROY
  PUBLIC :: TB_ID_DESTROY
  PUBLIC :: TB_BROOK_DESTROY

  PUBLIC :: TB_WRITE_INTEGER_R0
  PUBLIC :: TBB_WRITE_INTEGER_R0
  PUBLIC :: TB_WRITE_LOGICAL_R0
  PUBLIC :: TBB_WRITE_LOGICAL_R0
  PUBLIC :: TB_WRITE_DOUBLE_R0
  PUBLIC :: TBB_WRITE_DOUBLE_R0
  PUBLIC :: TB_WRITE_FLOAT_R0
  PUBLIC :: TBB_WRITE_FLOAT_R0
  PUBLIC :: TB_WRITE_STRING_R0
  PUBLIC :: TBB_WRITE_STRING_R0

  PUBLIC :: TB_WRITE_INTEGER_R1
  PUBLIC :: TBB_WRITE_INTEGER_R1
  PUBLIC :: TB_WRITE_LOGICAL_R1
  PUBLIC :: TBB_WRITE_LOGICAL_R1
  PUBLIC :: TB_WRITE_DOUBLE_R1
  PUBLIC :: TBB_WRITE_DOUBLE_R1
  PUBLIC :: TB_WRITE_STRING_R1
  PUBLIC :: TBB_WRITE_STRING_R1
  PUBLIC :: TB_WRITE_FLOAT_R1
  PUBLIC :: TBB_WRITE_FLOAT_R1

  PUBLIC :: TB_WRITE_INTEGER_R2
  PUBLIC :: TBB_WRITE_INTEGER_R2
  PUBLIC :: TB_WRITE_LOGICAL_R2
  PUBLIC :: TBB_WRITE_LOGICAL_R2
  PUBLIC :: TB_WRITE_DOUBLE_R2
  PUBLIC :: TBB_WRITE_DOUBLE_R2
  PUBLIC :: TB_WRITE_FLOAT_R2
  PUBLIC :: TBB_WRITE_FLOAT_R2

  PUBLIC :: TB_WRITE_INTEGER_R3
  PUBLIC :: TBB_WRITE_INTEGER_R3
  PUBLIC :: TB_WRITE_LOGICAL_R3
  PUBLIC :: TBB_WRITE_LOGICAL_R3
  PUBLIC :: TB_WRITE_DOUBLE_R3
  PUBLIC :: TBB_WRITE_DOUBLE_R3
  PUBLIC :: TB_WRITE_FLOAT_R3
  PUBLIC :: TBB_WRITE_FLOAT_R3

  PUBLIC :: TB_WRITE_INTEGER_R4
  PUBLIC :: TBB_WRITE_INTEGER_R4
  PUBLIC :: TB_WRITE_LOGICAL_R4
  PUBLIC :: TBB_WRITE_LOGICAL_R4
  PUBLIC :: TB_WRITE_DOUBLE_R4
  PUBLIC :: TBB_WRITE_DOUBLE_R4
  PUBLIC :: TB_WRITE_FLOAT_R4
  PUBLIC :: TBB_WRITE_FLOAT_R4

  PUBLIC :: TB_WRITE_INTEGER_R5
  PUBLIC :: TBB_WRITE_INTEGER_R5
  PUBLIC :: TB_WRITE_LOGICAL_R5
  PUBLIC :: TBB_WRITE_LOGICAL_R5
  PUBLIC :: TB_WRITE_DOUBLE_R5
  PUBLIC :: TBB_WRITE_DOUBLE_R5
  PUBLIC :: TB_WRITE_FLOAT_R5
  PUBLIC :: TBB_WRITE_FLOAT_R5

!  PUBLIC :: B_ASSIGN
!  PUBLIC :: B_ASSIGN_STRUCT_POINTER
  public :: TB_ASSIGN_TB_TO_TB
  public :: TB_ASSIGN_TB_TO_B
  public :: TB_ASSIGN_B_TO_TB
  public :: TB_ASSIGN_INT_TO_TB
  public :: TB_ASSIGN_INT_TO_B

  public :: TBrook_Endline
  public :: TB_B_Endline
  Public :: TB_ID_Endline

  interface Tbrook_Endline
     module procedure TB_B_Endline
     module procedure TB_ID_Endline
  end interface

  INTERFACE TBROOK_SET
     MODULE PROCEDURE BROOK_SET
     MODULE PROCEDURE TB_SET
  END INTERFACE

  INTERFACE TBROOK_WRITE

     MODULE PROCEDURE TB_WRITE_INTEGER_R0
     MODULE PROCEDURE TBB_WRITE_INTEGER_R0
     MODULE PROCEDURE TB_WRITE_LOGICAL_R0
     MODULE PROCEDURE TBB_WRITE_LOGICAL_R0
     MODULE PROCEDURE TB_WRITE_DOUBLE_R0
     MODULE PROCEDURE TBB_WRITE_DOUBLE_R0
     MODULE PROCEDURE TB_WRITE_FLOAT_R0
     MODULE PROCEDURE TBB_WRITE_FLOAT_R0
     MODULE PROCEDURE TB_WRITE_STRING_R0
     MODULE PROCEDURE TBB_WRITE_STRING_R0

     MODULE PROCEDURE TB_WRITE_INTEGER_R1
     MODULE PROCEDURE TBB_WRITE_INTEGER_R1
     MODULE PROCEDURE TB_WRITE_LOGICAL_R1
     MODULE PROCEDURE TBB_WRITE_LOGICAL_R1
     MODULE PROCEDURE TB_WRITE_DOUBLE_R1
     MODULE PROCEDURE TBB_WRITE_DOUBLE_R1
     MODULE PROCEDURE TB_WRITE_STRING_R1
     MODULE PROCEDURE TBB_WRITE_STRING_R1
     MODULE PROCEDURE TB_WRITE_FLOAT_R1
     MODULE PROCEDURE TBB_WRITE_FLOAT_R1

     MODULE PROCEDURE TB_WRITE_INTEGER_R2
     MODULE PROCEDURE TBB_WRITE_INTEGER_R2
     MODULE PROCEDURE TB_WRITE_LOGICAL_R2
     MODULE PROCEDURE TBB_WRITE_LOGICAL_R2
     MODULE PROCEDURE TB_WRITE_DOUBLE_R2
     MODULE PROCEDURE TBB_WRITE_DOUBLE_R2
     MODULE PROCEDURE TB_WRITE_FLOAT_R2
     MODULE PROCEDURE TBB_WRITE_FLOAT_R2

     MODULE PROCEDURE TB_WRITE_INTEGER_R3
     MODULE PROCEDURE TBB_WRITE_INTEGER_R3
     MODULE PROCEDURE TB_WRITE_LOGICAL_R3
     MODULE PROCEDURE TBB_WRITE_LOGICAL_R3
     MODULE PROCEDURE TB_WRITE_DOUBLE_R3
     MODULE PROCEDURE TBB_WRITE_DOUBLE_R3
     MODULE PROCEDURE TB_WRITE_FLOAT_R3
     MODULE PROCEDURE TBB_WRITE_FLOAT_R3

     MODULE PROCEDURE TB_WRITE_INTEGER_R4
     MODULE PROCEDURE TBB_WRITE_INTEGER_R4
     MODULE PROCEDURE TB_WRITE_LOGICAL_R4
     MODULE PROCEDURE TBB_WRITE_LOGICAL_R4
     MODULE PROCEDURE TB_WRITE_DOUBLE_R4
     MODULE PROCEDURE TBB_WRITE_DOUBLE_R4
     MODULE PROCEDURE TB_WRITE_FLOAT_R4
     MODULE PROCEDURE TBB_WRITE_FLOAT_R4

     MODULE PROCEDURE TB_WRITE_INTEGER_R5
     MODULE PROCEDURE TBB_WRITE_INTEGER_R5
     MODULE PROCEDURE TB_WRITE_LOGICAL_R5
     MODULE PROCEDURE TBB_WRITE_LOGICAL_R5
     MODULE PROCEDURE TB_WRITE_DOUBLE_R5
     MODULE PROCEDURE TBB_WRITE_DOUBLE_R5
     MODULE PROCEDURE TB_WRITE_FLOAT_R5
     MODULE PROCEDURE TBB_WRITE_FLOAT_R5
  END INTERFACE
  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE TB_ASSIGN_B_TO_B
     MODULE PROCEDURE TB_ASSIGN_BP_TO_B
     MODULE PROCEDURE TB_ASSIGN_TB_TO_TB
     MODULE PROCEDURE TB_ASSIGN_TB_TO_B
     MODULE PROCEDURE TB_ASSIGN_B_TO_TB
     MODULE PROCEDURE TB_ASSIGN_INT_TO_TB
     MODULE PROCEDURE TB_ASSIGN_INT_TO_B
  END INTERFACE
  interface tbu_make_file_entry
     module procedure tbu_fileentry_i0
     module procedure tbu_fileentry_r0
     module procedure tbu_fileentry_l0
     module procedure tbu_fileentry_d0
     module procedure tbu_fileentry_c0
     module procedure tbu_fileentry_i1
     module procedure tbu_fileentry_r1
     module procedure tbu_fileentry_l1
     module procedure tbu_fileentry_d1
     module procedure tbu_fileentry_c1
     module procedure tbu_fileentry_i2
     module procedure tbu_fileentry_r2
     module procedure tbu_fileentry_l2
     module procedure tbu_fileentry_d2
!!$     module procedure tbu_fileentry_c2
     module procedure tbu_fileentry_i3
     module procedure tbu_fileentry_r3
     module procedure tbu_fileentry_l3
     module procedure tbu_fileentry_d3
!!$     module procedure tbu_fileentry_c3
     module procedure tbu_fileentry_i4
     module procedure tbu_fileentry_r4
     module procedure tbu_fileentry_l4
     module procedure tbu_fileentry_d4
!!$     module procedure tbu_fileentry_c4
     module procedure tbu_fileentry_i5
     module procedure tbu_fileentry_r5
     module procedure tbu_fileentry_l5
     module procedure tbu_fileentry_d5
!!$     module procedure tbu_fileentry_c5
  end interface
CONTAINS

  SUBROUTINE TB_FLUSH(TB_ID, iStatus)
    integer, intent(in) :: TB_ID
    integer, intent(inout) :: iStatus

    iStatus = 0
    if ( initialized .and. TB_ID .gt. 0  .and. TB_ID .le. MAX_TBROOKS) then
       call Brook_Flush(TB(TB_ID)%B, iStatus)
    end if
    return
  END SUBROUTINE TB_FLUSH

  SUBROUTINE tbrook_dedupe_ids(B, UArray, Tag)
    implicit none
    type(Brook), target :: B
    integer, dimension(:) :: UArray
    character(LEN=*) :: Tag

    type(Brook), pointer :: P
    integer :: iSize
    integer :: i
    logical, dimension(MAX_TBROOKS) :: lArray
    logical :: b_flag
    
    lArray = .false.
    b_flag = .false.
    iSize = size(UArray)

    Tag = ' '

    do i = 1, iSize
       select case (UArray(i))
       case (OUT_TBROOK)
          b_flag = .true.
          if ( .not. lArray(OUT_TBROOK)) then
             P=> TBrook_Next(OUT_TBrook)
             if(ASSOCIATED(P)) B = B + P
          end if
          lArray(OUT_TBROOK)=.true.
       case (ERR_TBROOK)
          b_flag = .true.
          if ( .not. lArray(ERR_TBROOK)) then
             P=> TBrook_Next(ERR_TBrook)
             if(ASSOCIATED(P)) B = B + P
          end if
          lArray(ERR_TBROOK)=.true.
       case (TTY_TBROOK)
          if ( .not. lArray(TTY_TBROOK)) then
             B = B + TBrook_Get_Brook(TTY_TBrook)
          end if
          lArray(TTY_TBROOK)=.true.
       case (OUT_ERR_TBROOK)
          b_flag = .true.
          if ( .not. lArray(OUT_TBROOK)) then
             P=> TBrook_Next(OUT_TBrook)
             if(ASSOCIATED(P)) B = B + P
          end if
          lArray(OUT_TBROOK)=.true.
          if ( .not. lArray(ERR_TBROOK)) then
             P=> TBrook_Next(ERR_TBrook)
             if(ASSOCIATED(P)) B = B + P
          end if
          lArray(ERR_TBROOK)=.true.
       case (OUT_TTY_TBROOK)
          b_flag = .true.
          if ( .not. lArray(OUT_TBROOK)) then
             P=> TBrook_Next(OUT_TBrook)
             if(ASSOCIATED(P)) B = B + P
          end if
          lArray(OUT_TBROOK)=.true.
          if ( .not. lArray(TTY_TBROOK)) then
             B = B + TBrook_Get_Brook(TTY_TBrook)
          end if
          lArray(TTY_TBROOK)=.true.
       case (OUT_ERR_TTY_TBROOK)
          b_flag = .true.
          if ( .not. lArray(OUT_TBROOK)) then
             P=> TBrook_Next(OUT_TBrook)
             if(ASSOCIATED(P)) B = B + P
          end if
          lArray(OUT_TBROOK)=.true.
          if ( .not. lArray(ERR_TBROOK)) then
             P=> TBrook_Next(ERR_TBrook)
             if(ASSOCIATED(P)) B = B + P
          end if
          lArray(ERR_TBROOK)=.true.
          if ( .not. lArray(TTY_TBROOK)) then
             B = B + TBrook_Get_Brook(TTY_TBrook)
          end if
          lArray(TTY_TBROOK)=.true.
       case (AUX_TBROOK)
          if ( .not. lArray(AUX_TBROOK)) then
             B = B + TBrook_Get_Brook(AUX_TBrook)
          end if
          lArray(AUX_TBROOK)=.true.
       case (INT_TBROOK)
          if ( .not. lArray(INT_TBROOK)) then
             B = B + TBrook_Get_Brook(INT_TBrook)
          end if
          lArray(INT_TBROOK)=.true.
       end select
    end do
    if ( b_flag ) then
       B = B + BaseBrookAscii
    end if


!!$    P => B
!!$    Nullify(P)
!!$    do
!!$       if ( .not. associated(P)) exit
!!$       N => Brook_Next(P)
!!$       if ( .not. associated(N)) exit
!!$       if ( associated(P%Data) .and. associated(N%Data)) then
!!$          if(TRIM(Brook_File(P)) == TRIM(Brook_File(N))) then
!!$             P%Next => N%Next
!!$             deallocate(N%Data)
!!$             Nullify(N%Next)
!!$             deallocate(N)
!!$          else 
!!$             P=>Brook_Next(P)
!!$          end if
!!$       else 
!!$          P=>Brook_Next(P)
!!$       end if
!!$    end do

    TAG=' '
    do i = 1, MAX_TBROOKS
       if ( LArray(i) ) TAG = TRIM(TAG)//TRIM(TBrook_Get_Bare_Tag(i))
    end do

  END SUBROUTINE Tbrook_dedupe_ids



  SUBROUTINE TB_ID_DEDUPE(TB_ID, iStatus)
    ! Remove duplicate brooks from TB(TB_ID)
    integer, intent(in)    :: TB_ID
    integer, intent(inout) :: iStatus
    iStatus = 0
    if ( initialized .and. TB_ID .gt. 0  .and. TB_ID .le. MAX_TBROOKS) then
       call Brook_remove_Duplicates(TB(TB_ID)%B, iStatus)
    end if
    return
  END SUBROUTINE TB_ID_DEDUPE
            

  SUBROUTINE TB_ASSIGN_B_TO_B(LHS, RHS)
    use brook_module, only: brook, b_assign
    type(BROOK), target, intent(INOUT) :: LHS
    type(BROOK), target, intent(IN)    :: RHS
    integer :: iStatus
    iStatus = 0
    IF ( .not. initialized ) then
       call TBrook_Initialize(iStatus)
       if ( iStatus /= 0 ) return
    END IF
    call B_Assign(LHS, RHS)
    return
  END SUBROUTINE TB_ASSIGN_B_TO_B

  SUBROUTINE TB_ASSIGN_BP_TO_B(LHS, RHS)
    use brook_module, only: brook, b_assign_struct_pointer,  brook_pointer
    type(BROOK), target, intent(INOUT) :: LHS
    type(BROOK_POINTER), target, intent(IN)    :: RHS
    integer :: iStatus
    iStatus = 0
    IF ( .not. initialized ) then
       call TBrook_Initialize(iStatus)
       if ( iStatus /= 0 ) return
    END IF
    call B_Assign_Struct_pointer(LHS, RHS)
    return
  END SUBROUTINE TB_ASSIGN_BP_TO_B
    

  SUBROUTINE TB_ASSIGN_TB_TO_TB(LHS, RHS)
    type(TBROOK), target, intent(INOUT) :: LHS
    type(TBROOK), target, intent(IN) :: RHS
    integer :: iStatus

    iStatus = 0

    IF ( .not. initialized ) then
       call TBrook_Initialize(iStatus)
       if ( iStatus /= 0 ) return
    END IF
    
    call TBROOK_DESTROY(LHS, iStatus=iStatus)
    ALLOCATE(LHS%B, Stat=iStatus)
    LHS%TAG = RHS%TAG
    LHS%B = RHS%B

    return
  END SUBROUTINE TB_ASSIGN_TB_TO_TB

  SUBROUTINE TB_ASSIGN_INT_TO_TB(LHS, RHS)
    type(TBROOK), target, intent(INOUT) :: LHS
    integer, intent(IN) :: RHS
    integer :: iStatus

    iStatus = 0
    IF ( .not. initialized ) then
       call TBrook_Initialize(iStatus)
       if ( iStatus /= 0 ) return
    END IF
    if ( RHS > 0 .and. RHS <= MAX_TBROOKS) then
       call TB_ASSIGN_TB_TO_TB(LHS,TB(RHS))
    end if

    return
  END SUBROUTINE TB_ASSIGN_INT_TO_TB

  SUBROUTINE TB_ASSIGN_INT_TO_B(LHS, RHS)
    type(BROOK), target, intent(INOUT) :: LHS
    integer, intent(IN) :: RHS
    integer :: iStatus

    iStatus = 0
    IF ( .not. initialized ) then
       call TBrook_Initialize(iStatus)
       if ( iStatus /= 0 ) return
    END IF

    if ( RHS > 0 .and. RHS <= MAX_TBROOKS) then
       call TB_ASSIGN_TB_TO_B(LHS,TB(RHS))
    end if

    return
  END SUBROUTINE TB_ASSIGN_INT_TO_B

  SUBROUTINE TB_ASSIGN_B_TO_TB(LHS, RHS)
    type(TBROOK), target, intent(INOUT) :: LHS
    type(BROOK),  target, intent(IN)    :: RHS
    integer :: iStatus
    CHARACTER(LEN=TB_TAGLEN) :: tag
    iStatus = 0

    IF ( .not. initialized ) then
       call TBrook_Initialize(iStatus)
       if ( iStatus /= 0 ) return
    END IF
    tag = LHS%TAG
    call TBROOK_DESTROY(LHS, iStatus=iStatus)
    ALLOCATE(LHS%B, Stat=iStatus)
    LHS%TAG = tag
    LHS%B = RHS
    return
  END SUBROUTINE TB_ASSIGN_B_TO_TB

  SUBROUTINE TB_ASSIGN_TB_TO_B(LHS, RHS)
    type(BROOK),  target, intent(INOUT) :: LHS
    type(TBROOK), target, intent(IN)    :: RHS
    integer :: iStatus
    iStatus = 0
    IF ( .not. initialized ) then
       call TBrook_Initialize(iStatus)
       if ( iStatus /= 0 ) return
    END IF
    call TBROOK_DESTROY(LHS, iStatus=iStatus)
    LHS = RHS%B
    return
  END SUBROUTINE TB_ASSIGN_TB_TO_B
    
     LOGICAL FUNCTION I_SHALL_WRITE(SCOPE)
       !
       ! Determine who shall write to brooks.
       ! Only the parallel io processor and 
       ! local brooks are allowed to write to 
       ! to brooks
       use parallel_info_module, only: p_info
       INTEGER, OPTIONAL, INTENT(IN) :: SCOPE

       INTEGER :: myScope

       if (Present(SCOPE)) then
          myScope = SCOPE
       else 
          myScope = TB_SCOPE_GLOBAL
       end if

       if ( (myScope == TB_SCOPE_LOCAL) .or. (p_info%IOP) ) then
          I_SHALL_WRITE = .true.
       else 
          I_SHALL_WRITE = .false.
       end if
       return
     end FUNCTION I_SHALL_WRITE


       
  SUBROUTINE TB_Set(TB_id,           &
                    File,            &
                    Unit,            &
                    Form,            &
                    Position,        & 
                    StandardHeaders, &
                    Closeable,       &
                    IFormat,         &
                    DFormat,         &
                    LFormat,         &
                    AFormat,         &
                    XMLRoot,         &
                    iStatus)
    use brook_module, only: brook_set
    IMPLICIT NONE
    
    ! Argument List
    integer,                    intent(in)    :: TB_ID
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: File
    INTEGER,          OPTIONAL, INTENT(IN)    :: Unit
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: Form
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: Position
    LOGICAL,          OPTIONAL, INTENT(IN)    :: StandardHeaders
    LOGICAL,          OPTIONAL, INTENT(IN)    :: Closeable
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: XMLRoot
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: Iformat
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: Dformat
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: Lformat
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: Aformat
    INTEGER,                    INTENT(INOUT) :: iStatus

    if ( TB_ID .gt. 0  .and. TB_ID .le. MAX_TBROOKS) then
       call Brook_Set(TB(TB_ID)%B,                     &
                      File=File,                       &
                      Unit=Unit,                       &
                      Form=Form,                       &
                      Position=Position,               &
                      StandardHeaders=StandardHeaders, &
                      Closeable=CLoseable,             &
                      XMLRoot=XMLRoot,                 &
                      IFormat=IFormat,                 &
                      DFormat=DFormat,                 &
                      LFormat=LFormat,                 &
                      AFormat=AFormat,                 &
                      iStatus=iStatus)
    end if
    return
  END SUBROUTINE TB_Set


  Subroutine TB_ID_Endline(TB_ID, iStatus)
    integer, intent(in)    :: TB_ID
    integer, intent(inout) :: iStatus
    if ( TB_ID .gt. 0  .and. TB_ID .le. MAX_TBROOKS) then
       call TB_B_Endline(TBrook_Get_Brook(TB_ID), iStatus=iStatus)
    end if
  End Subroutine TB_ID_Endline

  Subroutine TB_B_Endline(B, iStatus)
    use brook_module, only: Brook_Iform, &
                            B_Iform_Binary
    type(Brook), intent(IN) :: B
    integer, intent(inout)  :: iStatus
    
    if ( BROOK_IFORM(B) /= B_IFORM_BINARY) then
       ! Only do this if we have a non binary brook
       call TBrook_Write(B, Variable='', ADVANCE=.true., iStatus=iStatus)
    end if
    return
  End Subroutine TB_B_Endline
    
     
    FUNCTION TBrook_GET_BROOK(ID) RESULT(ptr)
      ! Return a pointer to the brook of given TBROOK ID
      integer, intent(in)  :: ID
      Type(Brook), pointer :: ptr

      if ( initialized .and. ID .ge. 1 .and. ID .le. MAX_TBROOKS ) then
         ptr => TB(ID)%B
      else 
         ptr => NULL()
      end if
      RETURN
    END FUNCTION TBROOK_GET_BROOK

    FUNCTION TBrook_Next(ID) RESULT(ptr)
      ! Return a pointer to the brook of given TBROOK ID
      integer, intent(in)  :: ID
      Type(Brook), pointer :: ptr
      
      Nullify(ptr)
      if ( initialized .and. ID .ge. 1 .and. ID .le. MAX_TBROOKS ) then
         ptr => Brook_Next(TB(ID)%B)
      end if
      RETURN
    END FUNCTION TBROOK_NEXT

    INTEGER FUNCTION TBrook_GET_SCOPE(B, SCOPE)
      ! Return either TB_SCOPE_GLOBAL or TB_SCOPE_LOCAL
      
      Type(Brook), optional, intent(IN) :: B ! ignored for now
      integer,     optional, intent(in) :: SCOPE

      if ( Present(SCOPE) ) then
         TBrook_Get_Scope = SCOPE
      else 
         TBrook_Get_Scope = TB_SCOPE_GLOBAL
      end if
      return
    END FUNCTION TBrook_GET_SCOPE

    FUNCTION TBrook_GET_LASTTAG() 
      ! Return a pointer to the brook of given TBROOK ID
      Character(LEN=TB_TAGLEN)          :: TBrook_Get_LastTag

      TBrook_Get_LastTag = TRIM(TB_Last_Tag)
      RETURN
    END FUNCTION TBROOK_GET_LASTTAG

    Subroutine TBrook_SET_LASTTAG(ID,TAG)
      ! Return a pointer to the brook of given TBROOK ID
      integer, optional, intent(in) :: ID
      character(len=*), optional, intent(in) :: TAG

      if ( present(ID) ) then
         TB_Last_Tag = TBrook_Get_Tag(ID)
      else if (present(TAG) ) then
         TB_Last_Tag = TAG
      end if

      RETURN
    END Subroutine TBrook_SET_LASTTAG

    FUNCTION TBrook_GET_TAG(ID) Result(ptr)
      ! Return a pointer to the brook of given TBROOK ID
      integer,               intent(in) :: ID
      Character(LEN=TB_TAGLEN), pointer :: ptr
      character(LEN=TB_TAGLEN), target, save :: myTag

      myTag = '['//TRIM(TBrook_Get_Bare_Tag(ID))//']'
      ptr => myTag
      RETURN
    END FUNCTION TBROOK_GET_TAG
    FUNCTION TBrook_GET_BARE_TAG(ID) Result(ptr)
      ! Return a pointer to the brook of given TBROOK ID
      integer,               intent(in) :: ID
      Character(LEN=TB_TAGLEN), pointer :: ptr
      character(LEN=TB_TAGLEN), target, save :: myTag='UNKNOWN'

!!$      IF ( .not. initialized ) then
!!$         call TBrook_Initialize(iStatus)
!!$         if ( iStatus /= 0 ) return
!!$      END IF

      if ( initialized .and. ID .ge. 1 .and. ID .le. MAX_TBROOKS ) then
         ptr => TB(ID)%Tag
      else 
         ptr => myTag
      end if
      RETURN
    END FUNCTION TBROOK_GET_BARE_TAG

    SUBROUTINE TBROOK_CLEANUP(ISTATUS)
      use brook_module, only: Brook_Close,  &
                              Brook_Destroy
      ! Release any memory we have allocated,
      ! Close the main brook
      integer, intent(inout) :: iStatus
      integer :: i
      iStatus = 0
      if ( .not. initialized) return
      
      do i = 1, MAX_TBROOKS
         call Brook_Destroy(TB(i)%B, iStatus=iStatus)
         deallocate(TB(i)%B, STAT=iSTATUS)
      end do
      deallocate(TB, STAT=iSTATUS)
      deallocate(TB_OPENTAG, STAT=iSTATUS)
      
      call Brook_Close(BaseBrook, iSTATUS)
      call Brook_Destroy(BaseBrook, iSTATUS)

      return
    end SUBROUTINE TBROOK_CLEANUP

    SUBROUTINE TB_ID_CLOSE(ID, ISTATUS)
      use brook_module, only: Brook_Close
      INTEGER, intent(IN) :: ID
      integer, intent(INOUT) :: iStatus
      type(brook), pointer :: pB
      !
      if ( iStatus /= 0 ) return
      pB=>TBrook_Get_Brook(ID)
      call Brook_Close(pB, iStatus)
      return
    end SUBROUTINE TB_ID_CLOSE



    SUBROUTINE TB_BROOK_CLOSE(B, ISTATUS)
      use brook_module, only: Brook_Close
      type(brook), target, intent(IN) :: B
      integer, intent(INOUT) :: iStatus
      type(brook), pointer :: pB
      !
      if ( iStatus /= 0 ) return
      pB=>B
      call Brook_Close(pB, iStatus)
      return
    end SUBROUTINE TB_BROOK_CLOSE

    SUBROUTINE TB_CLOSE(T, ISTATUS)
      use brook_module, only: Brook_Close
      Type(TBrook), target, intent(in):: T
      integer, intent(INOUT) :: iStatus
      type(tbrook), pointer :: pT
      !
      pt=>T
      if ( iStatus /= 0 ) return
      call Brook_Close(pT%B, iStatus)
      return
    end SUBROUTINE TB_CLOSE

    SUBROUTINE TB_ID_DESTROY(ID, ISTATUS)
      use brook_module, only: Brook_Destroy
      INTEGER, intent(IN) :: ID
      integer, intent(INOUT) :: iStatus
      type(brook), pointer :: pB
      !

      if ( iStatus /= 0 ) return
      if ( ID .gt. 0  .and. ID .le. MAX_TBROOKS) then
         pB=>TBrook_Get_Brook(ID)
         call Brook_Destroy(pB, iStatus)
         TB(ID)%B=>NULL()
      end if
      return
    end SUBROUTINE TB_ID_DESTROY

    SUBROUTINE TB_BROOK_Destroy(B, ISTATUS)
      use brook_module, only: Brook_Destroy
      type(brook), target, intent(IN) :: B
      integer, intent(INOUT) :: iStatus
      type(brook), pointer :: pB
      !
      if ( iStatus /= 0 ) return
      pB=>B
      call Brook_Destroy(pB, iStatus)
      return
    end SUBROUTINE TB_BROOK_DESTROY

    INTEGER FUNCTION TB_UNIT(TB_ID)
      integer :: TB_ID
      TB_UNIT = -1
      if ( TB_ID .gt. 0  .and. TB_ID .le. MAX_TBROOKS) then
         TB_UNIT = BROOK_UNIT(TBROOK_GET_BROOK(TB_ID))
      end if
      return
    END FUNCTION TB_UNIT

    CHARACTER(LEN=B_STRLEN) FUNCTION TB_FILE(TB_ID)
      integer :: TB_ID
      TB_FILE = 'invalid file'
      if ( TB_ID .gt. 0  .and. TB_ID .le. MAX_TBROOKS) then
         TB_FILE = Brook_FILE(TB(TB_ID)%B)
      end if
      return
    END FUNCTION TB_FILE


    SUBROUTINE TB_DESTROY(T, ISTATUS)
      use brook_module, only: Brook_Destroy
      Type(TBrook), target, intent(inout):: T
      integer, intent(INOUT) :: iStatus
      type(tbrook), pointer :: pT
      !
      pt=>T

      if ( iStatus /= 0 ) return
      if ( associated(pT%B) ) call Brook_Destroy(pT%B, iStatus)
      pT%TAG = 'UNINITIALIZED'
      pT%B => NULL()
      return
    end SUBROUTINE TB_DESTROY


    SUBROUTINE TBROOK_INITIALIZE(ISTATUS)
      use brook_module, only: b_initialize, b_stdout
      INTEGER :: iStatus
      integer :: i
      
      ! Set error to 0
      iStatus = 0

      ! Do nothing if already initialized
      IF ( initialized ) RETURN

      call b_initialize(iStatus)
      IF ( iStatus /= 0 ) RETURN
      initialized = .true.

      BaseBrook=B_Stdout

      ALLOCATE(TB_OPENTAG(MAX_TBROOKS), STAT=iStatus)
      IF ( iStatus /= 0 ) RETURN
      TB_OPENTAG = .false.

      ! Allocate space for the data structures
      ALLOCATE(TB(MAX_TBROOKS), STAT=iStatus)
      IF ( iStatus /= 0 ) THEN
         DEALLOCATE(TB_OPENTAG)
         RETURN
      END IF

       ! Now copy the structure into the rest of the brooks so that

      do i = 1, MAX_TBROOKS
         allocate(TB(i)%B, STAT=iSTATUS)
      end do

       ! we can write in multiple forms and extract what we need from
       ! the different forms
       
       TB(OUT_TBROOK)%TAG = 'OUT'
       CALL TBrook_Set(TB(OUT_TBROOK)%B, form="ascii", iStatus=iStatus)

       TB(TTY_TBROOK)%TAG = 'TTY'
       TB(TTY_TBROOK)%B = B_Stdout
       CALL TBrook_Set(TB(TTY_TBROOK)%B, form="ascii", iStatus=iStatus)

       TB(ERR_TBROOK)%TAG = 'ERR'
       CALL TBrook_Set(TB(ERR_TBROOK)%B, form="ascii", iStatus=iStatus)

       TB(AUX_TBROOK)%TAG = 'AUX'
       CALL TBrook_Set(TB(AUX_TBROOK)%B, form="ascii", iStatus=iStatus)

       TB(INT_TBROOK)%TAG = 'INT'
       CALL TBrook_Set(TB(INT_TBROOK)%B, form="ascii", iStatus=iStatus)

       TB(OUT_TTY_TBROOK)%TAG = 'OUTTTY'
       TB(OUT_TTY_TBROOK)%B = B_Stdout
       CALL TBrook_Set(TB(OUT_TTY_TBROOK)%B, form="ascii", iStatus=iStatus)

       TB(OUT_ERR_TBROOK)%TAG = 'OUTERR'
       CALL TBrook_Set(TB(OUT_ERR_TBROOK)%B, form="ascii", iStatus=iStatus)

       TB(ERR_TTY_TBROOK)%TAG = 'ERRTTY'
       TB(ERR_TTY_TBROOK)%B = B_Stdout
       CALL TBrook_Set(TB(ERR_TTY_TBROOK)%B, form="ascii", iStatus=iStatus)

       TB(XML_TBROOK)%TAG = 'TXML'
       CALL TBrook_Set(TB(XML_TBROOK)%B, form="ascii", iStatus=iStatus)

       TB(OUT_ERR_TTY_TBROOK)%TAG = 'OUTERRTTY'
       TB(OUT_ERR_TTY_TBROOK)%B = B_Stdout
       CALL TBrook_Set(TB(OUT_ERR_TTY_TBROOK)%B, form="ascii", iStatus=iStatus)


       RETURN
     END SUBROUTINE TBROOK_INITIALIZE

     
     SUBROUTINE TBROOK_SETBASEFILE(File, ISTATUS)
       use brook_module, only: brook_unit
       CHARACTER(LEN=*), OPTIONAL :: File
       INTEGER :: iStatus
       INTEGER :: myUnit
       LOGICAL :: first_time = .false.
       
       ! Set error to 0
       iStatus = 0

       if (.not. first_time ) then
          myUnit = Brook_AssignUnit()
          first_time = .true.
       else 
          myUnit = Brook_Unit(BaseBrook)
       end if
          
       ! Do nothing if already initialized
       IF ( .not. initialized ) then
          call TBrook_Initialize(iStatus)
          if ( iStatus /= 0 ) return
       END IF

       call TBrook_Destroy(BaseBrook, iStatus=iStatus)
       if ( iStatus /= 0 ) return
       ! Set up all the base brook
       if ( present(FILE) ) then
          no_file_set = .false.
          CALL TBrook_Set(BaseBrook,                               &
               FILE=TRIM(FILE),                         &
               FORM="XML",                              &
               UNIT=myUnit,                             &
               XMLRoot="TruchasData",                   &
               iStatus=iStatus)
          CALL TBrook_Write(b=BaseBrook,             &
               Variable=TBrookVersion,  &
               XMLName="TBrookVersion", &
               iStatus=iStatus)
          BaseBrook%Next => NULL()

          BaseBrookAscii = BaseBrook
          call TBrook_set(BaseBrookAscii, FORM="ascii", iStatus=iStatus)

          ! Now copy the structure into the rest of the brooks so that
          ! we can write in multiple forms and extract what we need from
          ! the different forms

       end if
       RETURN
     END SUBROUTINE TBROOK_SETBASEFILE

     SUBROUTINE TB_ID_WXMLC_STRING(TB_ID, Comment, iStatus)
       integer, intent(in) :: TB_ID
       character(len=*), intent(in) :: Comment
       integer, intent(inout) :: iStatus

       call Brook_WriteXMLComment(TBrook_Get_Brook(TB_ID), Comment, iStatus)
       return
     END SUBROUTINE TB_ID_WXMLC_STRING

     SUBROUTINE TB_ID_WXMLC_STRINGARRAY(TB_ID, Comment, iStatus)
       integer, intent(in) :: TB_ID
       character(len=*), dimension(:), intent(in) :: Comment
       integer, intent(inout) :: iStatus

       call Brook_WriteXMLComment(TBrook_Get_Brook(TB_ID), Comment, iStatus)
       return
     END SUBROUTINE TB_ID_WXMLC_STRINGARRAY


     !The interface to the brook xml tag writers
     SUBROUTINE TBROOK_WRITEXMLTAG(B,             & ! Brook, target
          XMLTag,        & ! Name of tag to be opened
          XMLAttributes, & ! Attributes for tag to be opened
          XMLStringData, & ! String with data to be written inline
          Scope,         & ! Is this a local stream
          iStatus        & ! Did an error occur?
          )
       use brook_module, only: Brook_WriteXMLTag
       TYPE(BROOK),        TARGET, INTENT(IN)  :: B
       CHARACTER(LEN=*),           INTENT(IN)  :: XMLTag
       CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: XMLAttributes
       CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: XMLStringData
       INTEGER,          OPTIONAL, INTENT(IN)  :: SCOPE
       INTEGER,                    INTENT(OUT) :: iStatus

       iStatus = 0
       if ( I_SHALL_WRITE(SCOPE) ) then
          call Brook_WriteXMLTag(B=B, &
                                 XMLTag=XMLTag, &
                                 XMLAttributes=XMLAttributes, &
                                 XMLStringData=XMLStringData, &
                                 iStatus=iStatus)
       end if

       RETURN

     END SUBROUTINE TBROOK_WRITEXMLTAG

     SUBROUTINE tbrook_OpenXMLTag(B, XMLTag, XMLAttributes, Scope, iStatus)
       use brook_module, only: Brook_OpenXMLTag
       TYPE(BROOK),        TARGET, INTENT(IN) :: B
       CHARACTER(LEN=*),           INTENT(IN) :: XMLTag
       CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: XMLAttributes
       INTEGER,          OPTIONAL, INTENT(IN)  :: SCOPE
       INTEGER,                    INTENT(INOUT) :: iStatus

       ! Local Variables
       iStatus = 0

       if ( I_SHALL_WRITE(SCOPE) ) then
          call brook_OpenXMLTag(B=B, XMLTag=XMLTag, XMLAttributes=XMLAttributes, iStatus=iStatus)
       end if
       RETURN
     END SUBROUTINE tbrook_OpenXMLTag

     SUBROUTINE tbrook_CloseXMLTag(B, XMLTag, Scope, iStatus)
       use brook_module, only: Brook_CloseXMLTag
       TYPE(BROOK),        TARGET, INTENT(IN) :: B
       CHARACTER(LEN=*),           INTENT(IN) :: XMLTag
       INTEGER,          OPTIONAL, INTENT(IN) :: SCOPE
       INTEGER                                :: iStatus

       ! Local Variables
       iStatus = 0

       if ( I_SHALL_WRITE(SCOPE) ) then
          call brook_CloseXMLTag(B, XMLTag, iStatus)
       end if
       RETURN
     END SUBROUTINE tbrook_CloseXMLTag


  !=========================================
  ! SUBROUTINES Brook_WRITE_LOGICAL_*
  ! 
  ! Purpose: Write logical variables to output brooks
#define _TYPE_                  Logical
#define _DIMENSION_             0
#define _FUNCTION_NAME_         TB_WRITE_LOGICAL_R0
#define _BROOK_FUNC_            B_WRITE_LOGICAL_R0
#define _COLLATE_FUNC_          TB_COLLATE_LOGICAL
#define _FUNCTION_NAME_B_       TBB_WRITE_LOGICAL_R0
#include "tbrook_include.fpp"

#define _TYPE_                  Logical
#define _DIMENSION_             1
#define _FUNCTION_NAME_         TB_WRITE_LOGICAL_R1
#define _BROOK_FUNC_            B_WRITE_LOGICAL_R1
#define _COLLATE_FUNC_          TB_COLLATE_LOGICAL
#define _FUNCTION_NAME_B_       TBB_WRITE_LOGICAL_R1
#include "tbrook_include.fpp"

#define _TYPE_                  Logical
#define _DIMENSION_             2
#define _FUNCTION_NAME_         TB_WRITE_LOGICAL_R2
#define _BROOK_FUNC_            B_WRITE_LOGICAL_R2
#define _COLLATE_FUNC_          TB_COLLATE_LOGICAL
#define _FUNCTION_NAME_B_       TBB_WRITE_LOGICAL_R2
#include "tbrook_include.fpp"

#define _TYPE_                  Logical
#define _DIMENSION_             3
#define _FUNCTION_NAME_         TB_WRITE_LOGICAL_R3
#define _BROOK_FUNC_            B_WRITE_LOGICAL_R3
#define _COLLATE_FUNC_          TB_COLLATE_LOGICAL
#define _FUNCTION_NAME_B_       TBB_WRITE_LOGICAL_R3
#include "tbrook_include.fpp"

#define _TYPE_                  Logical
#define _DIMENSION_             4
#define _FUNCTION_NAME_         TB_WRITE_LOGICAL_R4
#define _BROOK_FUNC_            B_WRITE_LOGICAL_R4
#define _COLLATE_FUNC_          TB_COLLATE_LOGICAL
#define _FUNCTION_NAME_B_       TBB_WRITE_LOGICAL_R4
#include "tbrook_include.fpp"

#define _TYPE_                  Logical
#define _DIMENSION_             5
#define _FUNCTION_NAME_         TB_WRITE_LOGICAL_R5
#define _BROOK_FUNC_            B_WRITE_LOGICAL_R5
#define _COLLATE_FUNC_          TB_COLLATE_LOGICAL
#define _FUNCTION_NAME_B_       TBB_WRITE_LOGICAL_R5
#include "tbrook_include.fpp"

  ! END SUBROUTINES TB_WRITE_LOGICAL_*
  !=========================================

  !=========================================
  ! SUBROUTINES TB_WRITE_INTEGER_*
  ! Purpose: Write integer variables to output brooks
#define _TYPE_                  Integer
#define _DIMENSION_             0
#define _FUNCTION_NAME_         TB_WRITE_INTEGER_R0
#define _BROOK_FUNC_            B_WRITE_INTEGER_R0
#define _COLLATE_FUNC_          TB_COLLATE_INTEGER
#define _FUNCTION_NAME_B_       TBB_WRITE_INTEGER_R0
#include "tbrook_include.fpp"

#define _TYPE_                  Integer
#define _DIMENSION_             1
#define _FUNCTION_NAME_         TB_WRITE_INTEGER_R1
#define _BROOK_FUNC_            B_WRITE_INTEGER_R1
#define _COLLATE_FUNC_          TB_COLLATE_INTEGER
#define _FUNCTION_NAME_B_       TBB_WRITE_INTEGER_R1
#include "tbrook_include.fpp"

#define _TYPE_                  Integer
#define _DIMENSION_             2
#define _FUNCTION_NAME_         TB_WRITE_INTEGER_R2
#define _BROOK_FUNC_            B_WRITE_INTEGER_R2
#define _COLLATE_FUNC_          TB_COLLATE_INTEGER
#define _FUNCTION_NAME_B_       TBB_WRITE_INTEGER_R2
#include "tbrook_include.fpp"

#define _TYPE_                  Integer
#define _DIMENSION_             3
#define _FUNCTION_NAME_         TB_WRITE_INTEGER_R3
#define _BROOK_FUNC_            B_WRITE_INTEGER_R3
#define _COLLATE_FUNC_          TB_COLLATE_INTEGER
#define _FUNCTION_NAME_B_       TBB_WRITE_INTEGER_R3
#include "tbrook_include.fpp"

#define _TYPE_                  Integer
#define _DIMENSION_             4
#define _FUNCTION_NAME_         TB_WRITE_INTEGER_R4
#define _BROOK_FUNC_            B_WRITE_INTEGER_R4
#define _COLLATE_FUNC_          TB_COLLATE_INTEGER
#define _FUNCTION_NAME_B_       TBB_WRITE_INTEGER_R4
#include "tbrook_include.fpp"

#define _TYPE_                  Integer
#define _DIMENSION_             5
#define _FUNCTION_NAME_         TB_WRITE_INTEGER_R5
#define _BROOK_FUNC_            B_WRITE_INTEGER_R5
#define _COLLATE_FUNC_          TB_COLLATE_INTEGER
#define _FUNCTION_NAME_B_       TBB_WRITE_INTEGER_R5
#include "tbrook_include.fpp"
  ! END SUBROUTINES TB_WRITE_INTEGER_*
  !=========================================

  !=========================================
  ! SUBROUTINES TB_WRITE_DOUBLE_*
  ! Purpose: Write double precision variables to output brooks
#define _TYPE_                  double precision
#define _DIMENSION_             0
#define _FUNCTION_NAME_         TB_WRITE_DOUBLE_R0
#define _BROOK_FUNC_            B_WRITE_DOUBLE_R0
#define _COLLATE_FUNC_          TB_COLLATE_DOUBLE
#define _FUNCTION_NAME_B_       TBB_WRITE_DOUBLE_R0
#include "tbrook_include.fpp"

#define _TYPE_                  double precision
#define _DIMENSION_             1
#define _FUNCTION_NAME_         TB_WRITE_DOUBLE_R1
#define _BROOK_FUNC_            B_WRITE_DOUBLE_R1
#define _COLLATE_FUNC_          TB_COLLATE_DOUBLE
#define _FUNCTION_NAME_B_       TBB_WRITE_DOUBLE_R1
#include "tbrook_include.fpp"

#define _TYPE_                  double precision
#define _DIMENSION_             2
#define _FUNCTION_NAME_         TB_WRITE_DOUBLE_R2
#define _BROOK_FUNC_            B_WRITE_DOUBLE_R2
#define _COLLATE_FUNC_          TB_COLLATE_DOUBLE
#define _FUNCTION_NAME_B_       TBB_WRITE_DOUBLE_R2
#include "tbrook_include.fpp"

#define _TYPE_                  double precision
#define _DIMENSION_             3
#define _FUNCTION_NAME_         TB_WRITE_DOUBLE_R3
#define _BROOK_FUNC_            B_WRITE_DOUBLE_R3
#define _COLLATE_FUNC_          TB_COLLATE_DOUBLE
#define _FUNCTION_NAME_B_       TBB_WRITE_DOUBLE_R3
#include "tbrook_include.fpp"

#define _TYPE_                  double precision
#define _DIMENSION_             4
#define _FUNCTION_NAME_         TB_WRITE_DOUBLE_R4
#define _BROOK_FUNC_            B_WRITE_DOUBLE_R4
#define _COLLATE_FUNC_          TB_COLLATE_DOUBLE
#define _FUNCTION_NAME_B_       TBB_WRITE_DOUBLE_R4
#include "tbrook_include.fpp"

#define _TYPE_                  double precision
#define _DIMENSION_             5
#define _FUNCTION_NAME_         TB_WRITE_DOUBLE_R5
#define _BROOK_FUNC_            B_WRITE_DOUBLE_R5
#define _COLLATE_FUNC_          TB_COLLATE_DOUBLE
#define _FUNCTION_NAME_B_       TBB_WRITE_DOUBLE_R5
#include "tbrook_include.fpp"
  ! END SUBROUTINES TB_WRITE_DOUBLE_*
  !=========================================
  !=========================================
  ! SUBROUTINES TB_WRITE_FLOAT_*
  ! Purpose: Write real variables to output brooks
#define _TYPE_                  real
#define _DIMENSION_             0
#define _FUNCTION_NAME_         TB_WRITE_FLOAT_R0
#define _BROOK_FUNC_            B_WRITE_FLOAT_R0
#define _COLLATE_FUNC_          TB_COLLATE_FLOAT
#define _FUNCTION_NAME_B_       TBB_WRITE_FLOAT_R0
#include "tbrook_include.fpp"

#define _TYPE_                  real
#define _DIMENSION_             1
#define _FUNCTION_NAME_         TB_WRITE_FLOAT_R1
#define _BROOK_FUNC_            B_WRITE_FLOAT_R1
#define _COLLATE_FUNC_          TB_COLLATE_FLOAT
#define _FUNCTION_NAME_B_       TBB_WRITE_FLOAT_R1
#include "tbrook_include.fpp"

#define _TYPE_                  real
#define _DIMENSION_             2
#define _FUNCTION_NAME_         TB_WRITE_FLOAT_R2
#define _BROOK_FUNC_            B_WRITE_FLOAT_R2
#define _COLLATE_FUNC_          TB_COLLATE_FLOAT
#define _FUNCTION_NAME_B_       TBB_WRITE_FLOAT_R2
#include "tbrook_include.fpp"

#define _TYPE_                  real
#define _DIMENSION_             3
#define _FUNCTION_NAME_         TB_WRITE_FLOAT_R3
#define _BROOK_FUNC_            B_WRITE_FLOAT_R3
#define _COLLATE_FUNC_          TB_COLLATE_FLOAT
#define _FUNCTION_NAME_B_       TBB_WRITE_FLOAT_R3
#include "tbrook_include.fpp"

#define _TYPE_                  real
#define _DIMENSION_             4
#define _FUNCTION_NAME_         TB_WRITE_FLOAT_R4
#define _BROOK_FUNC_            B_WRITE_FLOAT_R4
#define _COLLATE_FUNC_          TB_COLLATE_FLOAT
#define _FUNCTION_NAME_B_       TBB_WRITE_FLOAT_R4
#include "tbrook_include.fpp"

#define _TYPE_                  real
#define _DIMENSION_             5
#define _FUNCTION_NAME_         TB_WRITE_FLOAT_R5
#define _BROOK_FUNC_            B_WRITE_FLOAT_R5
#define _COLLATE_FUNC_          TB_COLLATE_FLOAT
#define _FUNCTION_NAME_B_       TBB_WRITE_FLOAT_R5
#include "tbrook_include.fpp"
  ! END SUBROUTINES TB_WRITE_FLOAT_*
  !=========================================

  !=========================================
  ! SUBROUTINES TB_WRITE_STRING_*
  ! Purpose: Write string variables to output brooks
  ! Note only single arrays are allowed
#define _TYPE_                  Character(len=*)
#define _DIMENSION_             0
#define _FUNCTION_NAME_         TB_WRITE_STRING_R0
#define _BROOK_FUNC_            B_WRITE_STRING_R0
#define _COLLATE_FUNC_          TB_COLLATE_STRING
#define _FUNCTION_NAME_B_       TBB_WRITE_STRING_R0
#define _IS_CHARACTER_TYPE_
#include "tbrook_include.fpp"

#define _TYPE_                  Character(len=*)
#define _DIMENSION_             1
#define _FUNCTION_NAME_         TB_WRITE_STRING_R1
#define _BROOK_FUNC_            B_WRITE_STRING_R1
#define _COLLATE_FUNC_          TB_COLLATE_STRING
#define _FUNCTION_NAME_B_       TBB_WRITE_STRING_R1
#define _IS_CHARACTER_TYPE_
#include "tbrook_include.fpp"

  ! END SUBROUTINES TB_WRITE_CHARACTER_*
  !=========================================

#define _TBU_FILEENTRY_FUNC_ TBU_FILEENTRY_i0  
#define _DIMENSION_       0
#define _TYPE_            integer
#define _TYPE_CODE_       'i'
#include "tbu_include.fpp"

#define _TBU_FILEENTRY_FUNC_ TBU_FILEENTRY_i1  
#define _DIMENSION_       1
#define _TYPE_            integer
#define _TYPE_CODE_       'i'
#include "tbu_include.fpp"

#define _TBU_FILEENTRY_FUNC_ TBU_FILEENTRY_i2  
#define _DIMENSION_       2
#define _TYPE_            integer
#define _TYPE_CODE_       'i'
#include "tbu_include.fpp"

#define _TBU_FILEENTRY_FUNC_ TBU_FILEENTRY_i3
#define _DIMENSION_       3
#define _TYPE_            integer
#define _TYPE_CODE_       'i'
#include "tbu_include.fpp"

#define _TBU_FILEENTRY_FUNC_ TBU_FILEENTRY_i4
#define _DIMENSION_       4
#define _TYPE_            integer
#define _TYPE_CODE_       'i'
#include "tbu_include.fpp"

#define _TBU_FILEENTRY_FUNC_ TBU_FILEENTRY_i5
#define _DIMENSION_       5
#define _TYPE_            integer
#define _TYPE_CODE_       'i'
#include "tbu_include.fpp"

#define _TBU_FILEENTRY_FUNC_ TBU_FILEENTRY_r0  
#define _DIMENSION_       0
#define _TYPE_            real
#define _TYPE_CODE_       'r'
#include "tbu_include.fpp"

#define _TBU_FILEENTRY_FUNC_ TBU_FILEENTRY_r1  
#define _DIMENSION_       1
#define _TYPE_            real
#define _TYPE_CODE_       'r'
#include "tbu_include.fpp"

#define _TBU_FILEENTRY_FUNC_ TBU_FILEENTRY_r2  
#define _DIMENSION_       2
#define _TYPE_            real
#define _TYPE_CODE_       'r'
#include "tbu_include.fpp"

#define _TBU_FILEENTRY_FUNC_ TBU_FILEENTRY_r3
#define _DIMENSION_       3
#define _TYPE_            real
#define _TYPE_CODE_       'r'
#include "tbu_include.fpp"

#define _TBU_FILEENTRY_FUNC_ TBU_FILEENTRY_r4
#define _DIMENSION_       4
#define _TYPE_            real
#define _TYPE_CODE_       'r'
#include "tbu_include.fpp"

#define _TBU_FILEENTRY_FUNC_ TBU_FILEENTRY_r5
#define _DIMENSION_       5
#define _TYPE_            real
#define _TYPE_CODE_       'r'
#include "tbu_include.fpp"

#define _TBU_FILEENTRY_FUNC_ TBU_FILEENTRY_l0  
#define _DIMENSION_       0
#define _TYPE_            logical
#define _TYPE_CODE_       'l'
#include "tbu_include.fpp"

#define _TBU_FILEENTRY_FUNC_ TBU_FILEENTRY_l1  
#define _DIMENSION_       1
#define _TYPE_            logical
#define _TYPE_CODE_       'l'
#include "tbu_include.fpp"

#define _TBU_FILEENTRY_FUNC_ TBU_FILEENTRY_l2  
#define _DIMENSION_       2
#define _TYPE_            logical
#define _TYPE_CODE_       'l'
#include "tbu_include.fpp"

#define _TBU_FILEENTRY_FUNC_ TBU_FILEENTRY_l3
#define _DIMENSION_       3
#define _TYPE_            logical
#define _TYPE_CODE_       'l'
#include "tbu_include.fpp"

#define _TBU_FILEENTRY_FUNC_ TBU_FILEENTRY_l4
#define _DIMENSION_       4
#define _TYPE_            logical
#define _TYPE_CODE_       'l'
#include "tbu_include.fpp"

#define _TBU_FILEENTRY_FUNC_ TBU_FILEENTRY_l5
#define _DIMENSION_       5
#define _TYPE_            logical
#define _TYPE_CODE_       'l'
#include "tbu_include.fpp"

#define _TBU_FILEENTRY_FUNC_ TBU_FILEENTRY_d0  
#define _DIMENSION_       0
#define _TYPE_            double precision
#define _TYPE_CODE_       'd'
#include "tbu_include.fpp"

#define _TBU_FILEENTRY_FUNC_ TBU_FILEENTRY_d1  
#define _DIMENSION_       1
#define _TYPE_            double precision
#define _TYPE_CODE_       'd'
#include "tbu_include.fpp"

#define _TBU_FILEENTRY_FUNC_ TBU_FILEENTRY_d2  
#define _DIMENSION_       2
#define _TYPE_            double precision
#define _TYPE_CODE_       'd'
#include "tbu_include.fpp"

#define _TBU_FILEENTRY_FUNC_ TBU_FILEENTRY_d3
#define _DIMENSION_       3
#define _TYPE_            double precision
#define _TYPE_CODE_       'd'
#include "tbu_include.fpp"

#define _TBU_FILEENTRY_FUNC_ TBU_FILEENTRY_d4
#define _DIMENSION_       4
#define _TYPE_            double precision
#define _TYPE_CODE_       'd'
#include "tbu_include.fpp"

#define _TBU_FILEENTRY_FUNC_ TBU_FILEENTRY_d5
#define _DIMENSION_       5
#define _TYPE_            double precision
#define _TYPE_CODE_       'd'
#include "tbu_include.fpp"


#define _TBU_FILEENTRY_FUNC_ TBU_FILEENTRY_c0  
#define _DIMENSION_       0
#define _TYPE_            character(len=*)
#define _TYPE_CODE_       'c'
#include "tbu_include.fpp"

#define _TBU_FILEENTRY_FUNC_ TBU_FILEENTRY_c1  
#define _DIMENSION_       1
#define _TYPE_            character(len=*)
#define _TYPE_CODE_       'c'
#include "tbu_include.fpp"

END MODULE tbrook_module
