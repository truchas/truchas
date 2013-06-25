MODULE PGSLib_GS_Setup_MODULE
  use pgslib_Instrument
  use PGSlib_Stats
  use pgslib_timing_module
  IMPLICIT NONE 
  SAVE
  PRIVATE
  PUBLIC :: PGSLib_Setup_Trace, PGSLib_Setup_Basic_Trace


  ! $Id: pgslib_gs_setup_module.F,v 1.2 2000/11/02 21:05:55 ferrell Exp $

  ! Provides Setup_Trace, which establishes the Requestor<->Node communication pattern, or
  ! whatever communication is indicated by the inputs.  
  ! The terminology (in the case of elementss and nodes):
  !  Requestor = Elements, since requests an operation (gather or scatter)
  !  Complior  = Nodes, since complies with request.


  INTERFACE PGSLib_Setup_Trace
     MODULE PROCEDURE PGSLib_Setup_Trace_1D
     MODULE PROCEDURE PGSLib_Setup_Trace_2D
  END INTERFACE

  INTERFACE PGSLib_Setup_Supplement
     MODULE PROCEDURE PGSLib_Setup_Supplement_F
  END INTERFACE

  INTERFACE PGSLib_Setup_Duplicate
     MODULE PROCEDURE PGSLib_Setup_Duplicate_F
  END INTERFACE

! Interfaces to the C routines which support this module.

  INTERFACE PGSLib_Setup_N_Duplicate
     subroutine PGSLib_Setup_N_Duplicate_C(N_Duplicate, GS_Trace)
       implicit none
       integer :: N_Duplicate
       integer, DIMENSION(*) :: GS_Trace
     END subroutine PGSLib_Setup_N_Duplicate_C
  END INTERFACE

  INTERFACE PGSLib_Setup_Duplicate_Buffer
     subroutine pgslib_setup_duplicate_buffer_C(GS_Trace)
       implicit none
       integer, DIMENSION(*)  :: GS_Trace
     END subroutine pgslib_setup_duplicate_buffer_C
  END INTERFACE
     
  INTERFACE PGSLib_Prep_Supplement
     subroutine PGSLib_Prep_Supplement_C(N_Supplement, PEs, GS_Trace)
       implicit none
       integer :: N_Supplement
       integer, dimension(*) :: PEs
       integer, DIMENSION(*) :: GS_Trace
     END subroutine PGSLib_Prep_Supplement_C
  END INTERFACE

  !======================================================================
CONTAINS
  
  !====================================================================
  !          PGSLib_Setup_Trace
  ! PURPOSE:
  !         Return a pointer to a valid trace for the input communication
  !         pattern.  Also, renumber input index so that it points
  !         to either local data, or to supplement communication buffer.
  !         IF PE_ARRAY is present and allocated, we assume it contains
  !         the destination PE for each index.
  !         If PE_ARRAY is present but not allocated, fill it up and return
  !         with the PE number for each index.
  !         IF PE_ARRAY is not present, or present but not allocated, then
  !         SIZE_OF_DEST must be present.
  !
  !========================================================================
  function PGSLib_Setup_Trace_1D(INDEX,                &
                             SIZE_OF_DEST,            &
                             PE_ARRAY,                &
                             MASK)                    &
           RESULT(Trace)
    
    USE PGSLib_Type_MODULE
    USE PGSLib_Globals_MODULE
    USE PGSLib_Utility_MODULE
    USE PGSLib_Index_MODULE
    USE PGSLIB_GS_Util_MODULE
    use pgslib_Instrument

    IMPLICIT NONE 
    integer (PGSLib_Int_Type), INTENT(INOUT),      &       
                               dimension(:)        :: INDEX

    integer (PGSLib_Int_Type), INTENT(IN   ),      &
                               OPTIONAL            :: SIZE_OF_DEST

    integer (PGSLib_Int_Type), OPTIONAL,           &
                               POINTER,            &
                               dimension(:  )      :: PE_ARRAY
    logical (PGSLib_Log_Type), INTENT(IN   ),      &                                         
                               OPTIONAL,           &
                               dimension(:)        :: MASK

    type   (PGSLib_GS_Trace),  POINTER             :: Trace

    ! Local variables

    integer (PGSLib_Int_Type), POINTER,            &
         &                     dimension(:,:)      :: INDEX_2D

    integer (PGSLib_Int_Type), POINTER,            &
         &                     dimension(:,:)      :: PE_ARRAY_2D
    logical (PGSLib_Log_Type), POINTER,            &
                               dimension(:,:)      :: MASK_2D
    logical :: NEW_PE_ARRAY

    ! Nullify all local pointers, to put them into known state
    NULLIFY(INDEX_2D)
    NULLIFY(PE_ARRAY_2D)
    NULLIFY(MASK_2D)

    ! Move data into 2D arrays, then call 2D function
    ALLOCATE(INDEX_2D(1, SIZE(INDEX,1)))
    INDEX_2D(1,:) = INDEX
    if (PRESENT(PE_ARRAY)) then
       NEW_PE_ARRAY = .NOT. ASSOCIATED(PE_ARRAY)
       IF (NEW_PE_ARRAY) then
          NULLIFY(PE_ARRAY_2D)
       else
          ALLOCATE(PE_ARRAY_2D(1, SIZE(PE_ARRAY,1)))
          PE_ARRAY_2D(1,:) = PE_ARRAY
       END IF
    ELSE
       NEW_PE_ARRAY = .FALSE.
       NULLIFY(PE_ARRAY_2D)
    END if
    IF (PRESENT(MASK)) then
       ALLOCATE(MASK_2D(1, SIZE(MASK,1)))
       MASK_2D(1,:) = MASK
    END IF
    
    IF (PRESENT(PE_ARRAY) .AND. PRESENT(MASK)) then
       Trace => PGSLib_Setup_Trace(INDEX_2D, SIZE_OF_DEST, PE_ARRAY=PE_ARRAY_2D, MASK=MASK_2D)
    END IF

    IF (PRESENT(PE_ARRAY) .AND. .NOT. PRESENT(MASK)) then
       Trace => PGSLib_Setup_Trace(INDEX_2D, SIZE_OF_DEST, PE_ARRAY=PE_ARRAY_2D)
    END IF

    IF (.NOT. PRESENT(PE_ARRAY) .AND. PRESENT(MASK)) then
       Trace => PGSLib_Setup_Trace(INDEX_2D, SIZE_OF_DEST, MASK=MASK_2D)
    END IF

    IF (.NOT. PRESENT(PE_ARRAY) .AND. .NOT. PRESENT(MASK)) then
       Trace => PGSLib_Setup_Trace(INDEX_2D, SIZE_OF_DEST)
    END IF
    
    INDEX = INDEX_2D(1,:)
    DEALLOCATE(INDEX_2D)
    if (NEW_PE_ARRAY) then
       ALLOCATE(PE_ARRAY(SIZE(INDEX)))
       PE_ARRAY = PE_ARRAY_2D(1,:)
    end if
    if (ASSOCIATED(PE_ARRAY_2D)) DEALLOCATE(PE_ARRAY_2D)
    if (ASSOCIATED(MASK_2D))     DEALLOCATE(MASK_2D)

  end function PGSLib_Setup_Trace_1D



  function PGSLib_Setup_Trace_2D(INDEX,                &
                             SIZE_OF_DEST,            &
                             PE_ARRAY,                &
                             MASK)                    &
           RESULT(Trace)
    
    USE PGSLib_Type_MODULE
    USE PGSLib_Globals_MODULE
    USE PGSLib_Utility_MODULE
    USE PGSLib_Index_MODULE
    USE PGSLIB_GS_Util_MODULE

    IMPLICIT NONE 
    integer (PGSLib_Int_Type), INTENT(INOUT),      &                                         
         &                     dimension(:, :)     :: INDEX

    integer (PGSLib_Int_Type), INTENT(IN   ),      &
                               OPTIONAL            :: SIZE_OF_DEST

    integer (PGSLib_Int_Type), OPTIONAL,           &
         &                     POINTER,            &
         &                     dimension(:,:)      :: PE_ARRAY
    logical (PGSLib_Log_Type), INTENT(IN   ),      &                                         
                               OPTIONAL,           &
                               dimension(:, :)     :: MASK

    type   (PGSLib_GS_Trace),  POINTER             :: Trace

    ! Local variables
    INTEGER (PGSLib_Int_TYPE) :: NRequestors, CompliorsPerRequestor
    logical (PGSLib_Log_Type) :: Mask_Present, User_Data_Map
    integer (PGSLib_Int_Type), POINTER, dimension(:,:) :: Local_PE_ARRAY
    type (PGSLib_GID),         POINTER             :: GID
    integer (PGSLib_Int_Type), POINTER,            &
                               dimension(:)        :: Sup_Global_I
    integer (PGSLib_Int_Type), POINTER,            &
                               dimension(:)        :: Sup_Global_PE
    real (PGSLib_Single_Type)  :: Top_of_Setup_Time, After_Supp_Time


#ifdef USE_TIMERS_1
! Start the timer
    call Enter_Routine(Setup_Trace_STATISTICS())
    Top_of_Setup_Time = PGSLib_Read_Elapsed_Time(Setup_Trace_STATISTICS())
#endif

    ! NULLIFY all local pointers, to get them into known state
    NULLIFY(Local_PE_ARRAY)
    NULLIFY(GID)
    NULLIFY(Sup_Global_I)
    NULLIFY(Sup_Global_PE)

    ! If we have a mask, that is a permanent mask, so keep a record of that
    Mask_Present          = PRESENT(MASK)

    ! Next we make local INDEX from global INDEX.  At the same time
    ! we figure out which compliors we need (request) from other PE''s. That information
    ! is used to establish the Supplement buffers in the trace.

    ! If the user supplied a PE_ARRAY that array might already contain the
    ! mapping of INDEX to PE''s.
    if (PRESENT(PE_ARRAY)) then
       Local_PE_ARRAY => PE_ARRAY
    else
       NULLIFY(Local_PE_ARRAY)
    end if

    ! If we need to construct the PE_ARRAY, then we need SIZE_of_DEST.
    if (.NOT. PRESENT(PE_ARRAY)) then
       IF (.NOT. PRESENT(SIZE_OF_DEST)) then
          call pgslib_fatal_error("If PE_ARRAY is not supplied SIZE_OF_DEST " // &
                                   "must be supplied in PGSLib_Setup_Trace")
       end IF
    else
       IF (.NOT. ASSOCIATED(PE_ARRAY)) then
          IF (.NOT. PRESENT(SIZE_OF_DEST)) then
             call pgslib_fatal_error("If PE_ARRAY is not initialized SIZE_OF_DEST " // &
                                     "must be supplied in PGSLib_Setup_Trace")
          end IF
       end IF
    end if

    NULLIFY(GID)
    NULLIFY(Sup_Global_I)
    NULLIFY(Sup_Global_PE)
    call PGSLib_Setup_Supplement(INDEX, Local_PE_ARRAY, GID,      &
                                 Sup_Global_I, Sup_Global_PE, &
                                 SIZE_OF_DEST, MASK)

    IF (.NOT. PRESENT(PE_ARRAY) ) then
       DEALLOCATE(Local_PE_ARRAY)
    else
       if (.NOT. ASSOCIATED(PE_ARRAY)) PE_Array => Local_PE_ARRAY
       NULLIFY(Local_PE_ARRAY)
    end IF

#ifdef USE_TIMERS_1
    After_Supp_Time = PGSLib_Read_Elapsed_Time(Setup_Trace_STATISTICS())
#endif
    ! Establish the trace, which means communicating what
    ! are the supplements and what are the duplicates.
    Trace => PGSLib_Setup_Basic_Trace(Sup_Global_PE)

    ! If we had a mask, need to remember that
    Trace%Present_PMask = Mask_Present

#ifdef USE_TIMERS_1
    ! Start the Trace timer.
    ! Since we couldn't start it before we established the trace, we didn't get the
    ! time spent in setup_supplement.  So, add that on now.
    call PGSLib_Increment_Timer(Trace%Setup_Timer, (After_Supp_Time - Top_of_Setup_Time))
    call pgslib_start_timer(Trace%Setup_Timer)
#endif
    ! Now put some of the data into the trace fields
    ALLOCATE(Trace%Supplement_Global_Indices(PGSLib_Size_Of_Sup(Trace)))
    ALLOCATE(Trace%Supplement_Global_PEs(PGSLib_Size_Of_Sup(Trace)))
    ALLOCATE(Trace%Supplement_I_Data(PGSLib_Size_Of_Sup(Trace)))
    ALLOCATE(Trace%Supplement_R_Data(PGSLib_Size_Of_Sup(Trace)))
    ALLOCATE(Trace%Supplement_D_Data(PGSLib_Size_Of_Sup(Trace)))
    ALLOCATE(Trace%Supplement_L_Data(PGSLib_Size_Of_Sup(Trace)))

    Trace%Supplement_Global_Indices = Sup_Global_I
    Trace%Supplement_Global_PEs     = Sup_Global_PE
    DEALLOCATE(Sup_Global_I, Sup_Global_PE)

    ! Final Step: Communication to resolve who needs some of my Compliors.
    ! Since the trace has all the information, we need only pass in the trace.
    ! This is based on the assumption that we can resolve INDEX with the GID''s.
    ! The communication buffer is the Duplicate buffer, because any item
    ! in it is a duplicate of an item owned by this PE (although the values may
    ! differ in the middle of some communication operations).

    Call PGSLib_Setup_Duplicate(Trace, GID)

    call pgslib_release_gid(GID)
#ifdef USE_TIMERS_1
    ! Stop the timers
    call Exit_Routine(Setup_Trace_STATISTICS())
    call PGSLib_Stop_Timer(Trace%Setup_Timer)
#endif

  END Function PGSLib_Setup_Trace_2D


!======================================================================
!          PGSLib_Setup_Supplement(Requestor, Requestor_PE, Requestor_PE_Flag, Trace)
! PURPOSE:
!          Resolve which items (compliors) we need to get from other PE''s.
!          These go in a buffer called the Supplement (because these
!          items supplement the items owned by this PE).
!          This routine also ALLOCATES the Supplement in the Trace
!======================================================================
  subroutine PGSLib_Setup_Supplement_F(INDEX,        &
                                       PE_ARRAY,       &
                                       GID,            &
                                       Sup_Global_I,   &
                                       Sup_Global_PE,  &
                                       SIZE_OF_DEST,   &
                                       MASK)
    USE PGSLib_Type_MODULE
    USE PGSLib_Globals_MODULE
    USE PGSLib_Utility_MODULE
    USE PGSLib_Index_MODULE
    USE PGSLIB_GS_Util_MODULE
    implicit none 
    
    integer (PGSLib_Int_Type), INTENT(INOUT),      &                                         
         &                     dimension(:, :)     :: INDEX

    integer (PGSLib_Int_Type), POINTER,            &
                               dimension(:,:)      :: PE_ARRAY

    type (PGSLib_GID), POINTER                     :: GID
    integer (PGSLib_Int_Type), POINTER,            &
                               dimension(:)         :: Sup_Global_I
    integer (PGSLib_Int_Type), POINTER,            &
                               dimension(:)         :: Sup_Global_PE

    integer (PGSLib_Int_Type), INTENT(IN   ),      &
                               OPTIONAL            :: SIZE_OF_DEST

    logical (PGSLib_Log_Type), INTENT(IN   ),      &                                         
                               OPTIONAL,           &
                               dimension(:, :)     :: MASK

    ! Local variables
    integer (PGSLib_Int_Type) :: R1, R2, N_Supplement
    INTEGER (PGSLib_Int_TYPE) :: thisPE, PE, ierror
    logical (PGSLib_Log_Type) :: Mask_Present, NEW_PE_ARRAY
    type (PGSLib_Access_Table), POINTER :: Table

    ! Nullify all local pointers to put them into known state
    NULLIFY(Table)

#ifdef USE_TIMERS_1
    call Enter_Routine(Setup_Supplement_STATISTICS())
#endif    

    ierror = 0
    Mask_Present = PRESENT(MASK)

    ! We need to build a table of which Compliors are needed from other PEs
    Table => PGSLib_Init_offPE_Access_Table()

    ! If we have a PE array with stuff in it, we use that, otherwise we need
    ! to fill up the PE array.
    NEW_PE_ARRAY = .NOT. ASSOCIATED(PE_ARRAY)
    if (NEW_PE_ARRAY) then
       if (.NOT. PRESENT(SIZE_OF_DEST)) then
          call pgslib_fatal_error("If PE_ARRAY is not initialized SIZE_OF_DEST " // &
                                  "must be supplied in PGSLib_Setup_Supplement")
       end if
       ALLOCATE(PE_ARRAY(SIZE(INDEX,1), SIZE(INDEX,2)))
       GID => PGSLib_Setup_GID(SIZE_OF_DEST)
    end if

    thisPE = PGSLib_Inquire_thisPE_Actual() 
    
#ifdef USE_TIMERS_1
    call Enter_ROUTINE(Supplement_G_TO_L_Statistics())
#endif
    DO R2 = 1, SIZE(INDEX,2)
       Do R1 = 1, SIZE(INDEX,1)

          ! If Mask, then if this requestor is inactive, skip to next one
          if (Mask_Present) then
             if (.NOT. MASK(R1, R2) ) CYCLE
          endif

          ! Only do this if no mask, or if Mask(CmplPtr, Req) == TRUE
          if (NEW_PE_ARRAY) then
             ! FORTRAN Counts from 1, so PEs count from 1 to nPE
             PE = PGSLib_PE_From_Global_Index(INDEX(R1,R2), GID) 
             PE_ARRAY(R1, R2) = PE 
          end if
          
          ! If we do not own this index, put it into table
          IF (PE_ARRAY(R1, R2) /= thisPE) then
             ierror =  PGSLib_Add_Item_To_Table(INDEX(R1,R2), PE_ARRAY(R1, R2), Table)
          END IF

       ENDDO
    ENDDO

    N_Supplement = PGSLib_Count_Items_In_Table(Table)
#ifdef USE_TIMERS_1
    call Exit_ROUTINE(Supplement_G_TO_L_Statistics())
#endif

#ifdef DEBUG_GATH_SCATT
    write(pgslib_out_string, *) 'N_Supplement = ', N_Supplement
    call pgslib_output(pgslib_out_string)
#endif
    ! Since we know how many items we need to get from other PEs we can allocate
    !   space for the communication buffers.

    ALLOCATE(Sup_Global_I (N_Supplement))
    ALLOCATE(Sup_Global_PE(N_Supplement))

    ! We need to get the global INDEX from the offPE_Access_Table and 
    ! put them into Supplement_Global_Indices.  At the same time
    ! we need to put the new index into the offPE_Access_Table so that
    ! it can be accessed by PGSLib_Item_Index_From_OffPE_Access_Table.
    ! This is done in a single call for performance reasons.
    IF (N_Supplement .GT. 0) THEN
       ierror= PGSLib_Items_From_Table(Sup_Global_I, Sup_Global_PE,  &
            &                       Table)
    ENDIF

#ifdef USE_TIMERS_1
    call Enter_ROUTINE(Supplement_L_From_G_Statistics())
#endif

    ! The final part to this step is to renumber INDEX to point to either
    ! the local Compliors or to point to the Supplement buffer.
    ! If the point to the Supplement they are made negative to flag that.
    DO R2 = 1, SIZE(INDEX,2)
       Do R1 = 1, SIZE(INDEX,1)
          ! If Mask, then if this requestor is inactive, skip to next one
          if (Mask_Present) then
             if (.NOT. MASK(R1, R2) ) CYCLE
          endif

          ! Only do this if no mask, or if Mask(R1, R2) == TRUE
          IF (PE_ARRAY(R1, R2) == thisPE) THEN
             INDEX(R1, R2) = PGSLib_Index_Local_From_Global(INDEX(R1,R2), GID, PE=thisPE)
          ELSE
             INDEX(R1, R2) = -PGSLib_Item_Index_From_Table(INDEX(R1,R2), PE_ARRAY(R1, R2), Sup_Global_I, Sup_Global_PE)
          ENDIF
       ENDDO
    ENDDO

    ! Done with the access table, since we figured got what we need from it.
    call pgslib_free_offPE_access_table(Table)


#ifdef USE_TIMERS_1
    call Exit_ROUTINE(Supplement_L_From_G_Statistics())
#endif

#ifdef USE_TIMERS_1
    call Exit_Routine(Setup_Supplement_STATISTICS())
#endif    

    RETURN
  END SUBROUTINE PGSLib_Setup_Supplement_F

!======================================================================
!          PGSLib_Setup_Duplicate(Trace)
! PURPOSE:
!          Establish the Duplicate buffer in the trace.  This is done
!          by sending information obtained from the Supplement buffer
!          to the appropriate destinations.  On each PE, this information is used to
!          determine what items have to be served up to other PE''s during a scatter,
!          and which items will be arriving during a gather.  That is contained
!          in the Duplicate communication buffer.
!======================================================================
  subroutine PGSLib_Setup_Duplicate_F(Trace, GID)
    USE PGSLib_Type_MODULE
    USE PGSLib_Globals_MODULE
    USE PGSLib_Utility_MODULE
    USE PGSLib_Index_MODULE
    USE PGSLib_GS_Comm_MODULE
    USE PGSLIB_GS_Util_MODULE
    implicit none 
    type (PGSLib_GS_Trace), INTENT(INOUT) :: Trace
    type (PGSLib_GID),      INTENT(IN)    :: GID
    ! Local variables
    integer :: thisPE, Dup

    ! Allocate space for data buffers.
    ALLOCATE(Trace%Duplicate_Indices(PGSLib_Size_Of_Dup(Trace)))
    ALLOCATE(Trace%Duplicate_I_DATA(PGSLib_Size_Of_Dup(Trace)))
    ALLOCATE(Trace%Duplicate_R_DATA(PGSLib_Size_Of_Dup(Trace)))
    ALLOCATE(Trace%Duplicate_D_DATA(PGSLib_Size_Of_Dup(Trace)))
    ALLOCATE(Trace%Duplicate_L_DATA(PGSLib_Size_Of_Dup(Trace)))

    ! Get the indices for the duplicates.  These will arrive as global indices,
    ! but we will fix that below.
    Trace%Duplicate_Indices = pgslib_scatter_buffer(PGSLib_Sup_Global_Indices(Trace), Trace)

    ! Finally, we have to turn the global indices we got into local indices.
    ! We need a GID for this.  We should have one by the time we get here.
    thisPE = PGSLib_Inquire_thisPE_Actual() 
    DO Dup=1, PGSLib_Size_Of_Dup(Trace)
       Trace%Duplicate_Indices(Dup) = PGSLib_Index_Local_From_Global(Trace%Duplicate_Indices(Dup), &
            &                                                        GID, PE=thisPE)
    ENDDO

    RETURN
  END subroutine PGSLib_Setup_Duplicate_F

  !=======================================================================
  !          PGSLib_Setup_Buffers
  !=======================================================================
  function PGSLib_Setup_Basic_Trace(PE_ARRAY) RESULT(Trace)
    
    ! Establish a trace for buffers
    
    USE PGSLib_Type_MODULE
    USE PGSLib_GS_Comm_MODULE
    use pgslib_gs_util_module
    implicit none 

    integer (PGSLib_Int_Type), intent(IN   ), &
                               dimension(:)   :: PE_ARRAY
    type (PGSLib_GS_Trace),    POINTER        :: Trace

    ! Local variables
    integer (PGSLib_Int_Type) :: N_Slaves

    ! First need to allocate a trace.  Use RESULT Trace, for short name.
    Trace => PGSLib_GS_Init_Trace()

#ifdef USE_TIMERS_1
    call Enter_Routine(SETUP_BUFFERS_STATISTICS())
    call PGSLib_Start_Timer(Trace%Setup_Timer)
#endif

    Trace%N_Supplement = SIZE(PE_ARRAY,1)
    ! C counts PE_ARRAY from 0 to nPE-1, FORTRAN from 1 to nPE
    call pgslib_Prep_supplement(SIZE(PE_ARRAY,1), PE_ARRAY-1, Trace%GS_Trace)

    call pgslib_setup_n_duplicate(N_Slaves, Trace%GS_Trace)
    Trace%N_Duplicate = N_Slaves

    call pgslib_setup_duplicate_buffer(Trace%GS_Trace)

    ! Indicate that trace has been setup.
    Trace%SetupP = .TRUE.

#ifdef USE_TIMERS_1
    call Exit_Routine(SETUP_BUFFERS_STATISTICS())
    call PGSLib_Stop_Timer(Trace%Setup_Timer)
#endif

    return
  end function PGSLib_Setup_Basic_Trace

END MODULE PGSLib_GS_Setup_MODULE
    

