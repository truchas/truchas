MODULE  PGSLib_Index_GID_MODULE 
  use pgslib_stats
  use pgslib_timing_module
 
  implicit none
  save
  private
  public :: PGSLib_Setup_GID
  public :: PGSLib_Deallocate_GID
  public :: PGSLib_Release_GID
  public :: PGSLib_Init_GID
  public :: PGSLib_Index_Local_From_Global
  public :: PGSLib_Index_Global_From_Local
  public :: PGSLib_PE_From_Global_Index

!======================================================================
!          PURPOSE
!            These routines are for setup and use of the Global
!            Index Descriptor (GID).
!======================================================================

  !$Id: pgslib_index_gid_module.F,v 1.1.1.1 2000/10/11 22:44:31 ferrell Exp $

  INTERFACE PGSLib_Setup_GID
     MODULE PROCEDURE PGSLib_Setup_GID_F
  END INTERFACE ! PGSLib_Setup_GID

  INTERFACE PGSLib_Deallocate_GID
     MODULE PROCEDURE PGSLib_Deallocate_GID_F
  END INTERFACE ! PGSLib_Deallocate_GID

  INTERFACE PGSLib_Release_GID
     MODULE PROCEDURE PGSLib_Release_GID_F
  END INTERFACE ! PGSLib_Release_GID

  INTERFACE PGSLib_Init_GID
     MODULE PROCEDURE PGSLib_Init_GID_F
  END INTERFACE ! PGSLib_Init_GID

  INTERFACE PGSLib_Index_Local_From_Global
     MODULE PROCEDURE PGSLib_Ind_Loc_From_Glob
  END INTERFACE ! PGSLib_Index_Local_From_Global

  INTERFACE PGSLib_PE_From_Global_Index
     MODULE PROCEDURE PGSLib_PE_From_Global_Index_S
     MODULE PROCEDURE PGSLib_PE_From_Global_Index_V1
     MODULE PROCEDURE PGSLib_PE_From_Global_Index_V2
  END INTERFACE ! PGSLib_PE_From_Global_Index

  INTERFACE PGSLib_Index_Global_From_Local
     MODULE PROCEDURE PGS_Ind_Glb_From_Loc_F
  END INTERFACE

  CONTAINS


!======================================================================
!======================================================================
!          PGSLib_Setup_GID(GID, LocalExtents)
!        PURPOSE
!          Establish a Global Index Descriptor (GID) for a data structure
!          Space is allocated, and the data fields are filled in.
!          The GID describes the distribution of a 1D data structure accross
!          the processing elements.  The distribution is fully specified
!          by the LocalExtents on each PE.  The GID is replicated to all
!          PEs.  
!======================================================================


    function PGSLib_Setup_GID_F(LocalExtent) RESULT(GID)
      USE PGSLib_TYPE_MODULE
      USE PGSLib_Globals_MODULE,   ONLY : PGSLib_PEInfo
      USE PGSLib_Utility_MODULE,   ONLY : PGSLib_Fatal_Error
      USE PGSLib_IO_MODULE,        ONLY : PGSLib_Bcast,      &
           &                              PGSLib_Collate
      implicit none
      INTEGER (PGSLib_Int_Type), INTENT(IN) :: LocalExtent
      Type (PGSLib_GID),         POINTER    :: GID

      ! Local variables
      Integer (PGSLib_Int_Type) :: N
      
      ! Allocate the GID, and initialize it.
      ALLOCATE(GID)
      call PGSLib_Init_GID(GID)
      
      ! Allocate space for the GID pieces 
      ALLOCATE(GID%LocalExtents(PGSLib_PEInfo%nPE))
      ALLOCATE(GID%LowIndex   (PGSLib_PEInfo%nPE))
      ALLOCATE(GID%HighIndex  (PGSLib_PEInfo%nPE))
      
      ! Determine global extents by collating data from all PE''s.
      ! Get the local node extents from each PE, add them up to get the global size.
      GID%ExtentThisPE = LocalExtent
      call PGSLib_Collate(GID%LocalExtents, GID%ExtentThisPE)
      call PGSLib_BCast  (GID%LocalExtents)
      GID%TotalExtent = SUM(GID%LocalExtents)
      
      ! The LowIndex and HighIndex are redundant in this version (since they can be
      !          constructed from LocalExtents), but it is useful to have them around.
      ! This code assumes that the data items are in global order (sequential on each PE)
      
      GID%LowIndex(1) = 1
      
      DO N = 1, PGSLib_PEInfo%nPE -1
         GID%HighIndex(N) =  GID%LowIndex(N) +  (GID%LocalExtents(N) - 1)
         GID%LowIndex(N+1) = GID%HighIndex(N) + 1
      ENDDO
      
      GID%HighIndex(PGSLib_PEInfo%nPE) = GID%LowIndex(PGSLib_PEInfo%nPE) + &
           &                            (GID%LocalExtents(PGSLib_PEInfo%nPE) - 1)
      
      
      ! HighIndex should be the last index  If not, we made an error
      IF (GID%HighIndex(PGSLib_PEInfo%nPE) /= GID%TotalExtent) THEN
         call pgslib_fatal_error("In pgslib_setup_gid, HighIndex /= TotalExtent ")
      ENDIF
      
      ! If we got this far, then the GID is setup.
      GID%SetupP = .TRUE.

      RETURN
    END function PGSLib_Setup_GID_F
    
    !=====================================================================
    !=====================================================================
    !         PGSLib_DEALLOCATE_GID(GID)
    ! PURPOSE
    !   This routine releases the storage for GID.  After this call GID
    !   may not be used again until after another call to PGSLIB_SETUP_GID.
    !======================================================================
    subroutine PGSLib_DEALLOCATE_GID_F(GID)
      USE PGSLib_Type_MODULE
      implicit none
      type (PGSLib_GID), POINTER :: GID

      call PGSLib_RELEASE_GID(GID)
      if (ASSOCIATED(GID)) Deallocate(GID)

      RETURN
    END subroutine PGSLib_DEALLOCATE_GID_F

    !=====================================================================
    !         PGSLib_Release_GID(GID)
    ! PURPOSE
    !   This routine releases the storage for GID.  It does not deallocate
    !   the GID itself.
    !======================================================================
    subroutine PGSLib_Release_GID_F(GID)
      USE PGSLib_Type_MODULE
      implicit none
      type (PGSLib_GID) :: GID

      if (ASSOCIATED(GID%LocalExtents)) deallocate(GID%LocalExtents)
      if (ASSOCIATED(GID%LowIndex    )) deallocate(GID%LowIndex    )
      if (ASSOCIATED(GID%HighIndex   )) deallocate(GID%HighIndex   )

      call PGSLib_Init_GID(GID)

      RETURN
    END subroutine PGSLib_Release_GID_F

    !=====================================================================
    !=====================================================================
    !         PGSLib_Init_GID(GID)
    ! PURPOSE
    !   This routine Inits the storage for GID.  After this call GID
    !   may be released or setup safely.n
    !======================================================================
    subroutine PGSLib_Init_GID_F(GID)
      USE PGSLib_Type_MODULE
      implicit none
      type (PGSLib_GID) :: GID

      GID%SetupP = .FALSE.
      NULLIFY(GID%LocalExtents)
      NULLIFY(GID%LowIndex    )
      NULLIFY(GID%HighIndex   )
      GID%TotalExtent  = -1
      GID%ExtentThisPE = -1

      RETURN
    END subroutine PGSLib_Init_GID_F

    !=====================================================================
    !=====================================================================
    !          integer PGSLib_Index_Local_From_Global(GlobalIndex, PE, GID)
    !
    ! PURPOSE
    ! This routine takes a global canonical index and returns
    ! the local index.  The PE is OPTIONAL.  However, if it is
    ! known, then the routine is often faster.
    ! This version assumes a simple GID with each PE containing a compact
    ! range of global indices.
    !=====================================================================

    FUNCTION PGSLib_Ind_Loc_From_Glob(GI, GID, PE) RESULT(LocalIndex)
      USE PGSLib_Type_MODULE
      USE PGSLib_Globals_MODULE,   ONLY : PGSLib_PEInfo
      USE PGSLib_Utility_MODULE,   ONLY : PGSLib_Fatal_Error,   &
           &                              PGSLib_Inquire_nPE
      IMPLICIT NONE 
      integer (PGSLib_Int_Type)                       :: LocalIndex
      integer (PGSLib_Int_Type), INTENT(IN)           :: GI
      type (PGSLib_GID),         INTENT(IN)           :: GID
      integer (PGSLib_Int_Type), INTENT(IN), OPTIONAL :: PE

      ! Local variables
      integer (PGSLib_Int_Type) :: lpe


      ! FORTRAN counts processors from 1 nPE
      IF (PRESENT(PE)) THEN
         IF (  (GI .GE. GID%LowIndex(PE)) .AND. &
              &(GI .LE. GID%HighIndex(PE)) ) THEN
            LocalIndex = GI - GID%LowIndex(PE) + 1
         ELSE  ! Error, since the GI is not on the PE the user specified.
            Call PGSLib_Fatal_Error("Wrong PE specified in PGSLib_Index_Local_From_Global")
         ENDIF
         
      ELSE
         ! Need to make this a binary search rather than linear search
         DO lPE=1, PGSLib_Inquire_nPE()   ! If loop gets this far there has been an error
            IF (  (GI .GE. GID%LowIndex(lPE)) .AND. &
                 &(GI .LE. GID%HighIndex(lPE)) ) EXIT
         ENDDO
         
         IF (lPE .GT. PGSLib_Inquire_nPE()) THEN  ! Index was not found
            call pgslib_Fatal_Error("Failed to find local index in PGSLib_Index_Local_From_Global")
         ELSE
            LocalIndex = GI - GID%LowIndex(lPE) + 1
         ENDIF
      ENDIF
      
      RETURN
    END FUNCTION PGSLib_Ind_Loc_From_Glob
               
            
    !=====================================================================
    !=====================================================================
    !          integer PGSLib_PE_From_Global_Index(GI, GID)
    !
    ! PURPOSE
    ! Returns the PE which holds the item pointed to by global index GI.
    ! The global index must be in Global Canonical Index format.
    ! This version assumes a simple GID with each PE containing a compact
    ! range of global indices.
    !=====================================================================

    function PGSLib_PE_From_Global_Index_S(INDEX, GID)
      USE PGSLib_Type_MODULE
      USE PGSLib_Globals_MODULE,   ONLY : PGSLib_PEInfo
      USE PGSLib_Utility_MODULE,   ONLY : PGSLib_Fatal_Error,   &
           &                              PGSLib_Inquire_nPE,   &
           &                              pgslib_out_string
      IMPLICIT NONE 
      integer (PGSLib_Int_Type)                       :: PGSLib_PE_From_Global_Index_S
      integer (PGSLib_Int_Type), INTENT(IN)           :: INDEX
      type (PGSLib_GID) , INTENT(IN)                  :: GID

      ! Local variables
      integer (PGSLib_Int_Type) :: lpe
      ! Need to make this a binary search rather than linear search
      DO lPE=1, PGSLib_Inquire_nPE()   ! If loop gets this far there has been an error
         IF (  (INDEX .GE. GID%LowIndex(lPE)) .AND. &
              &(INDEX .LE. GID%HighIndex(lPE)) ) EXIT
      ENDDO
         
      IF (lpe .GT. PGSLib_Inquire_nPE()) THEN  ! Index was not found
         write(pgslib_out_string,'("Failed to find global index",i8," in PGSLib_PE_From_Global_Index")') INDEX
         call pgslib_Fatal_Error(pgslib_out_string)
      ELSE
         PGSLib_PE_From_Global_Index_S = lPE
      ENDIF

    END function PGSLib_PE_From_Global_Index_S
      

    !====================================================================
    !            Vector 1D version
    !====================================================================

    function PGSLib_PE_From_Global_Index_V1(INDICES, SIZE_OF_DEST, GID)
      USE PGSLib_Type_MODULE
      USE PGSLib_Globals_MODULE,   ONLY : PGSLib_PEInfo
      USE PGSLib_Utility_MODULE,   ONLY : PGSLib_Fatal_Error,   &
           &                              PGSLib_Inquire_nPE,   &
           &                              pgslib_out_string
      IMPLICIT NONE 
      integer (PGSLib_Int_Type), INTENT(IN),          &
                                 DIMENSION(:)         :: INDICES
      integer (PGSLib_Int_Type), INTENT(IN),          &
                                 OPTIONAL             :: SIZE_OF_DEST
      type (PGSLib_GID) , INTENT(IN),                  &
                                 TARGET,               &
                                 OPTIONAL             :: GID

      integer (PGSLib_Int_Type), DIMENSION(SIZE(INDICES,1)) :: PGSLib_PE_From_Global_Index_V1

      ! Local variables
      type (PGSLib_GID), POINTER :: Local_GID
      integer (PGSLib_Int_Type)  :: i

      ! We need either a GID or SIZE_OF_DEST, but not both
      if (PRESENT(SIZE_OF_DEST) .AND. PRESENT(GID)) then
         call pgslib_Fatal_Error('Provide only one of SIZE_OF_DEST or GID in PGSLib_PE_From_Global_Index')
      end if
      
      if (.NOT.PRESENT(SIZE_OF_DEST) .AND. .NOT.PRESENT(GID)) then
         call pgslib_Fatal_Error('Provide either SIZE_OF_DEST or GID in PGSLib_PE_From_Global_Index')
      end if
      
      ! If we did not get a GID, we need to make one
      if (.NOT. PRESENT(GID)) then
         Local_GID => PGSLib_Setup_GID(SIZE_OF_DEST)
      else
         Local_GID => GID
      end if
      
      ! Now that we have a GID, we can get the PEs one at a time
      do i = 1, SIZE(INDICES,1)
         PGSLib_PE_From_Global_Index_V1(i) = PGSLib_PE_From_Global_Index(INDICES(I), Local_GID)
      end do
      
      ! If we made a local GID, then get rid of it.
      if (.NOT. PRESENT(GID)) call PGSLib_Deallocate_GID(Local_GID)

      RETURN
    END function PGSLib_PE_From_Global_Index_V1
      

    !====================================================================
    !            Vector 2D version
    !====================================================================

    function PGSLib_PE_From_Global_Index_V2(INDICES, SIZE_OF_DEST, GID)
      USE PGSLib_Type_MODULE
      USE PGSLib_Globals_MODULE,   ONLY : PGSLib_PEInfo
      USE PGSLib_Utility_MODULE,   ONLY : PGSLib_Fatal_Error,   &
           &                              PGSLib_Inquire_nPE,   &
           &                              pgslib_out_string
      IMPLICIT NONE 
      integer (PGSLib_Int_Type), INTENT(IN),          &
                                 DIMENSION(:,:)         :: INDICES
      integer (PGSLib_Int_Type), INTENT(IN),          &
                                 OPTIONAL             :: SIZE_OF_DEST
      type (PGSLib_GID) , INTENT(IN),                  &
                                 TARGET,               &
                                 OPTIONAL             :: GID

      integer (PGSLib_Int_Type), DIMENSION(SIZE(INDICES,1),&
                                           SIZE(INDICES,2)) :: PGSLib_PE_From_Global_Index_V2

      ! Local variables
      type (PGSLib_GID), POINTER :: Local_GID
      integer (PGSLib_Int_Type)  :: i

      ! We need either a GID or SIZE_OF_DEST, but not both
      if (PRESENT(SIZE_OF_DEST) .AND. PRESENT(GID)) then
         call pgslib_Fatal_Error('Provide only one of SIZE_OF_DEST or GID in PGSLib_PE_From_Global_Index')
      end if
      
      if (.NOT.PRESENT(SIZE_OF_DEST) .AND. .NOT.PRESENT(GID)) then
         call pgslib_Fatal_Error('Provide either SIZE_OF_DEST or GID in PGSLib_PE_From_Global_Index')
      end if
      
      ! If we did not get a GID, we need to make one
      if (.NOT. PRESENT(GID)) then
         Local_GID => PGSLib_Setup_GID(SIZE_OF_DEST)
      else
         Local_GID => GID
      end if
      
      ! Now that we have a GID, we can get the PEs one row at a time
      do i = 1, SIZE(INDICES,1)
         PGSLib_PE_From_Global_Index_V2(i,:) = PGSLib_PE_From_Global_Index(Indices(i,:), GID=Local_GID)
      end do
      
      ! If we made a local GID, then get rid of it.
      if (.NOT. PRESENT(GID)) call PGSLib_Deallocate_GID(Local_GID)

      RETURN
    END function PGSLib_PE_From_Global_Index_V2
      

    !=====================================================================
    !=====================================================================
    !          integer PGSLib_Index_Global_From_Local(LocalIndex, PE, GID)
    !
    ! PURPOSE
    ! This routine takes a  local index and returns
    ! the global canonical index.  The PE is OPTIONAL.  If it is not
    ! provided the local PE is assumed.
    ! This version assumes a simple GID with each PE containing a compact
    ! range of global indices.
    !=====================================================================
    function PGS_Ind_Glb_From_Loc_F(LocalIndex, GID, PE)
      USE PGSLib_Type_MODULE
      USE PGSLib_Globals_MODULE,   ONLY : PGSLib_PEInfo
      USE PGSLib_Utility_MODULE,   ONLY : PGSLib_Fatal_Error
      implicit none
      integer (PGSLib_Int_Type)                       :: PGS_Ind_Glb_From_Loc_F
      integer (PGSLib_Int_Type), intent(IN)           :: LocalIndex
      type  (PGSLib_GID), intent(IN)           :: GID
      integer (PGSLib_Int_Type), intent(IN), OPTIONAL :: PE
      
      ! Local variables
      
      IF (PRESENT(PE)) THEN
         PGS_Ind_Glb_From_Loc_F = GID%LowIndex(PE+1) + LocalIndex - 1
      ELSE
         PGS_Ind_Glb_From_Loc_F = GID%LowIndex(PGSLib_PEInfo%thisPE) + LocalIndex - 1
      END IF

      RETURN
    END function PGS_Ind_Glb_From_Loc_F

  END MODULE PGSLib_Index_GID_MODULE

         







