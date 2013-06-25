  MODULE UPDATE_BCS
    !-----------------------------------------------------------------------
    ! 
    ! Purpose:
    !   Update boundary condition atlases to reflect current conditions
    !
    ! PUBLIC:  Update_Radiation_BC(time)
    !          Update_Dirichlet_BC(time)
    !
    !    Jim Sicilian, CCS-2 (sicilian@lanl.gov)
    !    May 2003 and November 2004
    !
    !----------------------------------------------------------------------

    Private
    Public Update_Radiation_BC, Update_Dirichlet_BC, Update_HTC_External_BC

    CONTAINS

    SUBROUTINE Update_Radiation_BC(time)
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    !   Update Radiation Boundary Environment Temperatures
    !
    !    Jim Sicilian, CCS-2 (sicilian@lanl.gov)
    !
    !-----------------------------------------------------------------------
    use kind_module,            only: real_kind, int_kind
    use bc_operations,          only: BC_Operator, BC_Atlas, BC_Spec_Get_Operator, &
                                      BC_Op_Get_Atlas, BC_Get_Face, BC_Get_Cell,   &
                                      BC_Get_Offset, BC_Get_Values, &
                                      DATA_SIZE, BC_RADIATION_Op,&
                                      TEMPERATURE_BC, BC_Get_DOF
    use tabular_utilities,      only: tabular_linear_interp
    use truchas_logging_services

    implicit none

    ! arguments
    real(kind=real_kind), INTENT(IN)           :: time
    ! local variables
    real(real_kind)                            :: tlast
    integer(int_kind)                          :: k,n, status
    integer(int_kind)                          :: DegreesOfFreedom, InterpolationPoints

    ! BC-related local variables.
    integer(int_kind)                          :: NumberBdyPts, BdyPt, BdyCell
    integer(int_kind)                          :: BdyFace, BdyOffSet
    type(BC_Operator), pointer                 :: RADIATION_Operator
    type(BC_Atlas),    pointer                 :: RADIATION_Atlas
    integer(int_kind), pointer, dimension(:)   :: BdyFaceList
    integer(int_kind), pointer, dimension(:)   :: BdyCellList
    integer(int_kind), pointer, dimension(:)   :: BdyOffsetList
    real(real_kind),   pointer, dimension(:,:) :: BdyValueList

    real(kind=real_kind), allocatable, SAVE, dimension(:)  :: timelist, templist

    ! Apply the RADIATION Operator
    ! Get the RADIATION operator.
    RADIATION_Operator => BC_Spec_Get_Operator(TEMPERATURE_BC, BC_RADIATION_Op)
    
    ! Since we are doing overwrite, we can do all boundary faces in one loop.
    ! So we want to traverse the whole atlas.
    RADIATION_Atlas => BC_OP_Get_Atlas(RADIATION_Operator)
    
    ! We are using whole arrays, so faster to point at them than get a private copy.
    BdyFaceList     => BC_Get_Face(RADIATION_Atlas)
    BdyCellList     => BC_Get_Cell(RADIATION_Atlas)
    BdyOffsetList   => BC_Get_Offset(RADIATION_Atlas)
    BdyValueList    => BC_Get_Values(RADIATION_Atlas)
    DegreesOfFreedom = BC_Get_DOF(RADIATION_Atlas)
    
    InterpolationPoints = DegreesOfFreedom - 2

    if(.not.allocated(timelist)) then
       allocate(timelist(InterpolationPoints), STAT = status)
       if (status /= 0) call TLS_panic ('Update_Radiation_BC: timelist(InterpolationPoints) allocate failed')
       allocate(templist(InterpolationPoints), STAT = status)
       if (status /= 0) call TLS_panic ('Update_Radiation_BC: templist(InterpolationPoints) allocate failed')
    endif

    NumberBdyPts = DATA_SIZE(RADIATION_Atlas)
    RAD_BDY_LOOP: do BdyPt = 1, NumberBdyPts
       BdyCell = BdyCellList(BdyPt)
       BdyFace = BdyFaceList(BdyPt)
       BdyOffSet = BdyOffsetList(BdyPt)
       ! create lists of time and temperature values
       k=3
       n=1
       tlast = -1.0
       COUNT_VALUES:  do while (k<InterpolationPoints)
           if(BdyValueList(k,BdyOffset) <= tlast) then
               exit COUNT_VALUES
           endif
           timelist(n) = BdyValueList(k, BdyOffset)
           tlast = timelist(n)
           templist(n) = BdyValueList(k+1, BdyOffset)
           n=n+1
           k=k+2
       end do COUNT_VALUES
       if(n-1>1) then
           BdyValueList(2, BdyOffset) = tabular_linear_interp(time, timelist(1:n-1), templist(1:n-1))
       else
!  needed to update model for sensitivity analysis
           BdyValueList(2, BdyOffset) = BdyValueList(4, BdyOffset) 
       endif
    end do RAD_BDY_LOOP

!    call Output_Atlas_To_XML(RADIATION_Atlas)
    return
    end subroutine Update_Radiation_Bc

    SUBROUTINE Update_Dirichlet_BC(time)
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    !   Update Dirichlet Boundary Temperatures
    !
    !    Jim Sicilian, CCS-2 (sicilian@lanl.gov)
    !
    !-----------------------------------------------------------------------
    use kind_module,            only: real_kind, int_kind
    use bc_operations,          only: BC_Operator, BC_Atlas, BC_Spec_Get_Operator, &
                                      BC_Op_Get_Atlas, BC_Get_Face, BC_Get_Cell,   &
                                      BC_Get_Offset, BC_Get_Values, &
                                      DATA_SIZE, BC_DIRICHLET_Op,&
                                      TEMPERATURE_BC, BC_Get_DOF
    use tabular_utilities,      only: tabular_linear_interp
    use truchas_logging_services

    implicit none

    ! arguments
    real(kind=real_kind), INTENT(IN)           :: time
    ! local variables
    real(real_kind)                            :: tlast
    integer(int_kind)                          :: k,n, status
    integer(int_kind)                          :: DegreesOfFreedom, InterpolationPoints

    ! BC-related local variables.
    integer(int_kind)                          :: NumberBdyPts, BdyPt, BdyCell
    integer(int_kind)                          :: BdyFace, BdyOffSet
    type(BC_Operator), pointer                 :: DIRICHLET_Operator
    type(BC_Atlas),    pointer                 :: DIRICHLET_Atlas
    integer(int_kind), pointer, dimension(:)   :: BdyFaceList
    integer(int_kind), pointer, dimension(:)   :: BdyCellList
    integer(int_kind), pointer, dimension(:)   :: BdyOffsetList
    real(real_kind),   pointer, dimension(:,:) :: BdyValueList

    real(kind=real_kind), allocatable, SAVE, dimension(:)  :: timelist, templist

    ! Apply the DIRICHLET Operator
    ! Get the Dirichlet operator.
    DIRICHLET_Operator => BC_Spec_Get_Operator(TEMPERATURE_BC, BC_DIRICHLET_Op)
    
    ! Since we are doing overwrite, we can do all boundary faces in one loop.
    ! So we want to traverse the whole atlas.
    DIRICHLET_Atlas => BC_OP_Get_Atlas(DIRICHLET_Operator)
    
    ! We are using whole arrays, so faster to point at them than get a private copy.
    BdyFaceList     => BC_Get_Face(DIRICHLET_Atlas)
    BdyCellList     => BC_Get_Cell(DIRICHLET_Atlas)
    BdyOffsetList   => BC_Get_Offset(DIRICHLET_Atlas)
    BdyValueList    => BC_Get_Values(DIRICHLET_Atlas)
    DegreesOfFreedom = BC_Get_DOF(DIRICHLET_Atlas)
    
    InterpolationPoints = DegreesOfFreedom - 1

    if(.not.allocated(timelist)) then
       allocate(timelist(InterpolationPoints), STAT = status)
       if (status /= 0) call TLS_panic ('Update_Dirichlet_BC: timelist(InterpolationPoints) allocate failed')
       allocate(templist(InterpolationPoints), STAT = status)
       if (status /= 0) call TLS_panic ('Update_Dirichlet_BC: templist(InterpolationPoints) allocate failed')
    endif

    NumberBdyPts = DATA_SIZE(DIRICHLET_Atlas)
    RAD_BDY_LOOP: do BdyPt = 1, NumberBdyPts
       BdyCell = BdyCellList(BdyPt)
       BdyFace = BdyFaceList(BdyPt)
       BdyOffSet = BdyOffsetList(BdyPt)
       ! create lists of time and temperature values
       k=2
       n=1
       tlast = -1.0
       COUNT_VALUES:  do while (k<InterpolationPoints)
           if(BdyValueList(k,BdyOffset) <= tlast) then
               exit COUNT_VALUES
           endif
           timelist(n) = BdyValueList(k, BdyOffset)
           tlast = timelist(n)
           templist(n) = BdyValueList(k+1, BdyOffset)
           n=n+1
           k=k+2
       end do COUNT_VALUES
       if(n-1>1) then
           BdyValueList(1, BdyOffset) = tabular_linear_interp(time, timelist(1:n-1), templist(1:n-1))
       else
!  needed to update model for sensitivity analysis
           BdyValueList(1, BdyOffset) = BdyValueList(3, BdyOffset)
       endif
    end do RAD_BDY_LOOP

!    call Output_Atlas_To_XML(DIRICHLET_Atlas)
    return
    end subroutine Update_Dirichlet_Bc

    SUBROUTINE Update_HTC_External_BC(time)
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    !   Update Dirichlet Boundary Temperatures
    !
    !    Jim Sicilian, CCS-2 (sicilian@lanl.gov)
    !
    !-----------------------------------------------------------------------
    use kind_module,            only: real_kind, int_kind
    use bc_operations,          only: BC_Operator, BC_Atlas, BC_Spec_Get_Operator, &
                                      BC_Op_Get_Atlas, BC_Get_Face, BC_Get_Cell,   &
                                      BC_Get_Offset, BC_Get_Values, &
                                      DATA_SIZE, BC_HTC_EXTERNAL_Op,&
                                      TEMPERATURE_BC, BC_Get_DOF
    use tabular_utilities,      only: tabular_linear_interp
    use truchas_logging_services

    implicit none

    ! arguments
    real(kind=real_kind), INTENT(IN)           :: time
    ! local variables
    real(real_kind)                            :: tlast
    integer(int_kind)                          :: k,n, status
    integer(int_kind)                          :: DegreesOfFreedom, InterpolationPoints

    ! BC-related local variables.
    integer(int_kind)                          :: NumberBdyPts, BdyPt, BdyCell
    integer(int_kind)                          :: BdyFace, BdyOffSet
    type(BC_Operator), pointer                 :: HTC_External_Operator
    type(BC_Atlas),    pointer                 :: HTC_External_Atlas
    integer(int_kind), pointer, dimension(:)   :: BdyFaceList
    integer(int_kind), pointer, dimension(:)   :: BdyCellList
    integer(int_kind), pointer, dimension(:)   :: BdyOffsetList
    real(real_kind),   pointer, dimension(:,:) :: BdyValueList

    real(kind=real_kind), allocatable, SAVE, dimension(:)  :: timelist, templist

    ! Apply the HTC_External Operator
    ! Get the HTC_External operator.
    HTC_External_Operator => BC_Spec_Get_Operator(TEMPERATURE_BC, BC_HTC_EXTERNAL_Op)
    
    ! Since we are doing overwrite, we can do all boundary faces in one loop.
    ! So we want to traverse the whole atlas.
    HTC_External_Atlas => BC_OP_Get_Atlas(HTC_External_Operator)
    
    ! We are using whole arrays, so faster to point at them than get a private copy.
    BdyFaceList     => BC_Get_Face(HTC_External_Atlas)
    BdyCellList     => BC_Get_Cell(HTC_External_Atlas)
    BdyOffsetList   => BC_Get_Offset(HTC_External_Atlas)
    BdyValueList    => BC_Get_Values(HTC_External_Atlas)
    DegreesOfFreedom = BC_Get_DOF(HTC_External_Atlas)
    
    InterpolationPoints = DegreesOfFreedom - 2

    if(.not.allocated(timelist)) then
       allocate(timelist(InterpolationPoints), STAT = status)
       if (status /= 0) call TLS_panic ('Update_HTC_External_BC: timelist(InterpolationPoints) allocate failed')
       allocate(templist(InterpolationPoints), STAT = status)
       if (status /= 0) call TLS_panic ('Update_HTC_External_BC: templist(InterpolationPoints) allocate failed')
    endif

    NumberBdyPts = DATA_SIZE(HTC_External_Atlas)
    RAD_BDY_LOOP: do BdyPt = 1, NumberBdyPts
       BdyCell = BdyCellList(BdyPt)
       BdyFace = BdyFaceList(BdyPt)
       BdyOffSet = BdyOffsetList(BdyPt)
       ! create lists of time and temperature values
       k=3
       n=1
       tlast = -1.0
       COUNT_VALUES:  do while (k<InterpolationPoints)
           if(BdyValueList(k,BdyOffset) <= tlast) then
               exit COUNT_VALUES
           endif
           timelist(n) = BdyValueList(k, BdyOffset)
           tlast = timelist(n)
           templist(n) = BdyValueList(k+1, BdyOffset)
           n=n+1
           k=k+2
       end do COUNT_VALUES
       if(n-1>1) then
           BdyValueList(2, BdyOffset) = tabular_linear_interp(time, timelist(1:n-1), templist(1:n-1))
       else
!  needed to update model for sensitivity analysis
           BdyValueList(2, BdyOffset) = BdyValueList(4, BdyOffset)
       endif
    end do RAD_BDY_LOOP

!    call Output_Atlas_To_XML(HTC_External_Atlas)
    return
  end subroutine Update_HTC_External_BC

  END MODULE UPDATE_BCS
