MODULE GS_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define the interfaces and global variables 
  !   for the gather/scatter support
  !
  !=======================================================================
  use gather_module,  only: Gather_BoundaryData
  use gs_util,        only: EE_GS_INIT,      &
                            EN_GS_INIT,      &
                            NN_GS_INIT
  use ee_gather_module, only: EE_Gather
  use en_gather_module, only: EN_Gather,      &
                              EN_MIN_Gather,  &
                              EN_MAX_Gather,  &
                              EN_OR_Gather
  use en_scatter_module, only: EN_SUM_Scatter, &
                               EN_MIN_Scatter, &
                               EN_MAX_Scatter, &
                               EN_OR_Scatter

  use nn_gather_module, only: NN_Gather, NN_Gather_BoundaryData

  implicit none
  private

  ! Public procedures
  Public :: EE_Gather,      &
            EN_Gather,      &
            NN_Gather,      &
            EN_MIN_Gather,  &
            EN_MAX_Gather,  &
            EN_OR_Gather,   &
            EN_SUM_Scatter, &
            EN_MIN_Scatter, &
            EN_MAX_Scatter, &
            EN_OR_Scatter,  &
            EE_GS_INIT,     &
            EN_GS_INIT,     &
            NN_GS_INIT,     &
            Gather_BoundaryData, &
            NN_Gather_BoundaryData

END MODULE GS_MODULE
