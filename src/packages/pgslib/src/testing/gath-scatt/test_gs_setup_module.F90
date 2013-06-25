MODULE TEST_GS_SETUP_MODULE
  ! Initilize the communication patterns needed.

  ! $Id: test_gs_setup_module.F,v 1.1.1.1 2000/10/11 22:44:27 ferrell Exp $

!  USE PGSLib_MODULE
!  USE TEST_Globals_MODULE
  IMPLICIT NONE
  SAVE
  PRIVATE
  PUBLIC :: TEST_EN_GS_Setup

CONTAINS

!======================================================================
!          TEST_EN_GS_SETUP
! PURPOSE:
!          Establish comm pattern for element<->node communication
!          The trace, EN_Trace, is in TEST_Globals_MODULE
!======================================================================

  subroutine TEST_EN_GS_SETUP(elements_nodes, NNodes, Mask)
    USE PGSLib_MODULE
    USE TEST_Globals_MODULE
    implicit none
    integer, INTENT(INOUT), dimension(:,:):: elements_nodes
    logical, INTENT(IN   ), optional, dimension(:,:):: mask
    integer, INTENT(INOUT)                :: NNodes

    NULLIFY(EN_Trace )
    EN_Trace => PGSLib_Setup_Trace(elements_nodes, NNodes, MASK=MASK)
    
    RETURN
  END subroutine TEST_EN_GS_SETUP
END MODULE TEST_GS_SETUP_MODULE
    
