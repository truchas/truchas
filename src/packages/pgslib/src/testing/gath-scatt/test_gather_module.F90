MODULE Test_Gather_MODULE
  !======================================================================
  ! PURPOSE 
  !   Test the gather routines provided by PGSLib
  !   This file is meant to be called from a driver.
  !======================================================================

  ! $Id: test_gather_module.F,v 1.1.1.1 2000/10/11 22:44:27 ferrell Exp $

  USE PGSLib_MODULE
  USE Test_Globals_MODULE, ONLY: EN_Trace
  IMPLICIT NONE 
  SAVE
  PRIVATE
  PUBLIC:: Test_EN_Gather
  

  INTERFACE Test_EN_Gather
     MODULE PROCEDURE Test_EN_Gather_INT
     MODULE PROCEDURE Test_EN_Gather_REAL
     MODULE PROCEDURE Test_EN_Gather_DOUBLE
!$VECTOR$!     MODULE PROCEDURE Test_EN_Gather_V_V_INT
!$VECTOR$!     MODULE PROCEDURE Test_EN_Gather_V_V_REAL
!$VECTOR$!     MODULE PROCEDURE Test_EN_Gather_V_V_DOUBLE
  END INTERFACE

CONTAINS
  
  ! In the subroutines below, the trace is EN_Trace, which is in the Glboals module.
  ! This avoids having to pass it around.
  
!======================================================================
!          Test_EN_Gather_INT
! This subroutine gathers from nodes to elements.  It calls
!      PGSLib_gather for the PE<->PE communication.
!======================================================================
  
  subroutine Test_EN_Gather_INT(elem_dest, nodes_src, elements_nodes, MASK)
    implicit none
    
    integer (PGSLib_INT_TYPE), intent(INOUT),                                 dimension(:, :) :: elem_dest
    integer (PGSLib_INT_TYPE), intent(INOUT),            dimension(SIZE(elem_dest,1), SIZE(elem_dest,2)) :: elements_nodes
    logical (PGSLib_Log_Type), intent(IN),  OPTIONAL, dimension(SIZE(elem_dest,1), SIZE(elem_dest,2)) :: Mask
    integer (PGSLib_INT_TYPE), intent(IN),                            dimension(:)  :: nodes_src
    

    RETURN
  END subroutine Test_EN_Gather_INT
  
!======================================================================
!          Test_EN_Gather_REAL
! This subroutine gathers from nodes to elements.  It calls
!      PGSLib_gather for the PE<->PE communication.
!======================================================================
      
  subroutine Test_EN_Gather_REAL(elem_dest, nodes_src, elements_nodes, MASK)
    implicit none
    
    real   (PGSLib_Real_TYPE), intent(INOUT),                                 dimension(:, :) :: elem_dest
    integer (PGSLib_INT_TYPE), intent(INOUT),            dimension(SIZE(elem_dest,1), SIZE(elem_dest,2)) :: elements_nodes
    logical (PGSLib_Log_Type), intent(IN),  OPTIONAL, dimension(SIZE(elem_dest,1), SIZE(elem_dest,2)) :: Mask
    real   (PGSLib_Real_TYPE), intent(IN),                            dimension(:)  :: nodes_src
    
    call pgslib_gather(elem_dest, nodes_src, elements_nodes, TRACE=EN_Trace, MASK=MASK)

    RETURN
  END subroutine Test_EN_Gather_REAL
  
!======================================================================
!          Test_EN_Gather_DOUBLE
! This subroutine gathers from nodes to elements.  It calls
!      PGSLib_gather for the PE<->PE communication.
!======================================================================      

  subroutine Test_EN_Gather_DOUBLE(elem_dest, nodes_src, elements_nodes, MASK)
    implicit none
    
    real (PGSLib_Double_TYPE), intent(INOUT),                                 dimension(:, :) :: elem_dest
    integer (PGSLib_INT_TYPE), intent(INOUT),            dimension(SIZE(elem_dest,1), SIZE(elem_dest,2)) :: elements_nodes
    logical (PGSLib_Log_Type), intent(IN),  OPTIONAL, dimension(SIZE(elem_dest,1), SIZE(elem_dest,2)) :: Mask
    real (PGSLib_Double_TYPE), intent(IN),                            dimension(:)  :: nodes_src
    
    call pgslib_gather(elem_dest, nodes_src, elements_nodes, TRACE=EN_Trace, MASK=MASK)


    RETURN
  END subroutine Test_EN_Gather_DOUBLE
  
!$VECTOR$!!======================================================================
!$VECTOR$!!          Test_EN_Gather_V_V_INT
!$VECTOR$!! This subroutine gathers from nodes to elements.  It calls
!$VECTOR$!!      PGSLib_gather for the PE<->PE communication.
!$VECTOR$!!======================================================================
!$VECTOR$!  
!$VECTOR$!  subroutine Test_EN_Gather_V_V_INT(elem_dest, nodes_vec_src, &
!$VECTOR$!       &                      elements_nodes, elements_nodes_item, elements_nodes_flag, &
!$VECTOR$!       &                      MASK, SRC_MASK)
!$VECTOR$!    implicit none
!$VECTOR$!    
!$VECTOR$!    integer (PGSLib_INT_TYPE), intent(INOUT),           &
!$VECTOR$!         &                     dimension( :              , :                ) :: elem_dest
!$VECTOR$!    integer (PGSLib_INT_TYPE), intent(IN),            &
!$VECTOR$!         &                     dimension(SIZE(elem_dest,1), SIZE(elem_dest,2)) :: elements_nodes
!$VECTOR$!    integer (PGSLib_INT_TYPE), intent(IN),            &
!$VECTOR$!         &                     dimension(SIZE(elem_dest,1), SIZE(elem_dest,2)) :: elements_nodes_item
!$VECTOR$!    integer (PGSLib_Log_Type), intent(IN),            &
!$VECTOR$!         &                     dimension(SIZE(elem_dest,2)) :: elements_nodes_flag
!$VECTOR$!    logical (PGSLib_Log_Type), intent(IN),  OPTIONAL, &
!$VECTOR$!         &                     dimension(SIZE(elem_dest,1), SIZE(elem_dest,2)) :: Mask
!$VECTOR$!    integer (PGSLib_INT_TYPE), intent(IN),            &
!$VECTOR$!         &                     dimension( :               , :               )  :: nodes_vec_src
!$VECTOR$!    logical (PGSLib_INT_TYPE), intent(IN),  OPTIONAL, dimension(SIZE(nodes_vec_src,1))  :: SRC_MASK
!$VECTOR$!    
!$VECTOR$!    call pgslib_gather(elem_dest, nodes_vec_src, &
!$VECTOR$!         &               elements_nodes, elements_nodes_item, elements_nodes_flag, EN_Trace, MASK, SRC_MASK)
!$VECTOR$!
!$VECTOR$!    RETURN
!$VECTOR$!  END subroutine Test_EN_Gather_V_V_INT
!$VECTOR$!  
!$VECTOR$!!======================================================================
!$VECTOR$!!          Test_EN_Gather_V_V_REAL
!$VECTOR$!! This subroutine gathers from nodes to elements.  It calls
!$VECTOR$!!      PGSLib_gather for the PE<->PE communication.
!$VECTOR$!!======================================================================
!$VECTOR$!      
!$VECTOR$!  subroutine Test_EN_Gather_V_V_REAL(elem_dest, nodes_vec_src, &
!$VECTOR$!       &                      elements_nodes, elements_nodes_item,elements_nodes_flag, &
!$VECTOR$!       &                      MASK, SRC_MASK)
!$VECTOR$!    implicit none
!$VECTOR$!    
!$VECTOR$!    real   (PGSLib_Real_TYPE), intent(INOUT),           &
!$VECTOR$!         &                                            dimension(:, :) :: elem_dest
!$VECTOR$!    integer (PGSLib_INT_TYPE), intent(IN),            &
!$VECTOR$!         &                                            dimension(SIZE(elem_dest,1), SIZE(elem_dest,2)) :: elements_nodes
!$VECTOR$!    integer (PGSLib_INT_TYPE), intent(IN),            &
!$VECTOR$!         &                                            dimension(SIZE(elem_dest,1), SIZE(elem_dest,2)) :: elements_nodes_item
!$VECTOR$!    integer (PGSLib_Log_Type), intent(IN),            &
!$VECTOR$!         &                                            dimension(SIZE(elem_dest,2)) :: elements_nodes_flag
!$VECTOR$!    logical (PGSLib_Log_Type), intent(IN),  OPTIONAL, &
!$VECTOR$!         &                                            dimension(SIZE(elem_dest,1), SIZE(elem_dest,2)) :: Mask
!$VECTOR$!    real   (PGSLib_Real_TYPE), intent(IN),            &
!$VECTOR$!         &                                            dimension( :               , :               )  :: nodes_vec_src
!$VECTOR$!    logical (PGSLib_INT_TYPE), intent(IN),  OPTIONAL, dimension(SIZE(nodes_vec_src,1)              )  :: SRC_MASK
!$VECTOR$!    
!$VECTOR$!    call pgslib_gather(elem_dest, nodes_vec_src, &
!$VECTOR$!         &               elements_nodes, elements_nodes_item, elements_nodes_flag, EN_Trace, MASK, SRC_MASK)
!$VECTOR$!    
!$VECTOR$!
!$VECTOR$!    RETURN
!$VECTOR$!  END subroutine Test_EN_Gather_V_V_REAL
!$VECTOR$!  
!$VECTOR$!!======================================================================
!$VECTOR$!!          Test_EN_Gather_V_V_DOUBLE
!$VECTOR$!! This subroutine gathers from nodes to elements.  It calls
!$VECTOR$!!      PGSLib_gather for the PE<->PE communication.
!$VECTOR$!!======================================================================      
!$VECTOR$!
!$VECTOR$!  subroutine Test_EN_Gather_V_V_DOUBLE(elem_dest, nodes_vec_src, &
!$VECTOR$!                              elements_nodes, elements_nodes_item, elements_nodes_flag, &
!$VECTOR$!       &                      MASK, SRC_MASK)
!$VECTOR$!    implicit none
!$VECTOR$!    
!$VECTOR$!    real   (PGSLib_Double_TYPE), intent(INOUT),         &
!$VECTOR$!         &                                            dimension(:, :) :: elem_dest
!$VECTOR$!    integer (PGSLib_INT_TYPE), intent(IN),            &
!$VECTOR$!         &                                            dimension(SIZE(elem_dest,1), SIZE(elem_dest,2)) :: elements_nodes
!$VECTOR$!    integer (PGSLib_INT_TYPE), intent(IN),            &
!$VECTOR$!         &                                            dimension(SIZE(elem_dest,1), SIZE(elem_dest,2)) :: elements_nodes_item
!$VECTOR$!    integer (PGSLib_Log_Type), intent(IN),            &
!$VECTOR$!         &                                            dimension(SIZE(elem_dest,2)) :: elements_nodes_flag
!$VECTOR$!    logical (PGSLib_Log_Type), intent(IN),  OPTIONAL, &
!$VECTOR$!         &                                            dimension(SIZE(elem_dest,1), SIZE(elem_dest,2)) :: Mask
!$VECTOR$!    real   (PGSLib_Double_TYPE), intent(IN),          &
!$VECTOR$!         &                                            dimension( :               , :               )  :: nodes_vec_src
!$VECTOR$!    logical (PGSLib_INT_TYPE), intent(IN),  OPTIONAL, dimension(SIZE(nodes_vec_src,1)              )  :: SRC_MASK
!$VECTOR$!    
!$VECTOR$!    
!$VECTOR$!    call pgslib_gather(elem_dest, nodes_vec_src, &
!$VECTOR$!         &               elements_nodes, elements_nodes_item, elements_nodes_flag, EN_Trace, MASK, SRC_MASK)
!$VECTOR$!    
!$VECTOR$!    RETURN
!$VECTOR$!  END subroutine Test_EN_Gather_V_V_DOUBLE
  
END MODULE Test_Gather_MODULE
