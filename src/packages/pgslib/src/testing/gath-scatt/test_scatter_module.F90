MODULE TEST_Scatter_MODULE
  USE PGSLib_MODULE
  USE TEST_Globals_MODULE, ONLY: EN_Trace
  IMPLICIT NONE 
  SAVE
  PRIVATE
  PUBLIC :: TEST_EN_Scatter_SUM, TEST_EN_Scatter_OR, TEST_EN_Scatter_AND
  PUBLIC :: TEST_EN_Scatter_MAX

  ! $Id: test_scatter_module.F,v 1.1.1.1 2000/10/11 22:44:27 ferrell Exp $

  INTERFACE TEST_EN_Scatter_SUM
     MODULE PROCEDURE TEST_EN_Scatter_SUM_INT
     MODULE PROCEDURE TEST_EN_Scatter_SUM_REAL
     MODULE PROCEDURE TEST_EN_Scatter_SUM_DOUBLE
  END INTERFACE

  INTERFACE TEST_EN_Scatter_MAX
     MODULE PROCEDURE TEST_EN_Scatter_MAX_DOUBLE
  END INTERFACE

  INTERFACE TEST_EN_Scatter_OR
     MODULE PROCEDURE TEST_EN_Scatter_OR_LOG
  END INTERFACE

  INTERFACE TEST_EN_Scatter_AND
     MODULE PROCEDURE TEST_EN_Scatter_AND_LOG
  END INTERFACE

CONTAINS

! In the subroutines below, the trace is EN_Trace, which is in the Glboals module.
! This avoids having to pass it around.

!======================================================================
!          TEST_EN_Scatter_SUM_INT
! This subroutine scatter-adds from elements to nodes.  It calls
!      PGSLib_scatter for the PE<->PE communication.
!======================================================================      
  subroutine TEST_EN_Scatter_SUM_INT(nodes_dest, elem_src, elements_nodes, MASK)
    implicit none
    
    integer (PGSLib_INT_TYPE), intent(INOUT),                            dimension(:) :: nodes_dest
    integer (PGSLib_INT_TYPE), intent(IN),                                          dimension(:, :) :: elem_src
    integer (PGSLib_INT_TYPE), intent(INOUT),            dimension(SIZE(elem_src,1), SIZE(elem_src,2)) :: elements_nodes
    logical (PGSLib_Log_Type), intent(IN),  OPTIONAL, dimension(SIZE(elem_src,1), SIZE(elem_src,2)) :: MASK

    call pgslib_scatter_SUM(nodes_dest, elem_src, elements_nodes, TRACE=EN_Trace, MASK=MASK)

    RETURN
  end subroutine TEST_EN_Scatter_SUM_INT
  
!======================================================================
!          TEST_EN_Scatter_SUM_REAL
! This subroutine scatter-adds from elements to nodes.  It calls
!      PGSLib_scatter for the PE<->PE communication.
!======================================================================
  
  subroutine TEST_EN_Scatter_SUM_REAL(nodes_dest, elem_src, elements_nodes, MASK)
    implicit none

    real   (PGSLib_Real_Type), intent(INOUT),                            dimension(:) :: nodes_dest
    real   (PGSLib_Real_Type), intent(IN),                                          dimension(:, :) :: elem_src
    integer (PGSLib_INT_TYPE), intent(INOUT),            dimension(SIZE(elem_src,1), SIZE(elem_src,2)) :: elements_nodes
    logical (PGSLib_Log_Type), intent(IN),  OPTIONAL, dimension(SIZE(elem_src,1), SIZE(elem_src,2)) :: MASK

    
    call pgslib_scatter_SUM(nodes_dest, elem_src, elements_nodes, TRACE=EN_Trace, MASK=MASK)

    RETURN
  end subroutine TEST_EN_Scatter_SUM_REAL

!======================================================================
!          TEST_EN_Scatter_SUM_DOUBLE
! This subroutine scatter-adds from elements to nodes.  It calls
!      PGSLib_scatter for the PE<->PE communication.
!======================================================================
      
  subroutine TEST_EN_Scatter_SUM_DOUBLE(nodes_dest, elem_src, elements_nodes, MASK)
    implicit none

    real   (PGSLib_Double_Type), intent(INOUT),                            dimension(:) :: nodes_dest
    real   (PGSLib_Double_Type), intent(IN),                                          dimension(:, :) :: elem_src
    integer (PGSLib_INT_TYPE), intent(INOUT),            dimension(SIZE(elem_src,1), SIZE(elem_src,2)) :: elements_nodes
    logical (PGSLib_Log_Type), intent(IN),  OPTIONAL, dimension(SIZE(elem_src,1), SIZE(elem_src,2)) :: MASK

    call pgslib_scatter_SUM(nodes_dest, elem_src, elements_nodes, TRACE=EN_Trace, MASK=MASK)

    RETURN
  end subroutine TEST_EN_Scatter_SUM_DOUBLE

!======================================================================
!          TEST_EN_Scatter_MAX_DOUBLE
! This subroutine scatter-adds from elements to nodes.  It calls
!      PGSLib_scatter for the PE<->PE communication.
!======================================================================
      
  subroutine TEST_EN_Scatter_MAX_DOUBLE(nodes_dest, elem_src, elements_nodes, MASK)
    implicit none

    real   (PGSLib_Double_Type), intent(INOUT),                            dimension(:) :: nodes_dest
    real   (PGSLib_Double_Type), intent(IN),                                          dimension(:, :) :: elem_src
    integer (PGSLib_INT_TYPE), intent(INOUT),            dimension(SIZE(elem_src,1), SIZE(elem_src,2)) :: elements_nodes
    logical (PGSLib_Log_Type), intent(IN),  OPTIONAL, dimension(SIZE(elem_src,1), SIZE(elem_src,2)) :: MASK

    call pgslib_scatter_MAX(nodes_dest, elem_src, elements_nodes, TRACE=EN_Trace, MASK=MASK)

    RETURN
  end subroutine TEST_EN_Scatter_MAX_DOUBLE

!======================================================================
!          TEST_EN_Scatter_OR_LOG
! This subroutine scatter-ORs from elements to nodes.  It calls
!      PGSLib_scatter for the PE<->PE communication.
!======================================================================
      
  subroutine TEST_EN_Scatter_OR_LOG(nodes_dest, elem_src, elements_nodes, MASK)
    implicit none

    logical (PGSLib_Log_Type), intent(INOUT),                            dimension(:) :: nodes_dest
    logical (PGSLib_Log_Type), intent(IN),                                          dimension(:, :) :: elem_src
    integer (PGSLib_INT_TYPE), intent(INOUT),            dimension(SIZE(elem_src,1), SIZE(elem_src,2)) :: elements_nodes
    logical (PGSLib_Log_Type), intent(IN),  OPTIONAL, dimension(SIZE(elem_src,1), SIZE(elem_src,2)) :: MASK

    call pgslib_scatter_OR(nodes_dest, elem_src, elements_nodes, TRACE=EN_Trace, MASK=MASK)
    
    RETURN
  end subroutine TEST_EN_Scatter_OR_LOG

!======================================================================
!          TEST_EN_Scatter_AND_LOG
! This subroutine scatter-ANDs from elements to nodes.  It calls
!      PGSLib_scatter for the PE<->PE communication.
!======================================================================
      
  subroutine TEST_EN_Scatter_AND_LOG(nodes_dest, elem_src, elements_nodes, MASK)
    implicit none

    logical (PGSLib_Log_Type), intent(INOUT),                            dimension(:) :: nodes_dest
    logical (PGSLib_Log_Type), intent(IN),                                          dimension(:, :) :: elem_src
    integer (PGSLib_INT_TYPE), intent(INOUT),            dimension(SIZE(elem_src,1), SIZE(elem_src,2)) :: elements_nodes
    logical (PGSLib_Log_Type), intent(IN),  OPTIONAL, dimension(SIZE(elem_src,1), SIZE(elem_src,2)) :: MASK

    call pgslib_scatter_and(nodes_dest, elem_src, elements_nodes, TRACE=EN_Trace, MASK=MASK)
    
    RETURN
  end subroutine TEST_EN_Scatter_AND_LOG

END MODULE TEST_Scatter_MODULE

