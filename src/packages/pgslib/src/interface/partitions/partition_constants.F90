MODULE PARTITION_CONSTANTS
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! $Id: partition_constants.F,v 1.1 2000/12/14 17:45:59 ferrell Exp $
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE
  PRIVATE
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!! Constant declarations
  integer, PUBLIC, parameter :: PARTITION_INVALID_SIZE       = -1
  integer, PUBLIC, parameter :: PARTITION_INVALID_PARTITIONS = -2
  integer, PUBLIC, parameter :: PARTITION_INVALID_PROCESSOR  = -3

  !!!!!!!!!! Identifies for graphs
  integer, PUBLIC, parameter :: GRAPH_HEAD = 2
  integer, PUBLIC, parameter :: GRAPH_TAIL = 1
  
  

END MODULE PARTITION_CONSTANTS
