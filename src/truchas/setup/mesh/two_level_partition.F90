MODULE Two_Level_Partition
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Purpose:
  !    Data for two-level partitioning, particularly on-processor 
  !    partitioning.
  !
  ! Author: Robert Ferrell (ferrell@cpca.com)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use pgslib_module, only: A_Set, A_SET_PARTITIONING, A_PARTITIONED_GRAPH

  implicit none

  PRIVATE
  ! The sets which these are derived from
  type(A_SET),               SAVE, PUBLIC :: Set_Of_Cells
  type(A_SET),               SAVE, PUBLIC :: Set_Of_Nodes
  type(A_SET),               SAVE, PUBLIC :: Set_Of_Cell_Cell_Edges
  ! The two-level partitions
  type (A_SET_PARTITIONING), SAVE, PUBLIC :: Cell_Two_Level_Partitioning
!  type (A_SET_PARTITIONING), SAVE, PUBLIC :: Node_Two_Level_Partitioning
  type (A_SET_PARTITIONING), SAVE, PUBLIC :: Cell_Cell_Two_Level_Edges_Part

  ! And the graph partitions
  type (A_PARTITIONED_GRAPH),SAVE, PUBLIC :: Cell_Cell_Two_Level_Partitioned

END MODULE Two_Level_Partition
