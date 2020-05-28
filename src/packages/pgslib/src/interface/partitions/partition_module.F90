MODULE PARTITION_MODULE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PURPOSE
  !   Provide the data types and support routines for paritioned
  !   sets and graphs.
  !
  ! $Id: 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use Partition_Constants,        ONLY: GRAPH_HEAD, GRAPH_TAIL, &
                                        PARTITION_INVALID_PARTITIONS
  use Partition_Data_Types,       ONLY: A_SET,                  &
                                        A_SET_PARTITIONING, &
                                        PARTITION_TO_PROCESSOR_MAP, &
                                        INITIALIZE, &
                                        ALLOC, &
                                        FREE, &
                                        SET, &
                                        LOOKUP, &
                                        TRANSLATE, &
                                        Partition_Is_Available, &
                                        Get_Num_Partitions, &
                                        Get_Num_Partitions_Available, &
                                        Get_Start, &
                                        Get_End, &
                                        Get_Sizes, &
                                        Get_Start_Available, &
                                        Get_End_Available, &
                                        Get_Sizes_Available

  use graph_partition_data_types,  ONLY: &
                                         A_PARTITIONED_GRAPH, &
                                         INITIALIZE, &
                                         ALLOC, &
                                         FREE, &
                                         SET, &
                                         Get_Num_Edges_Available, &
                                         Get_Set_Partitioning, &
                                         Get_Head_Partitions, &
                                         Get_Tail_Partitions, &
                                         Get_Head_Local, &
                                         Get_Head_Available, &
                                         LOOKUP_HEAD, &
                                         LOOKUP_TAIL, &
                                         Head_Is_Local, &
                                         Head_Is_Available


END MODULE PARTITION_MODULE
