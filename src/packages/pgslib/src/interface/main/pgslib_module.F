MODULE PGSLib_MODULE
  !======================================================================
  ! PURPOSE
  !   This is the master module for pgslib.  It does nothing more
  !   than USE the constitent modules.
  !======================================================================

  ! $Id: pgslib_module.F,v 1.6 2002/09/12 20:52:34 lally Exp $


  USE PGSLib_Type_MODULE
  USE PGSLib_Globals_MODULE
  USE PGSLib_IO_MODULE,         ONLY: PGSLib_BCast,            &
                                      PGSLib_Dist,             &
                                      PGSLib_Collate

  USE PGSLib_Reductions_MODULE,   only: PGSLib_Global_MINVAL, PGSLib_Global_MAXVAL, &
                                        PGSLib_Global_SUM,    PGSLib_Global_COUNT,  &
                                        PGSLib_Global_ALL,    PGSLib_Global_ANY,    &
                                        PGSLib_Global_DOT_PRODUCT,                  &
                                        PGSLib_Global_MINLOC, PGSLib_Global_MAXLOC
  
  USE PGSLib_Utility_MODULE,    ONLY: PGSLib_Initialize,       &
                                                   PGSLib_Finalize,         &
                                                   PGSLib_Inquire_nPE,      &
                                                   PGSLib_Inquire_thisPE,   &
                                                   PGSLib_Inquire_thisPE_Actual,   &
                                                   PGSLib_Inquire_IO_ROOT_PE,&
                                                   PGSLib_Inquire_IO_P,&
                                                   PGSLib_UseGlobalServices,&
                                                   PGSLib_Output,           &
                                                   PGSLib_Flush_Output,     &
                                                   PGSLib_Error,            &
                                                   PGSLib_Abort,            &
                                                   PGSLib_Fatal_Error,      &
                                                   PGSLib_Check_Error,      &
                                                   PGSLib_Barrier,          &
                                                   PGSLib_Start_Timer,      &
                                                   PGSLib_Stop_Timer,       &
                                                   PGSLib_Memory_Size,      &
                                                   PGSLib_CL_MAX_TOKEN_LENGTH


  USE PGSLib_Instrument
  USE PGSLib_Stats

  USE PGSLib_GS_MODULE,         ONLY: PGSLib_GS_Init_Trace,    &
                                      PGSLib_GS_Release_Trace, &
                                      PGSLib_Deallocate_Trace, &
                                      PGSLib_GS_Trace_Setup_P, &
                                      PGSLib_Setup_Trace,      &
                                      PGSLib_Scatter_Buffer,   &
                                      PGSLib_Gather_Buffer,    &
                                      PGSLib_Buffer_Gather

  USE PGSLib_User_GS_MODULE,    ONLY: PGSLib_Gather,&
                                    PGSLib_Scatter_SUM, PGSLib_Scatter_MAX, PGSLib_Scatter_MIN, &
                                    PGSLib_Scatter_OR, PGSLib_Scatter_AND

  USE PGSLib_Index_MODULE,      ONLY: PGSlib_Setup_GID


  USE PGSLib_Scan_MODULE,       ONLY: PGSLib_SUM_PREFIX,        &
                                      PGSLib_SUM_SUFFIX,        &
                                      PGSLib_PARITY_PREFIX,     &
                                      PGSLib_PARITY_SUFFIX

  USE PGSLib_Sort_MODULE,       ONLY: PGSLib_GRADE_UP,          &
                                      PGSLib_GRADE_UP_LOCAL

  USE PGSLib_Misc_MODULE,       ONLY: PGSLib_Permute,             &
                                      PGSLib_Redistribute,        &
                                      PGSLib_Global_Pack,         &
                                      PGSLib_Global_CShift,       &
                                      PGSLib_Global_EOShift,      &
                                      PGSLib_Block_Size

  USE Partition_Module,           ONLY: A_SET,                  &
                                        A_SET_PARTITIONING, &
                                        PARTITION_TO_PROCESSOR_MAP, &
                                        INITIALIZE, &
                                        ALLOC, &
                                        FREE, &
                                        SET, &
                                        LOOKUP, &
                                        TRANSLATE, &
                                        Get_Num_Partitions, &
                                        Get_Num_Partitions_Available, &
                                        Get_Start, &
                                        Get_End, &
                                        Get_Sizes, &
                                        Get_Start_Available, &
                                        Get_End_Available, &
                                        Get_Sizes_Available, &
                                        A_PARTITIONED_GRAPH, &
                                        Get_Num_Edges_Available, &
                                        GRAPH_HEAD, GRAPH_TAIL, &
                                        Get_Set_Partitioning, &
                                        PARTITION_INVALID_PARTITIONS, &
                                        Partition_Is_Available, &
                                        Get_Head_Partitions, &
                                        Get_Tail_Partitions, &
                                        Get_Head_Local, &
                                        Get_Head_Available, &
                                        LOOKUP_TAIL, &
                                        LOOKUP_HEAD, &
                                        Head_Is_Local, &
                                        Head_Is_Available
  

END MODULE PGSLib_MODULE
