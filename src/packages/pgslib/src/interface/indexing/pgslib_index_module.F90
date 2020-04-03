MODULE PGSLib_Index_MODULE
  USE PGSLib_Index_Table_MODULE, ONLY : PGSLib_Init_offPE_Access_Table, &
       &                                PGSLib_Free_offPE_Access_Table, &
       &                                PGSLib_Add_Item_To_Table,       &
       &                                PGSLib_Count_Items_In_Table,    &
       &                                PGSLib_Items_From_Table,        &
       &                                PGSLib_Item_Index_From_Table

  USE PGSLib_Index_GID_MODULE,   ONLY : PGSlib_Setup_GID,               &
       &                                PGSLib_Deallocate_GID,          &
       &                                PGSLib_Release_GID,             &
       &                                PGSLib_Init_GID,                &
       &                                PGSLib_Index_Local_From_Global, &
       &                                PGSLib_Index_Global_From_Local, &
       &                                PGSLib_PE_From_Global_Index
END MODULE PGSLib_Index_MODULE
