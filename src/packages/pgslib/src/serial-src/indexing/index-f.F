!=====================================================================
! Serial F90 routines to support PGSLib_Index_Table_MODULE
! These routines have the same calling convention as
! the parallal routines, but they are for serial use only.
!=====================================================================

! $Id: index-f.F,v 1.1.1.1 2000/10/11 22:44:23 ferrell Exp $

subroutine PGSLib_Init_Access_Table_C(Access_Table)
  implicit none 
  integer, dimension(*) :: Access_Table  

  Access_Table(2) = 1
  Access_Table(3) = 0

  RETURN
END subroutine PGSLib_Init_Access_Table_C

subroutine PGSLib_Free_Access_Table_C(Access_Table)
  implicit none 
  integer, dimension(*) :: Access_Table  

  Access_Table(2) = 0
  Access_Table(3) = 0

  RETURN
END subroutine PGSLib_Free_Access_Table_C

subroutine PGSLib_Add_Item_To_Table_C(Item, PE, Access_Table, ierror)
  implicit none
  integer               :: Item, PE
  integer, dimension(*) :: Access_Table
  integer               :: ierror

  IF (Access_Table(2) == 1) Then
     ierror = 0
  ELSE
     ierror = -1
  ENDIF

  RETURN
END subroutine PGSLib_Add_Item_To_Table_C

subroutine PGSLib_Count_Items_In_Table_C(Count, Access_Table)
  implicit none 
  integer               ::   Count
  integer, dimension(*) :: Access_Table

  IF (Access_Table(2) == 1) Then
     Count = 0
  ELSE
     Count = -1
  ENDIF

  RETURN
END subroutine PGSLib_Count_Items_In_Table_C

subroutine PGSLib_Items_From_Table_C(Items, PEs, Count, Access_Table, ierror)
  implicit none 
  integer               :: Count, ierror
  integer, DIMENSION(Count):: Items, PEs
  integer, dimension(*) :: Access_Table

  IF (Access_Table(2) == 1) Then
     IF(SIZE(Items,1) .GT. 0) Items = 0
     IF(SIZE(PEs,1) .GT. 0) PEs = 0
     ierror = 0
     Access_Table(3) = 1
  ELSE
     Items = -1
     PEs   = -1
     ierror = -1
  ENDIF


  RETURN
END subroutine PGSLib_Items_From_Table_C

subroutine PGSLib_Item_Index_From_Table_C(Index, Item, PE, Access_Table)
  implicit none 
  integer  :: Index
  integer  :: Item, PE
  integer, dimension(*) :: Access_Table

  IF ((Access_Table(2) == 1)  .AND. (Access_Table(3) == 1)) Then
     Index = 0
  ELSE
     Index = -1
  ENDIF

  RETURN
END subroutine PGSLib_Item_Index_From_Table_C

