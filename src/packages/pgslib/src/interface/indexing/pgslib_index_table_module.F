! Routines for constructing/maintaing offPE_Access_Tables

! $Id: pgslib_index_table_module.F,v 1.3 2002/09/20 19:15:10 lally Exp $

MODULE PGSLib_Index_Table_MODULE
  use pgslib_stats
  use pgslib_timing_module
  IMPLICIT NONE 
  SAVE
  PRIVATE
  PUBLIC :: PGSLib_Init_offPE_Access_Table
  PUBLIC :: PGSLib_Free_offPE_Access_Table
  PUBLIC :: PGSLib_Add_Item_To_Table
  PUBLIC :: PGSLib_Count_Items_In_Table
  PUBLIC :: PGSLib_Items_From_Table
  PUBLIC :: PGSLib_Item_Index_From_Table
  

  INTERFACE PGSLib_Init_offPE_Access_Table
     MODULE PROCEDURE PGSLib_Init_Access_Table_F
  END INTERFACE

  INTERFACE PGSLib_Free_offPE_Access_Table
     MODULE PROCEDURE PGSLib_Free_Access_Table_F
  END INTERFACE

  INTERFACE PGSLib_Add_Item_To_Table
     MODULE PROCEDURE PGSLib_Add_Item_To_Table_F
  END INTERFACE

  INTERFACE PGSLib_Count_Items_In_Table
     MODULE PROCEDURE PGSLib_Count_Items_In_Table_F
  END INTERFACE

  INTERFACE PGSLib_Items_From_Table
     MODULE PROCEDURE PGSLib_Items_From_Table_F
  END INTERFACE

  INTERFACE PGSLib_Item_Index_From_Table
     MODULE PROCEDURE PGSLib_Item_Index_From_Table_F
  END INTERFACE
  ! Interfaces to the C support routines
  INTERFACE 
     subroutine PGSLib_Init_Access_Table_C(C_Table)
       implicit none 
       integer, DIMENSION(*) :: C_Table  
     END subroutine PGSLib_Init_Access_Table_C

     subroutine PGSLib_Free_Access_Table_C(C_Table)
       implicit none 
       integer, DIMENSION(*) :: C_Table  
     END subroutine PGSLib_Free_Access_Table_C

     subroutine PGSLib_Add_Item_To_Table_C(Item, PE, C_Table, ierror)
       implicit none
       integer, intent(IN):: Item, PE
       integer, DIMENSION(*) :: C_Table
       integer, intent(OUT):: ierror
     END subroutine PGSLib_Add_Item_To_Table_C

     subroutine PGSLib_Count_Items_In_Table_C(Count, C_Table)
       implicit none 
       integer, intent(OUT)::   Count
       integer, DIMENSION(*) :: C_Table
     END subroutine PGSLib_Count_Items_In_Table_C

     subroutine PGSLib_Items_From_Table_C(Items, PEs, Count, C_Table, ierror)
       implicit none 
       integer, intent(IN) :: Count
       integer, intent(OUT), DIMENSION(Count):: Items, PEs
       integer, DIMENSION(*) :: C_Table
       integer, intent(OUT)              :: ierror
     END subroutine PGSLib_Items_From_Table_C

     subroutine PGSLib_Item_Index_From_Table_C(Index, Item, PE, C_Table)
       implicit none 
       integer, intent(OUT)  :: Index
       integer, intent( IN)  :: Item, PE
       integer, DIMENSION(*) :: C_Table
     END subroutine PGSLib_Item_Index_From_Table_C
  END INTERFACE

type ItemPEPair
   integer :: Item
   integer :: PE
end type ItemPEPair

INTERFACe OPERATOR(==)
   module procedure itempepair_eq
end INTERFACe

! This is to workaround a Sun compiler bug
! Otherwise would just use ==
INTERFACE OPERATOR(.EQItemPEPair.)
   module procedure itempepair_eq
end INTERFACE

INTERFACe OPERATOR(>)
   module procedure itempepair_gt
end INTERFACe

CONTAINS

  !=====================================================================
  !          PGSLib_Init_Access_Table_F(Table)
  ! PURPOSE:
  !          Initialize the internal hash table used for resolving off-PE
  !          accesses.
  !          If a hash table is already active, it is destroyed
  !          and replaced by a new one.
  !=====================================================================

  function PGSLib_Init_Access_Table_F() RESULT( Table )
    USE PGSLib_Type_MODULE
    implicit none 
    type (PGSLib_Access_Table), POINTER :: Table

    ! Allocate the Fortran table
    ALLOCATE(Table)
    ! Allocate enough space to hold a pointer
    ALLOCATE(Table%C_Table(4))

    call PGSLib_Init_Access_Table_C(Table%C_Table)

    RETURN
  END function PGSLib_Init_Access_Table_F

  !=====================================================================
  !          PGSLib_Free_Access_Table_F(Table)
  ! PURPOSE:
  !          Deallocate memeory used by access table
  !          If a hash table is already active, it is destroyed
  !          and replaced by a new one.
  !=====================================================================

  subroutine PGSLib_Free_Access_Table_F(Table)
    USE PGSLib_Type_MODULE
    implicit none 
    type (PGSLib_Access_Table), POINTER:: Table

    if (ASSOCIATED(Table%C_Table)) then
       call PGSLib_Free_Access_Table_C(Table%C_Table)
    end if
    DEALLOCATE(Table%C_Table)
    DEALLOCATE(Table)

    RETURN
  END subroutine PGSLib_Free_Access_Table_F

  !=====================================================================
  !         PGSLib_Add_Item_To_Table(Item, PE, Table)
  ! PURPOSE:
  !         Put an item into the internal hash table.  PE
  !         is the hashing key.
  !         The table must have been initialized before this call.
  !         Each table may have only one table active at any one time.
  !         A return value /= 0 indiacates an error.
  !=====================================================================

  function PGSLib_Add_Item_To_Table_F(Item, PE, Table)
    USE PGSLib_Type_MODULE
    implicit none 
    integer (PGSLib_Int_Type)                :: PGSLib_Add_Item_To_Table_F
    integer (PGSLib_Int_Type), INTENT(IN)    :: Item, PE
    type    (PGSLib_Access_Table), INTENT(INOUT):: Table

    ! local variables
    integer (PGSLib_Int_Type) :: ierror

#ifdef USE_TIMERS_2
    call enter_routine(Add_Item_To_T_STATISTICS())
#endif
    ! C counts PEs from 0 to nPE-1, FORTRAN counts from 1.
    call PGSLib_Add_Item_To_Table_C(Item, PE-1, Table%C_table, ierror)

    PGSLib_Add_Item_To_Table_F = ierror

#ifdef USE_TIMERS_2
    call exit_routine(Add_Item_To_T_STATISTICS())
#endif
    RETURN
  END function PGSLib_Add_Item_To_Table_F

  !=====================================================================
  !          PGSLib_Count_Items_In_Table
  ! PURPOSE: 
  !          Count the number of items in the hash table.
  !          The return value is the total number of items in the table.
  !          The hash table must have been initialized before this call.
  !          The return value is the number of items in the table.  If the
  !          return count is <0 that signals an error.
  !=====================================================================

  function PGSLib_Count_Items_In_Table_F(Table)
    USE PGSLib_Type_MODULE
    implicit none 
    integer (PGSLib_Int_Type)                :: PGSLib_Count_Items_In_Table_F
    type    (PGSLib_Access_Table), INTENT(INOUT):: Table

    ! local variables
    integer (PGSLib_Int_Type) :: count

    call PGSLib_Count_Items_In_Table_C(Count, Table%C_Table)

    PGSLib_Count_Items_In_Table_F = count
    RETURN
  END function PGSLib_Count_Items_In_Table_F


  !=====================================================================
  !          PGSLib_Items_From_Table
  ! PURPOSE: 
  !          Extract all the items and their associated PEs from the
  !          offPE access table, an internal hash table.  The purpose of
  !          this call is to get all the items and their PEs into two
  !          buffers.  Each buffer must have length at least
  !          PGSLib_Count_Items_In_Table(Table).  The return value
  !          is <0 if there is an error.  The table must be initialized
  !          before this call.
  !          THIS CALL CHANGES THE CONTENTS OF THE TABLE.
  !          This call has the side-effect of preparing the table
  !          for calls to PGSLib_Item_Index_From_Table.
  !=====================================================================

  function PGSLib_Items_From_Table_F(Items, PEs, Table)
    USE PGSLib_Type_MODULE
    implicit none 
    integer (PGSLib_Int_Type)                :: PGSLib_Items_From_Table_F
    integer (PGSLib_Int_Type),  INTENT(OUT), DIMENSION(:) :: Items, PEs
    type    (PGSLib_Access_Table), INTENT(INOUT)             :: Table

    ! Local variables
    integer (PGSLib_Int_Type) :: ierror, count

    ierror = 0
    count = PGSLib_Count_Items_In_Table(Table)
    IF (Count >= 0) THEN
       IF ((Size(Items,1) >= Count) .AND. (Size(PEs  ,1) >= Count)) THEN
          call PGSLib_Items_From_Table_C(Items, PEs, Count, Table%C_Table, ierror)
          ! C counts from 0 to nPE-1, FORTRAN counts from 1, so need to increment
          PEs = PEs + 1
        ELSE
           ierror = -1
        ENDIF
     ELSE
        ierror = -1
     ENDIF

     PGSLib_Items_From_Table_F = ierror
    RETURN
  END function PGSLib_Items_From_Table_F

  function itempepair_eq(a,b)
    implicit none
    type (ItemPEPair), intent(in) :: a,b
    logical :: itempepair_eq
    ItemPEPair_eq = (a%item == b%item) .AND. (a%pe == b%pe)
  end function itempepair_eq

  function itempepair_gt(a,b)
    ! NOTE: We are assuming that lists have items in decreasing order 
    ! *BUT* PEs in increasing order.  Hencs, the reversed polarities below.
    implicit none
    type (ItemPEPair), intent(in) :: a,b
    logical :: itempepair_gt
    ItemPEPair_gt = (a%pe < b%pe) .OR. &
                    ( (a%pe == b%pe) .AND. (a%item > b%item))   
  end function itempepair_gt

  !=====================================================================
  !          PGSLib_Item_Index_From_Table
  ! PURPOSE: 
  !          Extract the index of Supplement buffer for the input item.
  !          This is used to setup redirection of the Requestor to point
  !          to the proper location in the Supplement buffer.
  !=====================================================================
  
  function PGSLib_Item_Index_From_Table_F(Item, Item_PE, ITEM_LIST, PE_List) RESULT(Index)
    use pgslib_error_module
    USE PGSLib_Type_MODULE
    use pgslib_utility_module
    implicit none 
    integer (PGSLib_Int_Type)                :: Index
    integer (PGSLib_Int_Type),  INTENT(IN)   :: Item, Item_PE
    integer (PGSLib_Int_Type),  INTENT(IN), dimension(:)   :: Item_List, PE_List
!    type    (PGSLib_Access_Table), INTENT(INOUT):: Table

    ! Local variables
    integer (PGSLib_Int_Type) :: Lo, Mid, Hi, Count
    type (ItemPEPair)         :: ItemPE
    logical :: Error
    character (LEN=1024) :: out_str

#ifdef USE_TIMERS_2
    call enter_routine(Item_Index_From_T_STATISTICS())
#endif
    
    ! C counts PEs from 0 to nPE-1
!    Call PGSLib_Item_Index_From_Table_C(Index, Item, PE-1, Table%C_Table)

    ! This version (in test) finds the item in the item list, then returns the 
    ! index

    Lo = 1
    Hi = SIZE(Item_List,1)
    Index = -1
    ! Count is simply to make sure we don''t loop for ever.
    ! Loop should terminate in Log(N) time.
    ! Binary search assumes that Item is in decreasing order by items, but increasing order by PE
    ItemPE = ItemPEPair(Item, Item_PE)

    do Count = 1, size(Item_List,1)
       Mid = (Hi - Lo)/2 + Lo
       if (ItemPE .EQItemPEPair. ItemPEPair(Item_list(Lo ), PE_list(Lo) )) then
          Index = Lo 
          Exit
       end if
       if (ItemPE .EQItemPEPair. ItemPEPair(Item_list(Mid ), PE_list(Mid) )) then
          Index = Mid 
          Exit
       end if
       if (ItemPE .EQItemPEPair. ItemPEPair(Item_list(Hi ), PE_list(Hi) )) then
          Index = Hi 
          Exit
       end if
       if ( itempe > ItemPEPair(Item_list(mid), PE_List(Mid)) ) then
          Hi = Mid
       else
          Lo = Mid
       end if
    end do


    ! if we did not find the item in the list, that is an error

    Error = .FALSE.
    if (( (Index < 1) .OR. (Index > SIZE(Item_List)) ) ) then
       Error = .TRUE.
    else
       if (Item /= Item_List(Index)) then
          Error = .TRUE.
       end if
    end if

    if (Error) then

       itempe = ItemPEPair(Item, Item_PE)
       do count = 1, SIZE(Item_List)
          if (itempe .EQItemPEPair. ItemPEPair(Item_List(Count), PE_list(Count))) exit
       end do

       write(out_str,'(" Could not find itempe (",i9,",",i5,") in list")') item, item_pe
       call pgslib_error(out_str, "Pgslib_Item_Index_From_Table")

       if (Count <= SIZE(Item_List)) then
          write(out_str,'(" Found index ", i9," but item is in slot ", i9)') index, Count
       else
          write(out_str,'(" Found index ", i9," but item is not in list.")') index
       end if
       call pgslib_error(out_str)

       write(out_str,'(" SIZE(Item_List) = ", i6, " Lo, Mid, Hi = ", i6,i6,i6)') SIZE(Item_List), Lo, Mid, Hi
       call pgslib_error(out_str)
       write(out_str,'(" Lo  = ", i8, " Item_List(Lo)  = ", i8)') Lo, Item_List(Lo)
       call pgslib_error(out_str)
       write(out_str,'(" Mid = ", i8, " Item_List(Mid) = ", i8)') Mid, Item_List(Mid)
       call pgslib_error(out_str)
       write(out_str,'(" Hi  = ", i8, " Item_List(Hi)  = ", i8)') Hi, Item_List(Hi)
       call pgslib_error(out_str)
       
       do count = 1, SIZE(ITem_List)
          write(*, 10) PGSLib_Inquire_ThisPE(), Count, Item_List(Count), PE_List(Count)
       end do

10     FORMAT(1x, " ", i5,": ", i8, i8, i8)

       call pgslib_fatal_error("Not continuing")
       index = -1
       
    end if



#ifdef USE_TIMERS_2
      call exit_routine(Item_Index_From_T_STATISTICS())
#endif

    RETURN
  END function PGSLib_Item_Index_From_Table_F
    

END MODULE PGSLib_Index_Table_MODULE
