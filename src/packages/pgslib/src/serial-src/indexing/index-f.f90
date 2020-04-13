!!
!! index-f.F
!!
!! Serial-only dummy routines that correspond to the MPI parallel
!! C functions from par-src/indexing/index-c.c.
!!
!! This is a modern rewrite of the original code by Robert Ferrell
!! using the C interoperability features of Fortran 2003.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! November 2014
!!

subroutine pgslib_init_access_table_c (Access_Table) &
    bind(c, name="pgslib_init_access_table_c")
  use,intrinsic :: iso_c_binding, only: c_ptr
  type(c_ptr) :: Access_Table
end subroutine

subroutine pgslib_free_access_table_c (Access_Table) &
    bind(c, name="pgslib_free_access_table_c")
  use,intrinsic :: iso_c_binding, only: c_ptr
  type(c_ptr) :: Access_Table
end subroutine

subroutine pgslib_add_item_to_table_c (Item, PE, Access_Table, ierror) &
    bind(c, name="pgslib_add_item_to_table_c")
  use,intrinsic :: iso_c_binding, only: c_ptr, c_int
  integer(c_int), intent(in) :: Item, PE
  type(c_ptr) :: Access_Table
  integer(c_int), intent(out) :: ierror
  ierror = 0
end subroutine

subroutine PGSLib_Count_Items_In_Table_C (Count, Access_Table) &
    bind(c, name="pgslib_count_items_in_table_c")
  use,intrinsic :: iso_c_binding, only: c_ptr, c_int
  integer(c_int), intent(out) :: Count
  type(c_ptr) :: Access_Table
  Count = 0
end subroutine

subroutine PGSLib_Items_From_Table_C (Items, PEs, Count, Access_Table, ierror) &
    bind(c, name="pgslib_items_from_table_c")
  use,intrinsic :: iso_c_binding, only: c_ptr, c_int
  integer(c_int), intent(in)  :: Count
  integer(c_int), intent(out) :: Items(Count), PEs(Count), ierror
  type(c_ptr) :: Access_Table
  Items = 0
  PEs = 0
  ierror = 0
end subroutine

subroutine PGSLib_Item_Index_From_Table_C (Index, Item, PE, Access_Table) &
    bind(c, name="pgslib_item_index_from_table_c")
  use,intrinsic :: iso_c_binding, only: c_ptr, c_int
  integer(c_int), intent(out) :: Index
  integer(c_int), intent(in)  :: Item, PE
  type(c_ptr) :: Access_Table
  Index = 0
end subroutine
