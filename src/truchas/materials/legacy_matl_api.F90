module legacy_matl_api

  use matl_module
  use matl_utilities
  implicit none
  private

  public :: max_slots, mat_slot, mat_slot_new, nmat
  public :: matl, gather_vof, slot_decrease, slot_increase, slot_set
  public :: matl_get_vof, matl_set_vof, read_matl_data, matl_get_cell_vof, update_matl

end module
