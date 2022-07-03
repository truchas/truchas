#include "f90_assert.fpp"

module legacy_matl_api

  use matl_module
  use matl_utilities
  implicit none
  private

  public :: nmat, gather_vof, slot_resize,  matl_init, matl_free ! from matl_module
  public :: matl_get_vof, matl_get_cell_vof, define_matl, read_matl_data ! from matl_utilities

end module legacy_matl_api
