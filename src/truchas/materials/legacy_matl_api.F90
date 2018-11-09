#include "f90_assert.fpp"

module legacy_matl_api

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use matl_module, only: nmat ! NEEDS TO BE MOVED ELSEWHERE
  use legacy_matl_adapter0_type
  use legacy_matl_adapter1_type
  use legacy_matl_adapter2_type
  implicit none
  private

  public :: nmat
  public :: gather_vof, matl_get_vof, matl_set_vof, update_matl, matl_get_cell_vof, read_matl_data
  public :: matl_redef, matl_enddef, matl_add_cell_mat, matl_set_cell_vof, matl_normalize_vof
  public :: matl_get_solid

  ! Different implementations of "matl"
  !type(legacy_matl_adapter0) :: this
  !type(legacy_matl_adapter1) :: this
  type(legacy_matl_adapter2) :: this

contains

  subroutine matl_redef
    call this%redef
  end subroutine

  subroutine matl_enddef
    call this%enddef
  end subroutine

  subroutine matl_add_cell_mat(cellid, matid)
    integer, intent(in) :: cellid, matid
    call this%add_cell_mat(cellid, matid)
  end subroutine

  subroutine gather_vof(m, vof)
    integer, intent(in) :: m
    real(r8), intent(out) :: vof(:)
    call this%gather_vof(m, vof)
  end subroutine

  subroutine matl_get_vof(vof)
    real(r8), intent(out) :: vof(:,:)
    call this%get_vof(vof)
  end subroutine

  subroutine matl_set_vof(vof)
    real(r8), intent(in) :: vof(:,:)
    call this%set_vof(vof)
  end subroutine

  subroutine update_matl(vof)
    real(r8), intent(in) :: vof(:,:)
    call this%update_vof(vof)
  end subroutine

  subroutine read_matl_data(unit, version)
    integer, intent(in) :: unit, version
    call this%read_data(unit, version)
  end subroutine

  subroutine matl_get_cell_vof(n, vof)
    integer, intent(in) :: n
    real(r8), intent(out) :: vof(:)
    call this%get_cell_vof(n, vof)
  end subroutine

  subroutine matl_set_cell_vof(n, m, vfrac)
    integer, intent(in) :: n, m
    real(r8), intent(in) :: vfrac
    call this%set_cell_vof(n, m, vfrac)
  end subroutine

  subroutine matl_normalize_vof
    call this%normalize_vof
  end subroutine

  subroutine matl_get_solid(mask)
    logical, intent(out) :: mask(:)
    call this%get_solid_mask(mask)
  end subroutine

end module
