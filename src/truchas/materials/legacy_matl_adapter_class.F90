module legacy_matl_adapter_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: legacy_matl_adapter
  contains
    procedure(redef), deferred :: redef
    procedure(redef), deferred :: enddef
    procedure(add_cell_mat), deferred :: add_cell_mat
    procedure(gather_vof), deferred :: gather_vof
    procedure(get_vof), deferred :: get_vof
    procedure(set_vof), deferred :: set_vof
    procedure(update_vof), deferred :: update_vof
    procedure(read_data), deferred :: read_data
    procedure(get_cell_vof), deferred :: get_cell_vof
    procedure(set_cell_vof), deferred :: set_cell_vof
    procedure(normalize_vof), deferred :: normalize_vof
    procedure(get_solid_mask), deferred :: get_solid_mask
  end type legacy_matl_adapter

  abstract interface
    subroutine redef(this)
      import legacy_matl_adapter
      class(legacy_matl_adapter), intent(inout) :: this
    end subroutine
    subroutine add_cell_mat(this, cellid, matid)
      import legacy_matl_adapter
      class(legacy_matl_adapter), intent(inout) :: this
      integer, intent(in) :: cellid, matid
    end subroutine
    subroutine gather_vof(this, m, vof)
      import legacy_matl_adapter, r8
      class(legacy_matl_adapter), intent(in) :: this
      integer, intent(in) :: m
      real(r8), intent(out) :: vof(:)
    end subroutine
    subroutine get_vof(this, vof)
      import legacy_matl_adapter, r8
      class(legacy_matl_adapter), intent(in) :: this
      real(r8), intent(out) :: vof(:,:)
    end subroutine
    subroutine set_vof(this, vof)
      import legacy_matl_adapter, r8
      class(legacy_matl_adapter), intent(inout) :: this
      real(r8), intent(in) :: vof(:,:)
    end subroutine
    subroutine update_vof(this, vof)
      import legacy_matl_adapter, r8
      class(legacy_matl_adapter), intent(inout) :: this
      real(r8), intent(in) :: vof(:,:)
    end subroutine
    subroutine read_data(this, unit, version)
      import legacy_matl_adapter
      class(legacy_matl_adapter), intent(inout) :: this
      integer, intent(in) :: unit, version
    end subroutine
    subroutine get_cell_vof(this, n, vof)
      import legacy_matl_adapter, r8
      class(legacy_matl_adapter), intent(in) :: this
      integer, intent(in) :: n
      real(r8), intent(out) :: vof(:)
    end subroutine
    subroutine set_cell_vof(this, n, m, vfrac)
      import legacy_matl_adapter, r8
      class(legacy_matl_adapter), intent(inout) :: this
      integer, intent(in) :: n, m
      real(r8), intent(in) :: vfrac
    end subroutine
    subroutine normalize_vof(this)
      import legacy_matl_adapter
      class(legacy_matl_adapter), intent(inout) :: this
    end subroutine
    subroutine get_solid_mask(this, mask)
      import legacy_matl_adapter
      class(legacy_matl_adapter), intent(in) :: this
      logical, intent(out) :: mask(:)
    end subroutine
  end interface

end module legacy_matl_adapter_class
