module volume_tracker_class

  use kinds, only: r8
  use unstr_mesh_type
  use parameter_list_type
  implicit none
  private

  type, abstract, public :: volume_tracker
  contains
    procedure(vt_init), deferred :: init
    procedure(vt_flux_volumes), deferred :: flux_volumes
  end type volume_tracker

  abstract interface
    subroutine vt_init(this, mesh, nrealfluid, nfluid, nmat, liq_matid, params)
      import :: volume_tracker, unstr_mesh, parameter_list
      class(volume_tracker), intent(out) :: this
      type(unstr_mesh), intent(in), target :: mesh
      integer, intent(in) :: nrealfluid, nfluid, nmat, liq_matid(:)
      type(parameter_list), intent(inout) :: params
    end subroutine vt_init

    subroutine vt_flux_volumes(this, vel, vof_n, vof, flux_vol, fluids, void, dt)
      import :: volume_tracker, r8
      class(volume_tracker), intent(inout) :: this
      real(r8), intent(in) :: vel(:), vof_n(:,:), dt
      real(r8), intent(out) :: flux_vol(:,:), vof(:,:)
      integer, intent(in) :: fluids, void
    end subroutine vt_flux_volumes
  end interface

end module volume_tracker_class
