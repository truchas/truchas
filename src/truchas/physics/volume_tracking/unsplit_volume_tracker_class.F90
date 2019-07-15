!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module unsplit_volume_tracker_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use parameter_list_type
  use cell_tagged_mm_volumes_type, only : cell_tagged_mm_volumes    
  implicit none
  private

  type, abstract, public :: unsplit_volume_tracker
    real(r8), public :: cutoff
  contains
    procedure(vt_init), deferred :: init
    procedure(vt_flux_volumes), deferred :: flux_volumes
    procedure(vt_set_inflow_material), deferred :: set_inflow_material
    procedure(vt_write_interface), deferred :: write_interface
  end type unsplit_volume_tracker

  abstract interface
    subroutine vt_init(this, mesh, nrealfluid, nfluid, nmat, liq_matid, params)
      import :: unsplit_volume_tracker, unstr_mesh, parameter_list
      class(unsplit_volume_tracker), intent(out) :: this
      type(unstr_mesh), intent(in), target :: mesh
      integer, intent(in) :: nrealfluid, nfluid, nmat, liq_matid(:)
      type(parameter_list), intent(inout) :: params
    end subroutine vt_init

    subroutine vt_flux_volumes(this, vel, vel_cc, vof_n, vof, flux_vol, fluids, void, dt, a_interface_band)
      import :: unsplit_volume_tracker, r8, cell_tagged_mm_volumes
      class(unsplit_volume_tracker), intent(inout) :: this
      real(r8), intent(in) :: vel(:), vel_cc(:,:), vof_n(:,:), dt
      type(cell_tagged_mm_volumes), intent(out) :: flux_vol(:)
      real(r8), intent(out) ::  vof(:,:)
      integer, intent(in) :: fluids, void
      integer, intent(in) :: a_interface_band(:)
    end subroutine vt_flux_volumes

    subroutine vt_set_inflow_material(this, mat, faces)
      import :: unsplit_volume_tracker
      class(unsplit_volume_tracker), intent(inout) :: this
      integer, intent(in) :: mat, faces(:)
    end subroutine vt_set_inflow_material

    subroutine vt_write_interface(this, t, dt, cycle_number)
      import unsplit_volume_tracker, r8
      class(unsplit_volume_tracker), intent(in) :: this
      real(r8), intent(in) :: t
      real(r8), intent(in) :: dt
      integer, intent(in) :: cycle_number
    end subroutine vt_write_interface
  end interface

end module unsplit_volume_tracker_class
