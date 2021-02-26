!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module volume_tracker_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use parameter_list_type
  implicit none
  private

  type, abstract, public :: volume_tracker
  contains
    procedure(vt_init), deferred :: init
    procedure(vt_flux_volumes), deferred :: flux_volumes
    procedure(vt_set_inflow_material), deferred :: set_inflow_material
    procedure :: donor_volumes_view
    procedure :: acceptor_fractions_view
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

    subroutine vt_set_inflow_material(this, mat, faces)
      import :: volume_tracker
      class(volume_tracker), intent(inout) :: this
      integer, intent(in) :: mat, faces(:)
    end subroutine
  end interface

contains

  function donor_volumes_view(this) result(p)
    class(volume_tracker), intent(in) :: this
    real(r8), pointer :: p(:,:)

    p => null()
  end function donor_volumes_view

  function acceptor_fractions_view(this) result(p)
    class(volume_tracker), intent(in) :: this
    real(r8), pointer :: p(:)

    p => null()
  end function acceptor_fractions_view
end module volume_tracker_class
