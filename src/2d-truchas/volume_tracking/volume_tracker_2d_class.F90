!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module volume_tracker_2d_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_2d_mesh_type
  implicit none
  private

  type, abstract, public :: volume_tracker_2d
  contains
    procedure(vt_init), deferred :: init
    procedure(vt_flux_volumes), deferred :: flux_volumes
    procedure(vt_set_inflow_material), deferred :: set_inflow_material
  end type volume_tracker_2d

  abstract interface
    subroutine vt_init(this, mesh, nrealfluid, nfluid, nmat)
      import :: volume_tracker_2d, unstr_2d_mesh
      class(volume_tracker_2d), intent(out) :: this
      type(unstr_2d_mesh), intent(in), target :: mesh
      integer, intent(in) :: nrealfluid, nfluid, nmat
    end subroutine vt_init

    subroutine vt_flux_volumes(this, vel, vof_n, vof, flux_vol, int_normal, fluids, &
        void, dt)
      import :: volume_tracker_2d, r8
      class(volume_tracker_2d), intent(inout) :: this
      real(r8), intent(in) :: vel(:), vof_n(:,:), dt
      real(r8), intent(out) :: flux_vol(:,:), vof(:,:), int_normal(:,:,:)
      integer, intent(in) :: fluids, void
    end subroutine vt_flux_volumes

    subroutine vt_set_inflow_material(this, mat, faces)
      import :: volume_tracker_2d
      class(volume_tracker_2d), intent(inout) :: this
      integer, intent(in) :: mat, faces(:)
    end subroutine
  end interface

end module volume_tracker_2d_class
