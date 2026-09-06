!!
!! Aditya K. Pandare <apandare@lanl.gov>, January 2020
!! SPDX-License-Identifier: BSD-3-Clause
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module t2d_volume_tracker_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use simulation_environment_type, only: simulation_environment
  use t2d_unstr_mesh_type
  implicit none
  private

  type, abstract, public :: t2d_volume_tracker
  contains
    procedure(vt_init), deferred :: init
    procedure(vt_flux_volumes), deferred :: flux_volumes
    procedure(vt_set_inflow_material), deferred :: set_inflow_material
  end type t2d_volume_tracker

  abstract interface
    subroutine vt_init(this, env, mesh, nrealfluid, nfluid, nmat, axisym, priority, cutoff)
      import :: t2d_volume_tracker, t2d_unstr_mesh, r8
      import :: simulation_environment
      class(t2d_volume_tracker), intent(out) :: this
      type(simulation_environment), intent(in) :: env
      type(t2d_unstr_mesh), intent(in), target :: mesh
      integer, intent(in) :: nrealfluid, nfluid, nmat
      logical, intent(in) :: axisym
      integer, intent(in) :: priority(:)
      real(r8), optional, intent(in) :: cutoff
    end subroutine vt_init

    subroutine vt_flux_volumes(this, env, vel, vof_n, vof, flux_vol, int_normal, fluids, &
        void, dt)
      import :: t2d_volume_tracker, simulation_environment, r8
      class(t2d_volume_tracker), intent(inout) :: this
      type(simulation_environment), intent(inout) :: env
      real(r8), intent(in) :: vel(:), vof_n(:,:), dt
      real(r8), intent(out) :: flux_vol(:,:), vof(:,:), int_normal(:,:,:)
      integer, intent(in) :: fluids, void
    end subroutine vt_flux_volumes

    subroutine vt_set_inflow_material(this, mat, faces)
      import :: t2d_volume_tracker
      class(t2d_volume_tracker), intent(inout) :: this
      integer, intent(in) :: mat, faces(:)
    end subroutine
  end interface

end module t2d_volume_tracker_class
