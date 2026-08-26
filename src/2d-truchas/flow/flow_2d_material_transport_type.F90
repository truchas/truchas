!!
!! FLOW_2D_MATERIAL_TRANSPORT_TYPE
!!
!! This module defines FLOW_2D_MATERIAL_TRANSPORT, the flow-side entry point
!! for material-resolved cell-face flux volumes. It adapts the flow solver's
!! face-normal velocities to the cface-oriented interface of the two-
!! dimensional volume trackers. The initial implementation supports one
!! completely filled liquid material and uses the selected volume tracker to
!! compute its fluxes.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module flow_2d_material_transport_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use simulation_environment_type
  use unstr_2d_mesh_type
  use volume_tracker_2d_class
  implicit none
  private

  type, public :: flow_2d_material_transport
    private
    type(unstr_2d_mesh), pointer :: mesh => null() ! unowned reference
    class(volume_tracker_2d), allocatable :: tracker
    real(r8), allocatable :: vfrac_in(:,:), vfrac_out(:,:)
    real(r8), allocatable :: cface_velocity(:), interface_normal(:,:,:)
    !! Material-by-cell-face volumes transported over the current step.
    real(r8), allocatable, public :: flux_volumes(:,:)
  contains
    procedure :: init
    procedure :: advance
  end type

contains

  subroutine init(this, env, mesh, num_material, algorithm)

    use simple_volume_tracker_type
    use geometric_volume_tracker_type

    class(flow_2d_material_transport), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(unstr_2d_mesh), target, intent(in) :: mesh
    integer, intent(in) :: num_material
    character(*), intent(in), optional :: algorithm

    character(:), allocatable :: tracker_algorithm

    ASSERT(num_material == 1)
    tracker_algorithm = 'simple'
    if (present(algorithm)) tracker_algorithm = trim(algorithm)
    select case (tracker_algorithm)
    case ('simple')
      allocate(simple_volume_tracker :: this%tracker)
    case ('geometric')
      allocate(geometric_volume_tracker :: this%tracker)
    case default
      ASSERT(.false.)
    end select

    this%mesh => mesh
    allocate(this%vfrac_in(num_material,mesh%ncell), this%vfrac_out(num_material,mesh%ncell), &
        this%flux_volumes(num_material,size(mesh%cface)), this%cface_velocity(size(mesh%cface)), &
        this%interface_normal(2,num_material,mesh%ncell))
    call this%tracker%init(env, mesh, num_material, num_material, num_material, .false.)
  end subroutine


  !! Advance material transport from T_N to T_NP1. The current flow model is
  !! a single completely filled liquid, so the returned volume fraction is
  !! not retained; the tracker is nevertheless used to construct the fluxes.
  subroutine advance(this, t_n, t_np1, velocity_fn)
    class(flow_2d_material_transport), intent(inout) :: this
    real(r8), intent(in) :: t_n, t_np1
    real(r8), intent(in) :: velocity_fn(:)

    integer :: c, i, f
    real(r8) :: dt

    ASSERT(associated(this%mesh))
    ASSERT(size(this%flux_volumes,1) == 1)
    ASSERT(size(velocity_fn) >= this%mesh%nface)
    dt = t_np1 - t_n
    ASSERT(dt > 0.0_r8)

    this%vfrac_in = 0.0_r8
    this%vfrac_in(1,:) = 1.0_r8
    do c = 1, this%mesh%ncell_onP
      do i = this%mesh%cstart(c), this%mesh%cstart(c+1)-1
        f = this%mesh%cface(i)
        this%cface_velocity(i) = velocity_fn(f)
        if (btest(this%mesh%cfpar(c), i-this%mesh%cstart(c)+1)) &
            this%cface_velocity(i) = -this%cface_velocity(i)
      end do
    end do
    call this%tracker%flux_volumes(this%cface_velocity, this%vfrac_in, this%vfrac_out, &
        this%flux_volumes, this%interface_normal, 1, 0, dt)
  end subroutine

end module flow_2d_material_transport_type
