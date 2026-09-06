!!
!! T2D_FLOW_MATERIAL_TRANSPORT_TYPE
!!
!! This module defines T2D_FLOW_MATERIAL_TRANSPORT, the flow-side entry point
!! for material-resolved cell-face flux volumes. It adapts the flow solver's
!! face-normal velocities to the cface-oriented interface of the two-
!! dimensional volume trackers. The caller provides the current reduced
!! material distribution and retains its authoritative state; this type holds
!! only the trial distribution produced by the tracker.
!!
!! The reduced distribution is ordered as real fluid phases, an optional VOID
!! slot, and an optional lumped SOLID slot. The layout maintains the phase and
!! parent-material mappings; the independent full reconstruction-priority
!! ordering is supplied by the caller.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module t2d_flow_material_transport_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use simulation_environment_type
  use t2d_unstr_mesh_type
  use t2d_volume_tracker_class
  implicit none
  private

  type, public :: t2d_flow_material_transport
    private
    type(t2d_unstr_mesh), pointer :: mesh => null() ! unowned reference
    class(t2d_volume_tracker), allocatable :: tracker
    integer :: nrealfluid, nfluid
    real(r8), allocatable :: vfrac_out(:,:)
    real(r8), allocatable :: cface_velocity(:), interface_normal(:,:,:)
    !! Material-by-cell-face volumes transported over the current step.
    real(r8), allocatable, public :: flux_volumes(:,:)
  contains
    procedure :: init
    procedure :: set_inflow_material
    procedure :: advance
    procedure :: get_trial_volume_fractions
  end type

contains

  subroutine init(this, env, mesh, nrealfluid, nfluid, nmat, algorithm, priority, cutoff)

    use t2d_simple_volume_tracker_type
    use t2d_geometric_volume_tracker_type

    class(t2d_flow_material_transport), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(t2d_unstr_mesh), target, intent(in) :: mesh
    integer, intent(in) :: nrealfluid, nfluid, nmat
    character(*), intent(in), optional :: algorithm
    integer, intent(in), optional :: priority(:)
    real(r8), intent(in), optional :: cutoff

    character(:), allocatable :: tracker_algorithm
    integer, allocatable :: tracker_priority(:)
    integer :: i

    ASSERT(nrealfluid >= 0)
    ASSERT(nrealfluid <= nfluid)
    ASSERT(nfluid <= nmat)
    allocate(tracker_priority(nmat))
    tracker_priority = [(i, i=1,nmat)]
    if (present(priority)) tracker_priority = priority
    ASSERT(all(tracker_priority > 0) .and. all(tracker_priority <= nmat))
    ASSERT(all([(count(tracker_priority == i), i=1,nmat)] == 1))
    tracker_algorithm = 'simple'
    if (present(algorithm)) tracker_algorithm = trim(algorithm)
    select case (tracker_algorithm)
    case ('simple')
      allocate(t2d_simple_volume_tracker :: this%tracker)
    case ('geometric')
      allocate(t2d_geometric_volume_tracker :: this%tracker)
    case default
      ASSERT(.false.)
    end select

    this%mesh => mesh
    this%nrealfluid = nrealfluid
    this%nfluid = nfluid
    allocate(this%vfrac_out(nmat,mesh%ncell), this%flux_volumes(nfluid,size(mesh%cface)), &
        this%cface_velocity(size(mesh%cface)), this%interface_normal(2,nmat,mesh%ncell))
    if (present(cutoff)) then
      call this%tracker%init(env, mesh, nrealfluid, nfluid, nmat, .false., tracker_priority, cutoff)
    else
      call this%tracker%init(env, mesh, nrealfluid, nfluid, nmat, .false., tracker_priority)
    end if
  end subroutine


  !! Set a fixed reduced-material slot for inflow through the given local
  !! boundary faces.  Slot zero retains the tracker's proportional-inflow
  !! default.
  subroutine set_inflow_material(this, material, faces)
    class(t2d_flow_material_transport), intent(inout) :: this
    integer, intent(in) :: material, faces(:)

    ASSERT(material >= 0 .and. material <= this%nfluid)
    call this%tracker%set_inflow_material(material, faces)
  end subroutine


  !! Advance material transport from T_N to T_NP1, retaining the resulting
  !! distribution as trial state until the caller accepts the coupled step.
  subroutine advance(this, env, t_n, t_np1, velocity_fn, vfrac_n)
    class(t2d_flow_material_transport), intent(inout) :: this
    type(simulation_environment), intent(inout) :: env
    real(r8), intent(in) :: t_n, t_np1
    real(r8), intent(in) :: velocity_fn(:)
    real(r8), intent(in) :: vfrac_n(:,:)

    integer :: c, i, f
    real(r8) :: dt

    ASSERT(associated(this%mesh))
    ASSERT(size(velocity_fn) >= this%mesh%nface)
    ASSERT(size(vfrac_n,1) == size(this%vfrac_out,1))
    ASSERT(size(vfrac_n,2) == size(this%vfrac_out,2))
    dt = t_np1 - t_n
    ASSERT(dt > 0.0_r8)

    do c = 1, this%mesh%ncell_onP
      do i = this%mesh%cstart(c), this%mesh%cstart(c+1)-1
        f = this%mesh%cface(i)
        this%cface_velocity(i) = velocity_fn(f)
        if (btest(this%mesh%cfpar(c), i-this%mesh%cstart(c)+1)) &
            this%cface_velocity(i) = -this%cface_velocity(i)
      end do
    end do
    call this%tracker%flux_volumes(env, this%cface_velocity, vfrac_n, this%vfrac_out, &
        this%flux_volumes, this%interface_normal, this%nrealfluid, this%nfluid-this%nrealfluid, dt)
  end subroutine


  !! Return a no-copy view of the distribution produced by the most recent
  !! advance call. It is trial state; the caller decides whether to adopt it.
  subroutine get_trial_volume_fractions(this, vfrac)
    class(t2d_flow_material_transport), target, intent(in) :: this
    real(r8), pointer, intent(out) :: vfrac(:,:)

    vfrac => this%vfrac_out
  end subroutine

end module t2d_flow_material_transport_type
