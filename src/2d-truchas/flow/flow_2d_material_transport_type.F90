!!
!! FLOW_2D_MATERIAL_TRANSPORT_TYPE
!!
!! This module defines FLOW_2D_MATERIAL_TRANSPORT, the flow-side component
!! that supplies material-resolved cell-face flux volumes. The initial
!! implementation is a one-liquid bridge: every transported volume belongs
!! to the sole liquid material and is obtained directly from face-normal
!! velocity. A later implementation will use a volume tracker to advance the
!! mobile portions of simulation-wide material composition.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module flow_2d_material_transport_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_2d_mesh_type
  implicit none
  private

  type, public :: flow_2d_material_transport
    private
    type(unstr_2d_mesh), pointer :: mesh => null() ! unowned reference
    !! Material-by-cell-face volumes transported over the current step.
    real(r8), allocatable, public :: flux_volumes(:,:)
  contains
    procedure :: init
    procedure :: advance
  end type

contains

  subroutine init(this, mesh, num_material)
    class(flow_2d_material_transport), intent(out) :: this
    type(unstr_2d_mesh), target, intent(in) :: mesh
    integer, intent(in) :: num_material

    ASSERT(num_material == 1)
    this%mesh => mesh
    allocate(this%flux_volumes(num_material, size(mesh%cface)))
  end subroutine


  !! Advance material transport from T_N to T_NP1. The one-liquid bridge does
  !! not alter material composition; it constructs the oriented liquid flux
  !! volumes that a volume tracker will later provide.
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

    this%flux_volumes = 0.0_r8
    do c = 1, this%mesh%ncell_onP
      do i = this%mesh%cstart(c), this%mesh%cstart(c+1)-1
        f = this%mesh%cface(i)
        this%flux_volumes(1,i) = dt*this%mesh%area(f)*velocity_fn(f)
        if (btest(this%mesh%cfpar(c), i-this%mesh%cstart(c)+1)) &
            this%flux_volumes(1,i) = -this%flux_volumes(1,i)
      end do
    end do
  end subroutine

end module flow_2d_material_transport_type
