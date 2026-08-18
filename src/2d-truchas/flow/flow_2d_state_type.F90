!!
!! FLOW_2D_STATE_TYPE
!!
!! This module defines FLOW_2D_STATE, the cell- and face-centered state of a
!! two-dimensional incompressible flow calculation. It stores velocity at
!! cell centers and as face-normal values, and dynamic pressure at cell
!! centers. Values are stored on the full local mesh; callers update owned
!! entries and use GATHER_OFFP before values are needed on ghost entities.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

module flow_2d_state_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_2d_mesh_type
  implicit none
  private

  type, public :: flow_2d_state
    private
    type(unstr_2d_mesh), pointer :: mesh => null()  ! unowned reference
    real(r8), allocatable, public :: vel_cc(:,:)  ! (2, local cell)
    real(r8), allocatable, public :: vel_fn(:)    ! local face-normal velocity
    real(r8), allocatable, public :: p_cc(:)      ! local dynamic pressure
  contains
    procedure :: init
    procedure :: set_zero
    procedure :: gather_offp
  end type

contains

  subroutine init(this, mesh)
    class(flow_2d_state), intent(out) :: this
    type(unstr_2d_mesh), target, intent(in) :: mesh

    this%mesh => mesh
    allocate(this%vel_cc(2, mesh%ncell), this%vel_fn(mesh%nface), this%p_cc(mesh%ncell))
    call this%set_zero()
  end subroutine


  subroutine set_zero(this)
    class(flow_2d_state), intent(inout) :: this

    this%vel_cc = 0.0_r8
    this%vel_fn = 0.0_r8
    this%p_cc = 0.0_r8
  end subroutine


  subroutine gather_offp(this)
    class(flow_2d_state), intent(inout) :: this

    call this%mesh%cell_imap%gather_offp(this%vel_cc)
    call this%mesh%face_imap%gather_offp(this%vel_fn)
    call this%mesh%cell_imap%gather_offp(this%p_cc)
  end subroutine

end module flow_2d_state_type
