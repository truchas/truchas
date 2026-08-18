!!
!! FLOW_2D_PROJECTION_TYPE
!!
!! This module defines FLOW_2D_PROJECTION, the pressure Poisson system used
!! by a two-dimensional incompressible-flow projection. The system represents
!! -div(rho^-1 grad(p)) with the first-order face-normal operators.
!!
!! A pressure reference is represented by PIN_FACE: a regular boundary face
!! is assigned zero Dirichlet pressure through its normal face coefficient.
!! This removes the all-Neumann nullspace without replacing a cell equation.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module flow_2d_projection_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_2d_mesh_type
  use flow_2d_operators_type
  use pcsr_matrix_type
  implicit none
  private

  type, public :: flow_2d_projection
    private
    type(unstr_2d_mesh), pointer :: mesh => null()  ! unowned reference
    type(flow_2d_operators), pointer :: operators => null()  ! unowned reference
    type(pcsr_matrix) :: matrix_
  contains
    procedure :: init
    procedure :: assemble
    procedure :: matrix
  end type

contains

  subroutine init(this, mesh, operators)
    class(flow_2d_projection), intent(out) :: this
    type(unstr_2d_mesh), target, intent(in) :: mesh
    type(flow_2d_operators), target, intent(in) :: operators

    type(pcsr_graph), pointer :: graph
    integer :: c, i, neighbor

    this%mesh => mesh
    this%operators => operators
    allocate(graph)
    call graph%init(mesh%cell_imap)
    do c = 1, mesh%ncell_onP
      call graph%add_edge(c, c)
      do i = mesh%cstart(c), mesh%cstart(c+1)-1
        neighbor = mesh%cnhbr(i)
        if (neighbor > 0) call graph%add_edge(c, neighbor)
      end do
    end do
    call graph%add_complete()
    call this%matrix_%init(graph, take_graph=.true.)
  end subroutine


  !! Assemble -div(rho^-1 grad(p)). PIN_FACE is optional and is a local,
  !! on-process boundary face at which p=0 is imposed to remove an otherwise
  !! all-Neumann pressure nullspace.
  subroutine assemble(this, inv_density_f, pin_face)
    class(flow_2d_projection), intent(inout) :: this
    real(r8), intent(in) :: inv_density_f(:)
    integer, optional, intent(in) :: pin_face

    integer :: c, i, f, neighbor, pin_face_
    real(r8) :: coefficient

    ASSERT(size(inv_density_f) >= this%mesh%nface)
    pin_face_ = 0
    if (present(pin_face)) pin_face_ = pin_face
    ASSERT(pin_face_ >= 0 .and. pin_face_ <= this%mesh%nface_onP)
    if (pin_face_ > 0) then
      ASSERT(this%mesh%fcell(2,pin_face_) == 0)
    end if

    call this%matrix_%set_all(0.0_r8)
    do c = 1, this%mesh%ncell_onP
      do i = this%mesh%cstart(c), this%mesh%cstart(c+1)-1
        f = this%mesh%cface(i)
        neighbor = this%mesh%cnhbr(i)
        if (neighbor == 0) cycle
        coefficient = this%mesh%area(f)*inv_density_f(f)/this%operators%normal_distance(f)
        call this%matrix_%add_to(c, c, coefficient)
        call this%matrix_%add_to(c, neighbor, -coefficient)
      end do
    end do

    if (pin_face_ > 0) then
      c = this%mesh%fcell(1,pin_face_)
      coefficient = this%mesh%area(pin_face_)*inv_density_f(pin_face_)/ &
          this%operators%normal_distance(pin_face_)
      call this%matrix_%add_to(c, c, coefficient)
    end if
  end subroutine


  function matrix(this)
    class(flow_2d_projection), intent(in), target :: this
    type(pcsr_matrix), pointer :: matrix

    matrix => this%matrix_
  end function

end module flow_2d_projection_type
