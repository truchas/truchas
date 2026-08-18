!!
!! FLOW_2D_PROJECTION_TYPE
!!
!! This module defines FLOW_2D_PROJECTION, the pressure Poisson system used
!! by a two-dimensional incompressible-flow projection. The system represents
!! -div(rho^-1 grad(p)) with the first-order face-normal operators.
!!
!! Pressure Dirichlet data are incorporated in both the matrix and right-hand
!! side. When all pressure boundaries are homogeneous Neumann conditions, a
!! regular boundary face is assigned zero Dirichlet pressure to remove the
!! nullspace without replacing a cell equation.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module flow_2d_projection_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_2d_mesh_type
  use flow_2d_operators_type
  use flow_2d_bc_type
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


  !! Assemble -div(rho^-1 grad(p)) and its pressure-Dirichlet contribution to
  !! RHS. BC must already have been evaluated at the required time.
  subroutine assemble(this, inv_density_f, bc, rhs)
    class(flow_2d_projection), intent(inout) :: this
    real(r8), intent(in) :: inv_density_f(:)
    type(flow_2d_bc), intent(in) :: bc
    real(r8), intent(out) :: rhs(:)

    integer :: c, i, f, neighbor, pin_face, n
    real(r8) :: coefficient

    ASSERT(size(inv_density_f) >= this%mesh%nface)
    ASSERT(size(rhs) == this%mesh%ncell_onP)
    pin_face = bc%pressure_pin_face()
    ASSERT(pin_face >= 0 .and. pin_face <= this%mesh%nface_onP)
    if (pin_face > 0) then
      ASSERT(this%mesh%fcell(2,pin_face) == 0)
    end if

    call this%matrix_%set_all(0.0_r8)
    rhs = 0.0_r8
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

    if (allocated(bc%pressure_dirichlet)) then
      do n = 1, size(bc%pressure_dirichlet%index)
        f = bc%pressure_dirichlet%index(n)
        if (f > this%mesh%nface_onP) cycle
        c = this%mesh%fcell(1,f)
        ASSERT(this%mesh%fcell(2,f) == 0)
        coefficient = this%mesh%area(f)*inv_density_f(f)/this%operators%normal_distance(f)
        call this%matrix_%add_to(c, c, coefficient)
        rhs(c) = rhs(c) + coefficient*bc%pressure_dirichlet%value(n)
      end do
    end if

    if (pin_face > 0) then
      c = this%mesh%fcell(1,pin_face)
      coefficient = this%mesh%area(pin_face)*inv_density_f(pin_face)/ &
          this%operators%normal_distance(pin_face)
      call this%matrix_%add_to(c, c, coefficient)
    end if
  end subroutine


  function matrix(this)
    class(flow_2d_projection), intent(in), target :: this
    type(pcsr_matrix), pointer :: matrix

    matrix => this%matrix_
  end function

end module flow_2d_projection_type
