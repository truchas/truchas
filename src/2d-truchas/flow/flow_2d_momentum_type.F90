!!
!! FLOW_2D_MOMENTUM_TYPE
!!
!! This module defines FLOW_2D_MOMENTUM, the block-structured finite-volume
!! operator for the cell-centered velocity predictor in two-dimensional
!! incompressible flow. The matrix has one 2-by-2 block per cell pair. It
!! assembles unsteady inertial and first-order viscous terms, and accumulates
!! a conservative donor-cell momentum-transport contribution to the RHS from
!! material flux volumes.
!!
!! Velocity Dirichlet and free-slip conditions are eliminated through the
!! boundary face fluxes. Free slip contributes a normal outer-product block,
!! which couples the two Cartesian velocity components.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module flow_2d_momentum_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_2d_mesh_type
  use flow_2d_operators_type
  use flow_2d_bc_type
  use pbsr_matrix_type
  implicit none
  private

  type, public :: flow_2d_momentum
    private
    type(unstr_2d_mesh), pointer :: mesh => null()  ! unowned reference
    type(flow_2d_operators), pointer :: operators => null()  ! unowned reference
    type(pbsr_matrix) :: matrix_
  contains
    procedure :: init
    procedure :: assemble
    procedure :: add_advective_rhs
    procedure :: matrix
  end type

contains

  subroutine init(this, mesh, operators)
    class(flow_2d_momentum), intent(out) :: this
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
    call this%matrix_%init(2, graph, take_graph=.true.)
  end subroutine


  !! Add -div(rho*u*u) to RHS using first-order donor-cell transport.
  !! FLUX_VOLUMES is material-by-cell-face and contains signed material
  !! volumes transported over the pending time step; positive values leave
  !! the associated cell. DENSITY gives the corresponding material densities.
  subroutine add_advective_rhs(this, density, velocity_cc, flux_volumes, bc, rhs)
    class(flow_2d_momentum), intent(in) :: this
    real(r8), intent(in) :: density(:), velocity_cc(:,:), flux_volumes(:,:)
    type(flow_2d_bc), intent(in) :: bc
    real(r8), intent(inout) :: rhs(:,:)

    integer :: c, i, f, neighbor, n
    real(r8) :: mass_flux

    ASSERT(size(density) > 0)
    ASSERT(size(velocity_cc,1) == 2 .and. size(velocity_cc,2) >= this%mesh%ncell)
    ASSERT(size(flux_volumes,1) >= size(density) .and. size(flux_volumes,2) == size(this%mesh%cface))
    ASSERT(size(rhs,1) == 2 .and. size(rhs,2) == this%mesh%ncell_onP)

    do c = 1, this%mesh%ncell_onP
      do i = this%mesh%cstart(c), this%mesh%cstart(c+1)-1
        f = this%mesh%cface(i)
        neighbor = this%mesh%cnhbr(i)
        mass_flux = dot_product(density, flux_volumes(:size(density),i))

        if (any(flux_volumes(:,i) > 0.0_r8)) then
          rhs(:,c) = rhs(:,c) - mass_flux*velocity_cc(:,c)
        else if (neighbor > 0) then
          rhs(:,c) = rhs(:,c) - mass_flux*velocity_cc(:,neighbor)
        end if
      end do
    end do

    do n = 1, size(bc%velocity_dirichlet%index)
      f = bc%velocity_dirichlet%index(n)
      if (f > this%mesh%nface_onP) cycle
      c = this%mesh%fcell(1,f)
      i = this%mesh%cstart(c) - 1 + findloc(this%mesh%cface(this%mesh%cstart(c):this%mesh%cstart(c+1)-1), f, dim=1)
      if (.not.any(flux_volumes(:,i) < 0.0_r8)) cycle
      mass_flux = dot_product(density, flux_volumes(:size(density),i))
      rhs(:,c) = rhs(:,c) - mass_flux*bc%velocity_dirichlet%value(:,n)
    end do

    do n = 1, size(bc%pressure_dirichlet%index)
      f = bc%pressure_dirichlet%index(n)
      if (f > this%mesh%nface_onP) cycle
      c = this%mesh%fcell(1,f)
      i = this%mesh%cstart(c) - 1 + findloc(this%mesh%cface(this%mesh%cstart(c):this%mesh%cstart(c+1)-1), f, dim=1)
      if (.not.any(flux_volumes(:,i) < 0.0_r8)) cycle
      mass_flux = dot_product(density, flux_volumes(:size(density),i))
      rhs(:,c) = rhs(:,c) - mass_flux*velocity_cc(:,c)
    end do
  end subroutine


  !! Assemble rho*volume*I - dt*div(mu grad(u)) and its velocity-Dirichlet
  !! contribution to RHS. BC must already have been evaluated at the required
  !! time.
  subroutine assemble(this, dt, density_c, viscosity_f, bc, rhs)
    class(flow_2d_momentum), intent(inout) :: this
    real(r8), intent(in) :: dt, density_c(:), viscosity_f(:)
    type(flow_2d_bc), intent(in) :: bc
    real(r8), intent(out) :: rhs(:,:)

    integer :: c, i, f, neighbor, n
    real(r8) :: coefficient, identity(2,2), normal(2), block(2,2)

    ASSERT(dt >= 0.0_r8)
    ASSERT(size(density_c) >= this%mesh%ncell)
    ASSERT(size(viscosity_f) >= this%mesh%nface)
    ASSERT(size(rhs,1) == 2 .and. size(rhs,2) == this%mesh%ncell_onP)

    identity = 0.0_r8
    identity(1,1) = 1.0_r8
    identity(2,2) = 1.0_r8
    call this%matrix_%set_all(0.0_r8)
    rhs = 0.0_r8

    do c = 1, this%mesh%ncell_onP
      call this%matrix_%add_to(c, c, density_c(c)*this%mesh%volume(c)*identity)
      do i = this%mesh%cstart(c), this%mesh%cstart(c+1)-1
        f = this%mesh%cface(i)
        neighbor = this%mesh%cnhbr(i)
        if (neighbor == 0) cycle
        coefficient = dt*viscosity_f(f)*this%mesh%area(f)/this%operators%normal_distance(f)
        call this%matrix_%add_to(c, c, coefficient*identity)
        call this%matrix_%add_to(c, neighbor, -coefficient*identity)
      end do
    end do

    if (allocated(bc%velocity_dirichlet)) then
      do n = 1, size(bc%velocity_dirichlet%index)
        f = bc%velocity_dirichlet%index(n)
        if (f > this%mesh%nface_onP) cycle
        c = this%mesh%fcell(1,f)
        ASSERT(this%mesh%fcell(2,f) == 0)
        coefficient = dt*viscosity_f(f)*this%mesh%area(f)/this%operators%normal_distance(f)
        call this%matrix_%add_to(c, c, coefficient*identity)
        rhs(:,c) = rhs(:,c) + coefficient*bc%velocity_dirichlet%value(:,n)
      end do
    end if

    if (allocated(bc%velocity_zero_normal)) then
      do n = 1, size(bc%velocity_zero_normal%index)
        f = bc%velocity_zero_normal%index(n)
        if (f > this%mesh%nface_onP) cycle
        c = this%mesh%fcell(1,f)
        ASSERT(this%mesh%fcell(2,f) == 0)
        coefficient = dt*viscosity_f(f)*this%mesh%area(f)/this%operators%normal_distance(f)
        normal = this%mesh%unit_normal(:,f)
        block = coefficient*reshape([normal(1)**2, normal(1)*normal(2), &
            normal(1)*normal(2), normal(2)**2], [2,2])
        call this%matrix_%add_to(c, c, block)
      end do
    end if
  end subroutine


  function matrix(this)
    class(flow_2d_momentum), intent(in), target :: this
    type(pbsr_matrix), pointer :: matrix

    matrix => this%matrix_
  end function

end module flow_2d_momentum_type
