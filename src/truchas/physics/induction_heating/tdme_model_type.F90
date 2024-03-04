!!
!! TDME_MODEL_TYPE
!!
!! This module provides a derived type that encapsulates data and methods for
!! a discrete form of time-domain Maxwell equations. The discretization uses
!! the mimetic Whitney finite element complex on a 3D simplicial mesh and the
!! trapezoid method in time. Derived type components include discrete boundary
!! condition data and equation coefficients. It also provides methods that
!! support time stepping, however time stepping itself is left to application
!! code.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>
!! Refactored February 2024
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module tdme_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use simpl_mesh_type
  use bndry_func1_class
  use mimetic_discretization
  implicit none
  private

  type, public :: tdme_model
    type(simpl_mesh), pointer :: mesh => null() ! spatial discretization
    real(r8) :: dt
    real(r8), allocatable :: eps(:), mu(:), sigma(:)
    class(bndry_func1), allocatable :: ebc  ! tangential E condition (nxE)
    class(bndry_func1), allocatable :: hbc  ! tangential H condition (nxH)
    real(r8), allocatable :: mtr1(:,:)
    real(r8), allocatable :: mtr2(:,:)
    real(r8), allocatable :: mtr3(:,:,:)
    real(r8), allocatable :: g0(:)  ! persistent local workspace
  contains
    procedure :: init
    procedure :: trap_res
    procedure :: advance_bfield
    procedure :: compute_joule_heat
  end type

contains

  subroutine init(this, mesh, eps, mu, sigma, dt, ebc, hbc)

    use upper_packed_matrix_procs, only: sym_matmul, upm_cong_prod

    class(tdme_model), intent(out), target :: this
    type(simpl_mesh), intent(in), target :: mesh
    real(r8), intent(in) :: eps(:), mu(:), sigma(:)
    real(r8), intent(in) :: dt
    class(bndry_func1), allocatable, intent(inout) :: ebc, hbc

    !! Local curl operator matrix
    real(r8), parameter :: curl(4,6) = reshape([0,  0,  1,  1, &
                                                0,  1,  0, -1, &
                                                0, -1, -1,  0, &
                                                1,  0,  0,  1, &
                                               -1,  0,  1,  0, &
                                                1,  1,  0,  0], shape=shape(curl))

    integer :: j
    real(r8) :: ctm2c(21), m1(21), m2(10)

    ASSERT(size(eps) == mesh%ncell)
    ASSERT(size(mu) == mesh%ncell)
    ASSERT(size(sigma) == mesh%ncell)

    this%mesh => mesh
    this%dt = dt

    call move_alloc(ebc, this%ebc)
    call move_alloc(hbc, this%hbc)

    allocate(this%eps(mesh%ncell), this%mu(mesh%ncell), this%sigma(mesh%ncell))
    this%eps = eps
    this%mu = mu
    this%sigma = sigma

    ! NB: The matrices of the discrete equations are kept in an unassembled
    ! factored form with the assembly (the outer factors) done in-effect by
    ! the residual calculation. Here we compute and store only the unassembled
    ! finite element matrices (the middle factors) of the discrete equations.

    !TODO: Is this more efficient than assembling the matrices? We end up
    ! having to do this anyway to mtr1 for the preconditioner.

    allocate(this%mtr1(21,mesh%ncell), this%mtr2(21,mesh%ncell), this%mtr3(6,4,mesh%ncell))

    do j = 1, mesh%ncell
      m1 = W1_matrix_WE(mesh, j)
      m2 = W2_matrix_WE(mesh, j)
      ctm2c = ((0.5_r8*dt)**2/mu(j)) * upm_cong_prod(4, 6, m2, curl)
      this%mtr1(:,j) = (eps(j) + 0.5_r8*dt*sigma(j)) * m1 + ctm2c
      this%mtr2(:,j) = (eps(j) - 0.5_r8*dt*sigma(j)) * m1 - ctm2c
      this%mtr3(:,:,j) = (dt/mu(j)) * transpose(sym_matmul(4, 6, m2, curl))
    end do

  end subroutine init

  !! Compute the residual of the trapezoid rule discretization of Ampere's
  !! equation given the values E0, B0 of the fields at the start of the step
  !! and the value E of the E-field at the end of the step. The discretization
  !! of Faraday's equation was used to eliminate the B-field at the end of the
  !! step in favor of these other quantities.

  subroutine trap_res(this, t0, e0, b0, e, r)

    use upper_packed_matrix_procs, only: sym_matmul

    class(tdme_model), intent(inout) :: this
    real(r8), intent(in) :: t0, e0(:), b0(:)
    real(r8), intent(inout) :: e(:) ! out due to imposing E BC
    real(r8), intent(out) :: r(:)

    integer :: j
    real(r8) :: t
    real(r8), allocatable :: e_orig(:)

    ASSERT(size(e0) == this%mesh%nedge)
    ASSERT(size(b0) == this%mesh%nface)
    ASSERT(size(e) == this%mesh%nedge)
    ASSERT(size(r) == this%mesh%nedge)

    t = t0 + this%dt

    !! Overwrite the E field with the nxE boundary data,
    !! saving the original values to restore them later.
    if (allocated(this%ebc)) then
      call this%ebc%compute(t)
      associate (index => this%ebc%index, value => this%ebc%value)
        allocate(e_orig(size(index)))
        do j = 1, size(index)
          e_orig(j) = e(index(j))
          e(index(j)) = value(j)
        end do
      end associate
      call this%mesh%edge_imap%gather_offp(e) !NB: ebc only has data for on-process edges
      !TODO: Should bndry_edge_group_builder include off-process edges too? Is it feasible?
      !      We can eliminate this parallel communication (and one below) if so.
    end if

    r = 0.0_r8
    do j = 1, this%mesh%ncell
      associate (cedge => this%mesh%cedge(:,j), cface => this%mesh%cface(:,j))
        r(cedge) = r(cedge) + sym_matmul(6, this%mtr2(:,j), e0(cedge)) &
                            - sym_matmul(6, this%mtr1(:,j), e(cedge)) &
                            + matmul(this%mtr3(:,:,j), b0(cface))
      end associate
    end do

    if (allocated(this%hbc)) then
      call this%hbc%compute(t0)
      this%g0 = this%hbc%value
      call this%hbc%compute(t)
      associate (index => this%hbc%index, g0 => this%g0, g => this%hbc%value)
        do j = 1, size(index)
          r(index(j)) = r(index(j)) + (0.5_r8*this%dt)*(g0(j) + g(j))
        end do
      end associate
    end if
    !TODO: can we drop the omit_edge_list from bndry_edge_group_builder, as the
    !      following projection achieves the same end, at the cost of doing
    !      some wasted computations on edges that are overwritten.

    !! Overwrite the residuals with the nxE BC residuals and restore the
    !! corresponding E elements to their original value.
    if (allocated(this%ebc)) then
      associate (index => this%ebc%index, value => this%ebc%value)
        do j = 1, size(index)
          r(index(j)) = value(j) - e_orig(j)
          e(index(j)) = e_orig(j)
        end do
      end associate
      call this%mesh%edge_imap%gather_offp(e) !NB: see above
    end if

    call this%mesh%edge_imap%gather_offp(r)

  end subroutine trap_res

  !! Return the B-field value B1 at the end of a step given the values
  !! E0, B0 at the start of the step and the E-field value E1 at the
  !! end of the step, according to the trapezoid time discretization
  !! of Faraday's equation.

  subroutine advance_bfield(this, e0, b0, e1, b1)
    class(tdme_model), intent(in) :: this
    real(r8), intent(in)  :: e0(:), b0(:), e1(:)
    real(r8), intent(out) :: b1(:)
    ASSERT(size(e0) == this%mesh%nedge)
    ASSERT(size(e1) == size(e0))
    ASSERT(size(b0) == this%mesh%nface)
    ASSERT(size(b1) == size(b0))
    b1 = b0 - (0.5_r8 * this%dt) * curl(this%mesh, e0 + e1)
  end subroutine

  !! Compute the cell-based Joule heats Q (a power density) that correspond
  !! to the given values E of the edge-based electric field. The size of Q
  !! may be either the number of on-process cells or all cells. However E
  !! must be given on all edges with synced values on the off-process edges
  !! to ensure correct results.

  subroutine compute_joule_heat(this, e, q)

    use upper_packed_matrix_procs, only: upm_quad_form

    class(tdme_model), intent(in) :: this
    real(r8), intent(in)  :: e(:)
    real(r8), intent(out) :: q(:)

    integer :: j
    real(r8) :: s

    ASSERT(size(e) == this%mesh%nedge)
    ASSERT(size(q) <= this%mesh%ncell)

    do j = 1, size(q)
      if (this%sigma(j) /= 0.0_r8) then
        s = upm_quad_form(W1_matrix_WE(this%mesh,j), e(this%mesh%cedge(:,j)))
        q(j) = this%sigma(j) * s / abs(this%mesh%volume(j)) ! want a source density
      else
        q(j) = 0.0_r8
      end if
    end do

  end subroutine compute_joule_heat

end module tdme_model_type
