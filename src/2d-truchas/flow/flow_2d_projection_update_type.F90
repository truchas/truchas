!!
!! FLOW_2D_PROJECTION_UPDATE_TYPE
!!
!! This module defines FLOW_2D_PROJECTION_UPDATE, the collocated pressure
!! correction in a two-dimensional incompressible-flow step. It constructs a
!! pressure-consistent predicted face velocity, solves for the pressure
!! correction, then corrects face velocity, cell velocity, and pressure.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module flow_2d_projection_update_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_2d_mesh_type
  use flow_2d_operators_type
  use flow_2d_projection_type
  use flow_2d_projection_solver_type
  use flow_2d_bc_type
  use flow_2d_state_type
  use flow_domain_types
  use parallel_communication, only: global_maxval
  implicit none
  private

  type, public :: flow_2d_projection_update
    private
    type(unstr_2d_mesh), pointer :: mesh => null()  ! unowned reference
    type(flow_2d_operators), pointer :: operators => null()  ! unowned reference
    type(flow_2d_projection), pointer :: projection => null()  ! unowned reference
    type(flow_2d_projection_solver), pointer :: solver => null()  ! unowned reference
    real(r8) :: body_acceleration(2) = 0.0_r8
    real(r8), allocatable :: grad_p_old(:,:), grad_p_new(:,:), velocity_work(:,:)
    real(r8), allocatable :: gravity_head(:,:)
    real(r8), allocatable :: derivative_f(:), delta_p(:), rhs(:), flux(:)
  contains
    procedure :: init
    procedure :: correct
    procedure :: project_velocity
  end type

contains

  subroutine init(this, mesh, operators, projection, solver, body_acceleration)
    class(flow_2d_projection_update), intent(out) :: this
    type(unstr_2d_mesh), target, intent(in) :: mesh
    type(flow_2d_operators), target, intent(in) :: operators
    type(flow_2d_projection), target, intent(in) :: projection
    type(flow_2d_projection_solver), target, intent(in) :: solver
    real(r8), optional, intent(in) :: body_acceleration(:)
    this%mesh => mesh
    this%operators => operators
    this%projection => projection
    this%solver => solver
    if (present(body_acceleration)) then
      ASSERT(size(body_acceleration) == 2)
      this%body_acceleration = body_acceleration
    end if
    allocate(this%grad_p_old(2,mesh%ncell), this%grad_p_new(2,mesh%ncell), &
        this%velocity_work(2,mesh%ncell), this%derivative_f(mesh%nface), &
        this%gravity_head(2,mesh%nface), &
        this%delta_p(mesh%ncell), this%rhs(mesh%ncell_onP), this%flux(mesh%ncell_onP))
  end subroutine


  !! Project STATE%VEL_CC onto the velocity boundary conditions and the
  !! discrete face-flux continuity constraint.  The pressure correction is
  !! used only as a projection multiplier; this procedure does not alter
  !! STATE%P_CC.
  subroutine project_velocity(this, dt, inv_density_c, inv_density_f, cell_t, face_t, bc, state, stat, solved)
    class(flow_2d_projection_update), intent(inout) :: this
    real(r8), intent(in) :: dt, inv_density_c(:), inv_density_f(:)
    integer, intent(in) :: cell_t(:), face_t(:)
    type(flow_2d_bc), intent(in) :: bc
    type(flow_2d_state), intent(inout) :: state
    integer, intent(out) :: stat
    logical, optional, intent(out) :: solved

    integer :: c, f, pin_face

    ASSERT(dt > 0.0_r8)
    ASSERT(size(inv_density_c) >= this%mesh%ncell)
    ASSERT(size(inv_density_f) >= this%mesh%nface)
    ASSERT(size(cell_t) >= this%mesh%ncell)
    ASSERT(size(face_t) >= this%mesh%nface)
    stat = 0
    if (present(solved)) solved = .false.
    call this%mesh%cell_imap%gather_offp(state%vel_cc)
    do c = 1, this%mesh%ncell_onP
      if (cell_t(c) > regular_t) state%vel_cc(:,c) = 0.0_r8
    end do
    call interpolate_velocity(this, state%vel_cc, cell_t, face_t, bc, state%vel_fn)
    call apply_velocity_boundary_conditions(this%mesh, bc, face_t, state%vel_fn)
    call this%mesh%face_imap%gather_offp(state%vel_fn)

    call this%projection%assemble(inv_density_f, cell_t, face_t, bc, this%rhs, &
        bc%pressure_correction_dirichlet%value)
    call this%solver%setup()
    call this%operators%divergence(state%vel_fn, this%flux)
    this%rhs = this%rhs - this%flux/dt
    do c = 1, this%mesh%ncell_onP
      if (cell_t(c) > regular_t) this%rhs(c) = 0.0_r8
    end do
    this%delta_p = 0.0_r8
    if (global_maxval(abs(this%rhs)) > 0.0_r8) then
      if (present(solved)) solved = .true.
      call this%solver%solve(this%rhs, this%delta_p(1:this%mesh%ncell_onP), stat)
      if (stat /= 0) return
    end if
    call this%mesh%cell_imap%gather_offp(this%delta_p)

    call pressure_derivative(this, this%delta_p, cell_t, face_t, bc, this%derivative_f, correction=.true.)
    pin_face = bc%pressure_pin_face(face_t)
    if (pin_face > 0) then
      c = this%mesh%fcell(1,pin_face)
      this%derivative_f(pin_face) = -this%delta_p(c)/this%operators%normal_distance(pin_face)
    end if
    do f = 1, this%mesh%nface_onP
      state%vel_fn(f) = state%vel_fn(f) - dt*inv_density_f(f)*this%derivative_f(f)
    end do
    call apply_velocity_boundary_conditions(this%mesh, bc, face_t, state%vel_fn)
    call this%mesh%face_imap%gather_offp(state%vel_fn)

    call gradient_correction(this, this%delta_p, cell_t, face_t, this%grad_p_new, bc)
    do c = 1, this%mesh%ncell_onP
      if (cell_t(c) > regular_t) then
        state%vel_cc(:,c) = 0.0_r8
      else
        state%vel_cc(:,c) = state%vel_cc(:,c) - dt*inv_density_c(c)*this%grad_p_new(:,c)
      end if
    end do
    call this%mesh%cell_imap%gather_offp(state%vel_cc)
  end subroutine


  !! Apply one incremental pressure correction. STATE%VEL_CC is the momentum
  !! predictor velocity on entry and the corrected velocity on return.
  subroutine correct(this, dt, inv_density_c, inv_density_f, density_delta_c, cell_t, face_t, bc, state, stat, initial, solved)
    class(flow_2d_projection_update), intent(inout) :: this
    real(r8), intent(in) :: dt, inv_density_c(:), inv_density_f(:), density_delta_c(:)
    integer, intent(in) :: cell_t(:), face_t(:)
    type(flow_2d_bc), intent(in) :: bc
    type(flow_2d_state), intent(inout) :: state
    integer, intent(out) :: stat
    logical, optional, intent(in) :: initial
    logical, optional, intent(out) :: solved

    integer :: c, f, pin_face
    logical :: initial_

    ASSERT(dt > 0.0_r8)
    stat = 0
    if (present(solved)) solved = .false.
    initial_ = .false.
    if (present(initial)) initial_ = initial
    ASSERT(size(inv_density_c) >= this%mesh%ncell)
    ASSERT(size(inv_density_f) >= this%mesh%nface)
    ASSERT(size(density_delta_c) >= this%mesh%ncell)
    ASSERT(size(cell_t) >= this%mesh%ncell)
    ASSERT(size(face_t) >= this%mesh%nface)
    do c = 1, this%mesh%ncell_onP
      if (cell_t(c) > regular_t) then
        state%vel_cc(:,c) = 0.0_r8
        state%p_cc(c) = 0.0_r8
      end if
    end do
    call this%mesh%cell_imap%gather_offp(state%p_cc)
    if (initial_) then
      !! Refresh GRAVITY_HEAD, but do not use the resulting gradient as the
      !! previous committed pressure gradient for the initial predictor.
      call gradient_pressure(this, state%p_cc, cell_t, face_t, this%grad_p_new, bc, inv_density_c, density_delta_c)
      this%grad_p_old = 0.0_r8
    else
      call gradient_pressure(this, state%p_cc, cell_t, face_t, this%grad_p_old, bc, inv_density_c, density_delta_c)
    end if
    this%velocity_work = state%vel_cc
    do c = 1, this%mesh%ncell_onP
      this%velocity_work(:,c) = this%velocity_work(:,c) + dt*inv_density_c(c)*this%grad_p_old(:,c)
    end do
    call this%mesh%cell_imap%gather_offp(this%velocity_work)
    do c = 1, this%mesh%ncell_onP
      if (cell_t(c) > regular_t) this%velocity_work(:,c) = 0.0_r8
    end do
    call interpolate_velocity(this, this%velocity_work, cell_t, face_t, bc, state%vel_fn)
    call pressure_derivative(this, state%p_cc, cell_t, face_t, bc, this%derivative_f)
    do f = 1, this%mesh%nface_onP
      state%vel_fn(f) = state%vel_fn(f) - dt*inv_density_f(f)*this%derivative_f(f)
    end do
    call apply_velocity_boundary_conditions(this%mesh, bc, face_t, state%vel_fn)
    call this%mesh%face_imap%gather_offp(state%vel_fn)

    call this%projection%assemble(inv_density_f, cell_t, face_t, bc, this%rhs, &
        bc%pressure_correction_dirichlet%value)
    call this%solver%setup()
    call this%operators%divergence(state%vel_fn, this%flux)
    this%rhs = this%rhs - this%flux/dt
    do c = 1, this%mesh%ncell_onP
      if (cell_t(c) > regular_t) this%rhs(c) = 0.0_r8
    end do
    this%delta_p = 0.0_r8
    if (global_maxval(abs(this%rhs)) > 0.0_r8) then
      if (present(solved)) solved = .true.
      call this%solver%solve(this%rhs, this%delta_p(1:this%mesh%ncell_onP), stat)
      if (stat /= 0) return
    end if
    call this%mesh%cell_imap%gather_offp(this%delta_p)

    call pressure_derivative(this, this%delta_p, cell_t, face_t, bc, this%derivative_f, correction=.true.)
    pin_face = bc%pressure_pin_face(face_t)
    if (pin_face > 0) then
      c = this%mesh%fcell(1,pin_face)
      this%derivative_f(pin_face) = -this%delta_p(c)/this%operators%normal_distance(pin_face)
    end if
    do f = 1, this%mesh%nface_onP
      state%vel_fn(f) = state%vel_fn(f) - dt*inv_density_f(f)*this%derivative_f(f)
    end do
    call apply_velocity_boundary_conditions(this%mesh, bc, face_t, state%vel_fn)
    call this%mesh%face_imap%gather_offp(state%vel_fn)

    state%p_cc = state%p_cc + this%delta_p
    do c = 1, this%mesh%ncell_onP
      if (cell_t(c) > regular_t) then
        state%p_cc(c) = 0.0_r8
        state%vel_cc(:,c) = 0.0_r8
      end if
    end do
    call this%mesh%cell_imap%gather_offp(state%p_cc)
    call gradient_pressure(this, state%p_cc, cell_t, face_t, this%grad_p_new, bc, inv_density_c, density_delta_c)
    do c = 1, this%mesh%ncell_onP
      if (cell_t(c) <= regular_t) then
        state%vel_cc(:,c) = state%vel_cc(:,c) - dt*inv_density_c(c)* &
            (this%grad_p_new(:,c) - this%grad_p_old(:,c))
      end if
    end do
    call this%mesh%cell_imap%gather_offp(state%vel_cc)
  end subroutine


  subroutine gradient_pressure(this, pressure, cell_t, face_t, gradient, bc, inv_density_c, density_delta_c)
    class(flow_2d_projection_update), intent(inout) :: this
    real(r8), intent(in) :: pressure(:)
    integer, intent(in) :: cell_t(:), face_t(:)
    real(r8), intent(out) :: gradient(:,:)
    type(flow_2d_bc), intent(in) :: bc
    real(r8), intent(in) :: inv_density_c(:), density_delta_c(:)

    integer :: f, c1, c2
    real(r8) :: rho

    this%gravity_head = 0.0_r8
    do f = 1, this%mesh%nface_onP
      c1 = this%mesh%fcell(1,f)
      c2 = this%mesh%fcell(2,f)
      if (cell_t(c1) > regular_t) cycle
      rho = 1.0_r8/inv_density_c(c1) + density_delta_c(c1)
      this%gravity_head(1,f) = -this%body_acceleration(1)* &
          (this%mesh%cell_centroid(1,c1)-this%mesh%face_centroid(1,f))*rho - &
          this%body_acceleration(2)* &
          (this%mesh%cell_centroid(2,c1)-this%mesh%face_centroid(2,f))*rho
      if (c2 > 0) then
        if (cell_t(c2) > regular_t) cycle
        rho = 1.0_r8/inv_density_c(c2) + density_delta_c(c2)
        this%gravity_head(2,f) = -this%body_acceleration(1)* &
            (this%mesh%cell_centroid(1,c2)-this%mesh%face_centroid(1,f))*rho - &
            this%body_acceleration(2)* &
            (this%mesh%cell_centroid(2,c2)-this%mesh%face_centroid(2,f))*rho
      end if
    end do
    call this%mesh%face_imap%gather_offp(this%gravity_head)
    if (allocated(bc%pressure_neumann)) then
      if (allocated(bc%pressure_dirichlet)) then
        call this%operators%gradient_cc(pressure, gradient, bc%pressure_neumann, bc%pressure_dirichlet, &
            this%gravity_head, cell_t, face_t)
      else
        call this%operators%gradient_cc(pressure, gradient, normal_flux_bc=bc%pressure_neumann, &
            gravity_head=this%gravity_head, cell_t=cell_t, face_t=face_t)
      end if
    else if (allocated(bc%pressure_dirichlet)) then
      call this%operators%gradient_cc(pressure, gradient, dirichlet_bc=bc%pressure_dirichlet, &
          gravity_head=this%gravity_head, cell_t=cell_t, face_t=face_t)
    else
      call this%operators%gradient_cc(pressure, gradient, gravity_head=this%gravity_head, &
          cell_t=cell_t, face_t=face_t)
    end if
  end subroutine


  subroutine gradient_correction(this, pressure, cell_t, face_t, gradient, bc)
    class(flow_2d_projection_update), intent(in) :: this
    real(r8), intent(in) :: pressure(:)
    integer, intent(in) :: cell_t(:), face_t(:)
    real(r8), intent(out) :: gradient(:,:)
    type(flow_2d_bc), intent(in) :: bc

    call this%operators%gradient_cc(pressure, gradient, bc%pressure_neumann, &
        bc%pressure_correction_dirichlet, cell_t=cell_t, face_t=face_t)
  end subroutine


  subroutine pressure_derivative(this, pressure, cell_t, face_t, bc, derivative, correction)
    class(flow_2d_projection_update), intent(inout) :: this
    real(r8), intent(in) :: pressure(:)
    integer, intent(in) :: cell_t(:), face_t(:)
    type(flow_2d_bc), intent(in) :: bc
    real(r8), intent(out) :: derivative(:)
    logical, optional, intent(in) :: correction

    logical :: correction_

    correction_ = .false.
    if (present(correction)) correction_ = correction
    if (correction_) then
      if (allocated(bc%pressure_neumann)) then
        if (allocated(bc%pressure_dirichlet)) then
          call this%operators%derivative_cf(pressure, derivative, bc%pressure_neumann, &
              bc%pressure_dirichlet, bc%pressure_correction_dirichlet%value, face_t=face_t)
        else
          call this%operators%derivative_cf(pressure, derivative, normal_flux_bc=bc%pressure_neumann, face_t=face_t)
        end if
      else if (allocated(bc%pressure_dirichlet)) then
        call this%operators%derivative_cf(pressure, derivative, &
            dirichlet_bc=bc%pressure_dirichlet, &
            dirichlet_value=bc%pressure_correction_dirichlet%value, face_t=face_t)
      else
        call this%operators%derivative_cf(pressure, derivative, face_t=face_t)
      end if
      return
    end if

    if (allocated(bc%pressure_neumann)) then
      if (allocated(bc%pressure_dirichlet)) then
        call this%operators%derivative_cf(pressure, derivative, bc%pressure_neumann, bc%pressure_dirichlet, &
            gravity_head=this%gravity_head, face_t=face_t)
      else
        call this%operators%derivative_cf(pressure, derivative, normal_flux_bc=bc%pressure_neumann, &
            gravity_head=this%gravity_head, face_t=face_t)
      end if
    else if (allocated(bc%pressure_dirichlet)) then
      call this%operators%derivative_cf(pressure, derivative, dirichlet_bc=bc%pressure_dirichlet, &
          gravity_head=this%gravity_head, face_t=face_t)
    else
      call this%operators%derivative_cf(pressure, derivative, gravity_head=this%gravity_head, face_t=face_t)
    end if
  end subroutine


  subroutine interpolate_velocity(this, velocity, cell_t, face_t, bc, velocity_f)
    class(flow_2d_projection_update), intent(in) :: this
    real(r8), intent(in) :: velocity(:,:)
    integer, intent(in) :: cell_t(:), face_t(:)
    type(flow_2d_bc), intent(in) :: bc
    real(r8), intent(out) :: velocity_f(:)

    integer :: f, c1, c2

    if (allocated(bc%velocity_zero_normal)) then
      if (allocated(bc%velocity_dirichlet)) then
        call this%operators%interpolate_cf(velocity, velocity_f, bc%velocity_zero_normal, bc%velocity_dirichlet)
      else
        call this%operators%interpolate_cf(velocity, velocity_f, zero_normal_bc=bc%velocity_zero_normal)
      end if
    else if (allocated(bc%velocity_dirichlet)) then
      call this%operators%interpolate_cf(velocity, velocity_f, dirichlet_bc=bc%velocity_dirichlet)
    else
      call this%operators%interpolate_cf(velocity, velocity_f)
    end if
    do f = 1, this%mesh%nface_onP
      c1 = this%mesh%fcell(1,f)
      c2 = this%mesh%fcell(2,f)
      if (face_t(f) == regular_void_t) then
        if (cell_t(c1) <= regular_t) then
          velocity_f(f) = dot_product(this%mesh%normal(:,f), velocity(:,c1))/this%mesh%area(f)
        else if (c2 > 0) then
          if (cell_t(c2) <= regular_t) then
            velocity_f(f) = dot_product(this%mesh%normal(:,f), velocity(:,c2))/this%mesh%area(f)
          else
            velocity_f(f) = 0.0_r8
          end if
        else
          velocity_f(f) = 0.0_r8
        end if
      else if (face_t(f) > regular_t) then
        velocity_f(f) = 0.0_r8
      end if
    end do
  end subroutine


  subroutine apply_velocity_boundary_conditions(mesh, bc, face_t, velocity_f)
    type(unstr_2d_mesh), intent(in) :: mesh
    type(flow_2d_bc), intent(in) :: bc
    integer, intent(in) :: face_t(:)
    real(r8), intent(inout) :: velocity_f(:)

    integer :: i, f

    if (allocated(bc%velocity_zero_normal)) then
      do i = 1, size(bc%velocity_zero_normal%index)
        f = bc%velocity_zero_normal%index(i)
        if (f > mesh%nface_onP) cycle
        if (face_t(f) /= regular_t) cycle
        velocity_f(f) = 0.0_r8
      end do
    end if
    if (allocated(bc%velocity_dirichlet)) then
      do i = 1, size(bc%velocity_dirichlet%index)
        f = bc%velocity_dirichlet%index(i)
        if (f > mesh%nface_onP) cycle
        if (face_t(f) /= regular_t) cycle
        velocity_f(f) = dot_product(mesh%unit_normal(:,f), bc%velocity_dirichlet%value(:,i))
      end do
    end if
  end subroutine

end module flow_2d_projection_update_type
