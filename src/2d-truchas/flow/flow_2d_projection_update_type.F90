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
    real(r8) :: thermal_expan_coef = 0.0_r8, expan_ref_temp = 0.0_r8
    real(r8), allocatable :: grad_p_old(:,:), grad_p_new(:,:), velocity_work(:,:)
    real(r8), allocatable :: gravity_head(:,:)
    real(r8), allocatable :: buoyancy_temperature(:)
    real(r8), allocatable :: derivative_f(:), delta_p(:), rhs(:), flux(:)
  contains
    procedure :: init
    procedure :: set_buoyancy_temperature
    procedure :: correct
    procedure :: project_velocity
  end type

contains

  subroutine init(this, mesh, operators, projection, solver, body_acceleration, &
      thermal_expan_coef, expan_ref_temp)
    class(flow_2d_projection_update), intent(out) :: this
    type(unstr_2d_mesh), target, intent(in) :: mesh
    type(flow_2d_operators), target, intent(in) :: operators
    type(flow_2d_projection), target, intent(in) :: projection
    type(flow_2d_projection_solver), target, intent(in) :: solver
    real(r8), optional, intent(in) :: body_acceleration(:)
    real(r8), optional, intent(in) :: thermal_expan_coef, expan_ref_temp

    this%mesh => mesh
    this%operators => operators
    this%projection => projection
    this%solver => solver
    if (present(body_acceleration)) then
      ASSERT(size(body_acceleration) == 2)
      this%body_acceleration = body_acceleration
    end if
    if (present(thermal_expan_coef)) then
      ASSERT(thermal_expan_coef >= 0.0_r8)
      this%thermal_expan_coef = thermal_expan_coef
    end if
    if (present(expan_ref_temp)) this%expan_ref_temp = expan_ref_temp
    allocate(this%grad_p_old(2,mesh%ncell), this%grad_p_new(2,mesh%ncell), &
        this%velocity_work(2,mesh%ncell), this%derivative_f(mesh%nface), &
        this%gravity_head(2,mesh%nface), &
        this%delta_p(mesh%ncell), this%rhs(mesh%ncell_onP), this%flux(mesh%ncell_onP))
    if (this%thermal_expan_coef > 0.0_r8) then
      allocate(this%buoyancy_temperature(mesh%ncell), source=0.0_r8)
    end if
  end subroutine


  subroutine set_buoyancy_temperature(this, temperature)
    class(flow_2d_projection_update), intent(inout) :: this
    real(r8), intent(in) :: temperature(:)

    if (.not.allocated(this%buoyancy_temperature)) return
    ASSERT(size(temperature) == this%mesh%ncell_onP)
    this%buoyancy_temperature(:this%mesh%ncell_onP) = temperature
    call this%mesh%cell_imap%gather_offp(this%buoyancy_temperature)
  end subroutine


  !! Project STATE%VEL_CC onto the velocity boundary conditions and the
  !! discrete face-flux continuity constraint.  The pressure correction is
  !! used only as a projection multiplier; this procedure does not alter
  !! STATE%P_CC.
  subroutine project_velocity(this, dt, inv_density_c, inv_density_f, bc, state, stat)
    class(flow_2d_projection_update), intent(inout) :: this
    real(r8), intent(in) :: dt, inv_density_c(:), inv_density_f(:)
    type(flow_2d_bc), intent(in) :: bc
    type(flow_2d_state), intent(inout) :: state
    integer, intent(out) :: stat

    integer :: c, f, pin_face

    ASSERT(dt > 0.0_r8)
    stat = 0
    call this%mesh%cell_imap%gather_offp(state%vel_cc)
    call interpolate_velocity(this, state%vel_cc, bc, state%vel_fn)
    call apply_velocity_boundary_conditions(this%mesh, bc, state%vel_fn)
    call this%mesh%face_imap%gather_offp(state%vel_fn)

    call this%projection%assemble(inv_density_f, bc, this%rhs, bc%pressure_correction_dirichlet%value)
    call this%solver%setup()
    call this%operators%divergence(state%vel_fn, this%flux)
    this%rhs = this%rhs - this%flux/dt
    this%delta_p = 0.0_r8
    if (global_maxval(abs(this%rhs)) > 0.0_r8) then
      call this%solver%solve(this%rhs, this%delta_p(1:this%mesh%ncell_onP), stat)
      if (stat /= 0) return
    end if
    call this%mesh%cell_imap%gather_offp(this%delta_p)

    call pressure_derivative(this, this%delta_p, bc, this%derivative_f, correction=.true.)
    pin_face = bc%pressure_pin_face()
    if (pin_face > 0) then
      c = this%mesh%fcell(1,pin_face)
      this%derivative_f(pin_face) = -this%delta_p(c)/this%operators%normal_distance(pin_face)
    end if
    do f = 1, this%mesh%nface_onP
      state%vel_fn(f) = state%vel_fn(f) - dt*inv_density_f(f)*this%derivative_f(f)
    end do
    call apply_velocity_boundary_conditions(this%mesh, bc, state%vel_fn)
    call this%mesh%face_imap%gather_offp(state%vel_fn)

    call gradient_correction(this, this%delta_p, this%grad_p_new, bc)
    do c = 1, this%mesh%ncell_onP
      state%vel_cc(:,c) = state%vel_cc(:,c) - dt*inv_density_c(c)*this%grad_p_new(:,c)
    end do
    call this%mesh%cell_imap%gather_offp(state%vel_cc)
  end subroutine


  !! Apply one incremental pressure correction. STATE%VEL_CC is the momentum
  !! predictor velocity on entry and the corrected velocity on return.
  subroutine correct(this, dt, inv_density_c, inv_density_f, bc, state, stat, initial)
    class(flow_2d_projection_update), intent(inout) :: this
    real(r8), intent(in) :: dt, inv_density_c(:), inv_density_f(:)
    type(flow_2d_bc), intent(in) :: bc
    type(flow_2d_state), intent(inout) :: state
    integer, intent(out) :: stat
    logical, optional, intent(in) :: initial

    integer :: c, f, pin_face
    logical :: initial_

    ASSERT(dt > 0.0_r8)
    stat = 0
    initial_ = .false.
    if (present(initial)) initial_ = initial
    ASSERT(size(inv_density_c) >= this%mesh%ncell)
    ASSERT(size(inv_density_f) >= this%mesh%nface)
    if (initial_) then
      !! Refresh GRAVITY_HEAD, but do not use the resulting gradient as the
      !! previous committed pressure gradient for the initial predictor.
      call gradient_pressure(this, state%p_cc, this%grad_p_new, bc, inv_density_c)
      this%grad_p_old = 0.0_r8
    else
      call gradient_pressure(this, state%p_cc, this%grad_p_old, bc, inv_density_c)
    end if
    this%velocity_work = state%vel_cc
    do c = 1, this%mesh%ncell_onP
      this%velocity_work(:,c) = this%velocity_work(:,c) + dt*inv_density_c(c)*this%grad_p_old(:,c)
    end do
    call this%mesh%cell_imap%gather_offp(this%velocity_work)
    call interpolate_velocity(this, this%velocity_work, bc, state%vel_fn)
    call pressure_derivative(this, state%p_cc, bc, this%derivative_f)
    do f = 1, this%mesh%nface_onP
      state%vel_fn(f) = state%vel_fn(f) - dt*inv_density_f(f)*this%derivative_f(f)
    end do
    call apply_velocity_boundary_conditions(this%mesh, bc, state%vel_fn)
    call this%mesh%face_imap%gather_offp(state%vel_fn)

    call this%projection%assemble(inv_density_f, bc, this%rhs, bc%pressure_correction_dirichlet%value)
    call this%solver%setup()
    call this%operators%divergence(state%vel_fn, this%flux)
    this%rhs = this%rhs - this%flux/dt
    this%delta_p = 0.0_r8
    if (global_maxval(abs(this%rhs)) > 0.0_r8) then
      call this%solver%solve(this%rhs, this%delta_p(1:this%mesh%ncell_onP), stat)
      if (stat /= 0) return
    end if
    call this%mesh%cell_imap%gather_offp(this%delta_p)

    call pressure_derivative(this, this%delta_p, bc, this%derivative_f, correction=.true.)
    pin_face = bc%pressure_pin_face()
    if (pin_face > 0) then
      c = this%mesh%fcell(1,pin_face)
      this%derivative_f(pin_face) = -this%delta_p(c)/this%operators%normal_distance(pin_face)
    end if
    do f = 1, this%mesh%nface_onP
      state%vel_fn(f) = state%vel_fn(f) - dt*inv_density_f(f)*this%derivative_f(f)
    end do
    call apply_velocity_boundary_conditions(this%mesh, bc, state%vel_fn)
    call this%mesh%face_imap%gather_offp(state%vel_fn)

    state%p_cc = state%p_cc + this%delta_p
    call this%mesh%cell_imap%gather_offp(state%p_cc)
    call gradient_pressure(this, state%p_cc, this%grad_p_new, bc, inv_density_c)
    do c = 1, this%mesh%ncell_onP
      state%vel_cc(:,c) = state%vel_cc(:,c) - dt*inv_density_c(c)* &
          (this%grad_p_new(:,c) - this%grad_p_old(:,c))
    end do
    call this%mesh%cell_imap%gather_offp(state%vel_cc)
  end subroutine


  subroutine gradient_pressure(this, pressure, gradient, bc, inv_density_c)
    class(flow_2d_projection_update), intent(inout) :: this
    real(r8), intent(in) :: pressure(:)
    real(r8), intent(out) :: gradient(:,:)
    type(flow_2d_bc), intent(in) :: bc
    real(r8), intent(in) :: inv_density_c(:)

    integer :: f, c1, c2
    real(r8) :: rho

    this%gravity_head = 0.0_r8
    do f = 1, this%mesh%nface_onP
      c1 = this%mesh%fcell(1,f)
      c2 = this%mesh%fcell(2,f)
      rho = 1.0_r8/inv_density_c(c1)
      if (allocated(this%buoyancy_temperature)) then
        rho = rho*(1.0_r8 - this%thermal_expan_coef*(this%buoyancy_temperature(c1) - &
            this%expan_ref_temp))
      end if
      this%gravity_head(1,f) = -this%body_acceleration(1)* &
          (this%mesh%cell_centroid(1,c1)-this%mesh%face_centroid(1,f))*rho - &
          this%body_acceleration(2)* &
          (this%mesh%cell_centroid(2,c1)-this%mesh%face_centroid(2,f))*rho
      if (c2 > 0) then
        rho = 1.0_r8/inv_density_c(c2)
        if (allocated(this%buoyancy_temperature)) then
          rho = rho*(1.0_r8 - this%thermal_expan_coef*(this%buoyancy_temperature(c2) - &
              this%expan_ref_temp))
        end if
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
            this%gravity_head)
      else
        call this%operators%gradient_cc(pressure, gradient, normal_flux_bc=bc%pressure_neumann, &
            gravity_head=this%gravity_head)
      end if
    else if (allocated(bc%pressure_dirichlet)) then
      call this%operators%gradient_cc(pressure, gradient, dirichlet_bc=bc%pressure_dirichlet, &
          gravity_head=this%gravity_head)
    else
      call this%operators%gradient_cc(pressure, gradient, gravity_head=this%gravity_head)
    end if
  end subroutine


  subroutine gradient_correction(this, pressure, gradient, bc)
    class(flow_2d_projection_update), intent(in) :: this
    real(r8), intent(in) :: pressure(:)
    real(r8), intent(out) :: gradient(:,:)
    type(flow_2d_bc), intent(in) :: bc

    call this%operators%gradient_cc(pressure, gradient, bc%pressure_neumann, &
        bc%pressure_correction_dirichlet)
  end subroutine


  subroutine pressure_derivative(this, pressure, bc, derivative, correction)
    class(flow_2d_projection_update), intent(inout) :: this
    real(r8), intent(in) :: pressure(:)
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
              bc%pressure_dirichlet, bc%pressure_correction_dirichlet%value)
        else
          call this%operators%derivative_cf(pressure, derivative, normal_flux_bc=bc%pressure_neumann)
        end if
      else if (allocated(bc%pressure_dirichlet)) then
        call this%operators%derivative_cf(pressure, derivative, &
            dirichlet_bc=bc%pressure_dirichlet, &
            dirichlet_value=bc%pressure_correction_dirichlet%value)
      else
        call this%operators%derivative_cf(pressure, derivative)
      end if
      return
    end if

    if (allocated(bc%pressure_neumann)) then
      if (allocated(bc%pressure_dirichlet)) then
        call this%operators%derivative_cf(pressure, derivative, bc%pressure_neumann, bc%pressure_dirichlet, &
            gravity_head=this%gravity_head)
      else
        call this%operators%derivative_cf(pressure, derivative, normal_flux_bc=bc%pressure_neumann, &
            gravity_head=this%gravity_head)
      end if
    else if (allocated(bc%pressure_dirichlet)) then
      call this%operators%derivative_cf(pressure, derivative, dirichlet_bc=bc%pressure_dirichlet, &
          gravity_head=this%gravity_head)
    else
      call this%operators%derivative_cf(pressure, derivative, gravity_head=this%gravity_head)
    end if
  end subroutine


  subroutine interpolate_velocity(this, velocity, bc, velocity_f)
    class(flow_2d_projection_update), intent(in) :: this
    real(r8), intent(in) :: velocity(:,:)
    type(flow_2d_bc), intent(in) :: bc
    real(r8), intent(out) :: velocity_f(:)

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
  end subroutine


  subroutine apply_velocity_boundary_conditions(mesh, bc, velocity_f)
    type(unstr_2d_mesh), intent(in) :: mesh
    type(flow_2d_bc), intent(in) :: bc
    real(r8), intent(inout) :: velocity_f(:)

    integer :: i, f

    if (allocated(bc%velocity_zero_normal)) then
      do i = 1, size(bc%velocity_zero_normal%index)
        f = bc%velocity_zero_normal%index(i)
        if (f > mesh%nface_onP) cycle
        velocity_f(f) = 0.0_r8
      end do
    end if
    if (allocated(bc%velocity_dirichlet)) then
      do i = 1, size(bc%velocity_dirichlet%index)
        f = bc%velocity_dirichlet%index(i)
        if (f > mesh%nface_onP) cycle
        velocity_f(f) = dot_product(mesh%unit_normal(:,f), bc%velocity_dirichlet%value(:,i))
      end do
    end if
  end subroutine

end module flow_2d_projection_update_type
