!!
!! FLOW_2D_MODEL_TYPE
!!
!! This module defines FLOW_2D_MODEL, the mesh-associated, isothermal
!! single-fluid model for a two-dimensional incompressible-flow calculation.
!! It owns the spatial operators, constant fluid properties, boundary
!! conditions, and discrete momentum and pressure-correction operators.  It
!! also stores the uniform body acceleration used by the effective pressure
!! gradient and, optionally, a Boussinesq thermal-expansion coefficient.  It
!! does not own a solution state or choose a time-integration algorithm.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module flow_2d_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use unstr_2d_mesh_type
  use flow_2d_operators_type
  use flow_2d_bc_type
  use flow_2d_momentum_type
  use flow_2d_projection_type
  use simulation_environment_type
  implicit none
  private

  type, public :: flow_2d_model
    private
    type(unstr_2d_mesh), pointer, public :: mesh => null()  ! unowned reference
    type(flow_2d_operators), pointer, public :: operators => null()
    type(flow_2d_bc), pointer, public :: bc => null()
    type(flow_2d_momentum), pointer, public :: momentum => null()
    type(flow_2d_projection), pointer, public :: projection => null()
    real(r8), allocatable, public :: density(:), density_c(:), inv_density_c(:), inv_density_f(:), viscosity_f(:)
    real(r8), public :: body_acceleration(2) = 0.0_r8
    real(r8), public :: thermal_expan_coef = 0.0_r8
    real(r8), public :: expan_ref_temp = 0.0_r8
    real(r8), allocatable :: buoyancy_temperature(:)
  contains
    procedure :: init
    procedure :: compute_bc
    procedure :: set_buoyancy_temperature
    procedure :: pressure_gradient
  end type

contains

  subroutine init(this, env, mesh, bc_params, density, viscosity, stat, errmsg, body_acceleration, &
      thermal_expan_coef, expan_ref_temp)
    class(flow_2d_model), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(unstr_2d_mesh), target, intent(inout) :: mesh
    type(parameter_list), target, intent(inout) :: bc_params
    real(r8), intent(in) :: density, viscosity
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    real(r8), optional, intent(in) :: body_acceleration(:)
    real(r8), optional, intent(in) :: thermal_expan_coef, expan_ref_temp

    stat = 0
    if (present(thermal_expan_coef)) then
      if (thermal_expan_coef < 0.0_r8) then
        stat = 1
        errmsg = 'thermal expansion coefficient must be nonnegative'
        return
      end if
      this%thermal_expan_coef = thermal_expan_coef
    end if
    if (present(expan_ref_temp)) this%expan_ref_temp = expan_ref_temp
    if (present(body_acceleration)) then
      if (size(body_acceleration) /= 2) then
        stat = 1
        errmsg = 'body acceleration must have two components'
        return
      end if
      this%body_acceleration = body_acceleration
    end if
    this%mesh => mesh
    allocate(this%operators, this%bc, this%momentum, this%projection)
    allocate(this%density(1), this%density_c(mesh%ncell), this%inv_density_c(mesh%ncell), &
        this%inv_density_f(mesh%nface), this%viscosity_f(mesh%nface))
    this%density = density
    this%density_c = this%density(1)
    this%inv_density_c = 1.0_r8/this%density(1)
    this%inv_density_f = 1.0_r8/this%density(1)
    this%viscosity_f = viscosity
    if (this%thermal_expan_coef > 0.0_r8) then
      allocate(this%buoyancy_temperature(mesh%ncell), source=0.0_r8)
    end if

    call this%operators%init(mesh)
    call this%bc%init(env, mesh, bc_params, stat, errmsg)
    if (stat /= 0) return
    call this%momentum%init(mesh, this%operators)
    call this%projection%init(mesh, this%operators)
  end subroutine


  !! Set the cell temperature used to evaluate the Boussinesq buoyancy term.
  !! The reference temperature determines the pressure potential associated
  !! with the constant part of the buoyancy force.
  subroutine set_buoyancy_temperature(this, temperature)
    class(flow_2d_model), intent(inout) :: this
    real(r8), intent(in) :: temperature(:)

    if (.not.allocated(this%buoyancy_temperature)) return
    ASSERT(size(temperature) == this%mesh%ncell_onP)
    this%buoyancy_temperature(:this%mesh%ncell_onP) = temperature
    call this%mesh%cell_imap%gather_offp(this%buoyancy_temperature)
  end subroutine


  subroutine compute_bc(this, time, dt, stat, errmsg)
    class(flow_2d_model), intent(inout) :: this
    real(r8), intent(in) :: time, dt
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    stat = 0
    call this%bc%compute(time, dt)
    call this%bc%check_velocity_flux(stat, errmsg)
  end subroutine


  !! Compute the cell pressure gradient using the physical pressure boundary
  !! conditions evaluated by COMPUTE_BC.
  subroutine pressure_gradient(this, pressure, gradient)
    class(flow_2d_model), intent(in) :: this
    real(r8), intent(in) :: pressure(:)
    real(r8), intent(out) :: gradient(:,:)

    integer :: f, c1, c2
    real(r8) :: gravity_head(2,this%mesh%nface), rho

    gravity_head = 0.0_r8
    do f = 1, this%mesh%nface_onP
      c1 = this%mesh%fcell(1,f)
      c2 = this%mesh%fcell(2,f)
      rho = this%density_c(c1)
      if (allocated(this%buoyancy_temperature)) then
        rho = rho*(1.0_r8 - this%thermal_expan_coef*(this%buoyancy_temperature(c1) - &
            this%expan_ref_temp))
      end if
      gravity_head(1,f) = -rho*dot_product(this%body_acceleration, &
          this%mesh%cell_centroid(:,c1) - this%mesh%face_centroid(:,f))
      if (c2 > 0) then
        rho = this%density_c(c2)
        if (allocated(this%buoyancy_temperature)) then
          rho = rho*(1.0_r8 - this%thermal_expan_coef*(this%buoyancy_temperature(c2) - &
              this%expan_ref_temp))
        end if
        gravity_head(2,f) = -rho*dot_product(this%body_acceleration, &
            this%mesh%cell_centroid(:,c2) - this%mesh%face_centroid(:,f))
      end if
    end do
    call this%mesh%face_imap%gather_offp(gravity_head)
    if (allocated(this%bc%pressure_neumann)) then
      if (allocated(this%bc%pressure_dirichlet)) then
        call this%operators%gradient_cc(pressure, gradient, this%bc%pressure_neumann, &
            this%bc%pressure_dirichlet, gravity_head)
      else
        call this%operators%gradient_cc(pressure, gradient, &
            normal_flux_bc=this%bc%pressure_neumann, gravity_head=gravity_head)
      end if
    else if (allocated(this%bc%pressure_dirichlet)) then
      call this%operators%gradient_cc(pressure, gradient, &
          dirichlet_bc=this%bc%pressure_dirichlet, gravity_head=gravity_head)
    else
      call this%operators%gradient_cc(pressure, gradient, gravity_head=gravity_head)
    end if
  end subroutine

end module flow_2d_model_type
