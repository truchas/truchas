!!
!! FLOW_2D_MODEL_TYPE
!!
!! This module defines FLOW_2D_MODEL, the mesh-associated, isothermal
!! single-fluid model for a two-dimensional incompressible-flow calculation.
!! It owns the spatial operators, constant fluid properties, boundary
!! conditions, and discrete momentum and pressure-correction operators.  It
!! also stores the uniform body acceleration used by the effective pressure
!! gradient.  It does not own a solution state or choose a time-integration
!! algorithm.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

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
  contains
    procedure :: init
    procedure :: compute_bc
    procedure :: pressure_gradient
  end type

contains

  subroutine init(this, env, mesh, bc_params, density, viscosity, stat, errmsg, body_acceleration)
    class(flow_2d_model), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(unstr_2d_mesh), target, intent(inout) :: mesh
    type(parameter_list), target, intent(inout) :: bc_params
    real(r8), intent(in) :: density, viscosity
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    real(r8), optional, intent(in) :: body_acceleration(:)

    stat = 0
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

    call this%operators%init(mesh)
    call this%bc%init(env, mesh, bc_params, stat, errmsg)
    if (stat /= 0) return
    call this%momentum%init(mesh, this%operators)
    call this%projection%init(mesh, this%operators)
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
    real(r8) :: gravity_head(2,this%mesh%nface)

    gravity_head = 0.0_r8
    do f = 1, this%mesh%nface_onP
      c1 = this%mesh%fcell(1,f)
      c2 = this%mesh%fcell(2,f)
      gravity_head(1,f) = -this%density_c(c1)*dot_product(this%body_acceleration, &
          this%mesh%cell_centroid(:,c1) - this%mesh%face_centroid(:,f))
      if (c2 > 0) gravity_head(2,f) = -this%density_c(c2)*dot_product(this%body_acceleration, &
          this%mesh%cell_centroid(:,c2) - this%mesh%face_centroid(:,f))
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
