!!
!! FLOW_2D_MODEL_TYPE
!!
!! This module defines FLOW_2D_MODEL, the mesh-associated material-property
!! model for a two-dimensional incompressible-flow calculation.
!! It owns the spatial operators, fluid properties, boundary
!! conditions, and discrete momentum and pressure-correction operators.  It
!! also stores the uniform body acceleration and the temperature-dependent
!! density and viscosity values used by the operators.  It does not own a
!! solution state or choose a time-integration algorithm.
!! In inviscid mode the viscosity property and its derived arrays are omitted,
!! and the momentum operator uses only the cell mass blocks.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module flow_2d_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use scalar_func_class
  use scalar_func_factories, only: alloc_const_scalar_func
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
    logical, public :: inviscid = .false.
    real(r8), allocatable, public :: density(:), density_c(:), density_delta_c(:), inv_density_c(:), &
        inv_density_f(:), viscosity_c(:), viscosity_f(:)
    real(r8), public :: body_acceleration(2) = 0.0_r8
    class(scalar_func), allocatable :: viscosity, density_delta
  contains
    procedure :: init
    procedure :: set_volume_fractions
    procedure :: compute_bc
    procedure :: set_buoyancy_temperature
    procedure :: pressure_gradient
    procedure :: assemble_momentum
  end type

contains

  subroutine init(this, env, mesh, bc_params, density, viscosity, stat, errmsg, body_acceleration, &
      viscosity_func, density_delta_func, inviscid)
    class(flow_2d_model), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(unstr_2d_mesh), target, intent(inout) :: mesh
    type(parameter_list), target, intent(inout) :: bc_params
    real(r8), intent(in) :: density(:)
    real(r8), optional, intent(in) :: viscosity
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    real(r8), optional, intent(in) :: body_acceleration(:)
    class(scalar_func), allocatable, optional, intent(inout) :: viscosity_func, density_delta_func
    logical, optional, intent(in) :: inviscid
    integer :: c

    stat = 0
    if (present(inviscid)) this%inviscid = inviscid
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
    allocate(this%density(size(density)), this%density_c(mesh%ncell), this%density_delta_c(mesh%ncell), &
        this%inv_density_c(mesh%ncell), this%inv_density_f(mesh%nface))
    if (.not.this%inviscid) allocate(this%viscosity_c(mesh%ncell), this%viscosity_f(mesh%nface))
    if (size(density) == 0 .or. any(density <= 0.0_r8)) then
      stat = 1
      errmsg = 'fluid density must be positive'
      return
    end if
    this%density = density
    this%density_c = this%density(1)
    this%density_delta_c = 0.0_r8
    this%inv_density_c = 1.0_r8/this%density(1)
    this%inv_density_f = 1.0_r8/this%density(1)
    if (.not.this%inviscid) then
      if (present(viscosity_func)) then
        call move_alloc(viscosity_func, this%viscosity)
      else if (present(viscosity)) then
        call alloc_const_scalar_func(this%viscosity, viscosity)
      else
        stat = 1
        errmsg = 'viscosity is required for a viscous flow model'
        return
      end if
    end if
    if (present(density_delta_func)) then
      call move_alloc(density_delta_func, this%density_delta)
    else
      call alloc_const_scalar_func(this%density_delta, 0.0_r8)
    end if

    call this%operators%init(mesh)
    call this%set_buoyancy_temperature([(0.0_r8, c=1,mesh%ncell_onP)])
    if (allocated(this%viscosity_c)) then
      if (any(this%viscosity_c(:mesh%ncell_onP) <= 0.0_r8)) then
        stat = 1
        errmsg = 'viscosity must be positive at the initial temperature'
        return
      end if
    end if
    call this%bc%init(env, mesh, bc_params, stat, errmsg)
    if (stat /= 0) return
    call this%momentum%init(mesh, this%operators, this%inviscid)
    call this%projection%init(mesh, this%operators)
  end subroutine


  !! Set the local material volume fractions used to form mixture density and
  !! the face inverse densities used by the pressure projection.  The leading
  !! rows correspond to the real fluid materials in DENSITY order.
  subroutine set_volume_fractions(this, vfrac)
    class(flow_2d_model), intent(inout) :: this
    real(r8), intent(in) :: vfrac(:,:)

    integer :: c1, c2, f
    real(r8) :: weight(2), rho_f

    ASSERT(size(vfrac,1) >= size(this%density))
    ASSERT(size(vfrac,2) >= this%mesh%ncell)
    this%density_c = matmul(this%density, vfrac(:size(this%density),:this%mesh%ncell))
    this%inv_density_c = 1.0_r8/this%density_c
    do f = 1, this%mesh%nface_onP
      c1 = this%mesh%fcell(1,f)
      c2 = this%mesh%fcell(2,f)
      if (c2 == 0) then
        rho_f = this%density_c(c1)
      else
        weight = this%mesh%volume([c1,c2])
        rho_f = dot_product(this%density_c([c1,c2]), weight)/sum(weight)
      end if
      this%inv_density_f(f) = 1.0_r8/rho_f
    end do
    call this%mesh%face_imap%gather_offp(this%inv_density_f)
  end subroutine


  !! Set the cell temperature used to evaluate the material properties and
  !! the Boussinesq buoyancy term.
  subroutine set_buoyancy_temperature(this, temperature)
    class(flow_2d_model), intent(inout) :: this
    real(r8), intent(in) :: temperature(:)

    ASSERT(size(temperature) == this%mesh%ncell_onP)
    block
      integer :: c
      real(r8) :: state(1)
      do c = 1, this%mesh%ncell_onP
        state(1) = temperature(c)
        this%density_delta_c(c) = this%density_delta%eval(state)
        if (allocated(this%viscosity)) this%viscosity_c(c) = this%viscosity%eval(state)
      end do
    end block
    call this%mesh%cell_imap%gather_offp(this%density_delta_c)
    if (.not.allocated(this%viscosity_c)) return
    call this%mesh%cell_imap%gather_offp(this%viscosity_c)
    block
      integer :: c1, c2, f
      do f = 1, this%mesh%nface_onP
        c1 = this%mesh%fcell(1,f)
        c2 = this%mesh%fcell(2,f)
        if (c2 == 0) then
          this%viscosity_f(f) = this%viscosity_c(c1)
        else
          this%viscosity_f(f) = 2.0_r8*this%viscosity_c(c1)*this%viscosity_c(c2) / &
              (this%viscosity_c(c1) + this%viscosity_c(c2))
        end if
      end do
    end block
    call this%mesh%face_imap%gather_offp(this%viscosity_f)
  end subroutine


  !! Assemble the momentum predictor operator for the current flow model.
  !! Inviscid flow has only the cell mass blocks; viscous flow also includes
  !! diffusion and velocity-boundary contributions.
  subroutine assemble_momentum(this, dt, rhs)
    class(flow_2d_model), intent(inout) :: this
    real(r8), intent(in) :: dt
    real(r8), intent(out) :: rhs(:,:)

    if (this%inviscid) then
      call this%momentum%assemble_inviscid(this%density_c, rhs)
    else
      call this%momentum%assemble(dt, this%density_c, this%viscosity_f, this%bc, rhs)
    end if
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
      rho = this%density_c(c1) + this%density_delta_c(c1)
      gravity_head(1,f) = -rho*dot_product(this%body_acceleration, &
          this%mesh%cell_centroid(:,c1) - this%mesh%face_centroid(:,f))
      if (c2 > 0) then
        rho = this%density_c(c2) + this%density_delta_c(c2)
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
