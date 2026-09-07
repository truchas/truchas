!!
!! T2D_FLOW_MODEL_TYPE
!!
!! This module defines T2D_FLOW_MODEL, the mesh-associated material-property
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

module t2d_flow_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use scalar_func_class
  use t2d_unstr_mesh_type
  use material_model_type
  use t2d_flow_operators_type
  use t2d_flow_bc_type
  use t2d_flow_momentum_type
  use t2d_flow_projection_type
  use t2d_flow_material_props_type
  use flow_domain_types
  use simulation_environment_type
  implicit none
  private

  type, public :: t2d_flow_model
    private
    type(t2d_unstr_mesh), pointer, public :: mesh => null()  ! unowned reference
    type(t2d_flow_operators), pointer, public :: operators => null()
    type(t2d_flow_bc), pointer, public :: bc => null()
    type(t2d_flow_momentum), pointer, public :: momentum => null()
    type(t2d_flow_projection), pointer, public :: projection => null()
    type(t2d_flow_material_props), public :: matl_props
    logical, public :: inviscid = .false.
    real(r8), public :: body_acceleration(2) = 0.0_r8
  contains
    procedure :: init
    procedure :: init_core
    procedure :: init_material
    procedure :: set_fluid_fraction_cutoff
    procedure :: set_min_face_fraction
    procedure :: set_volume_fractions
    procedure :: set_initial_material_state
    procedure :: set_pre_solidification_state
    procedure :: accept_material_state
    procedure :: compute_bc
    procedure :: set_buoyancy_temperature
    procedure :: assemble_momentum
  end type

contains

  subroutine init(this, env, mesh, bc_params, density, viscosity, stat, errmsg, body_acceleration, &
      viscosity_func, density_delta_func, inviscid)
    class(t2d_flow_model), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(t2d_unstr_mesh), target, intent(inout) :: mesh
    type(parameter_list), target, intent(inout) :: bc_params
    real(r8), intent(in) :: density(:)
    real(r8), optional, intent(in) :: viscosity
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    real(r8), optional, intent(in) :: body_acceleration(:)
    class(scalar_func), allocatable, optional, intent(inout) :: viscosity_func, density_delta_func
    logical, optional, intent(in) :: inviscid
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
    if (present(viscosity_func)) then
      if (present(density_delta_func)) then
        call this%matl_props%init(mesh, density, this%inviscid, stat, errmsg, viscosity_func=viscosity_func, &
            density_delta_func=density_delta_func)
      else
        call this%matl_props%init(mesh, density, this%inviscid, stat, errmsg, viscosity_func=viscosity_func)
      end if
    else if (present(viscosity)) then
      if (present(density_delta_func)) then
        call this%matl_props%init(mesh, density, this%inviscid, stat, errmsg, viscosity=viscosity, &
            density_delta_func=density_delta_func)
      else
        call this%matl_props%init(mesh, density, this%inviscid, stat, errmsg, viscosity=viscosity)
      end if
    else if (present(density_delta_func)) then
      call this%matl_props%init(mesh, density, this%inviscid, stat, errmsg, &
          density_delta_func=density_delta_func)
    else
      call this%matl_props%init(mesh, density, this%inviscid, stat, errmsg)
    end if
    if (stat /= 0) return
    call check_initial_properties(this, mesh, stat, errmsg)
    if (stat /= 0) return
    call this%init_core(env, mesh, bc_params, stat, errmsg)
  end subroutine


  !! Initialize the model's material properties directly from fluid phases.
  !! The model core must already have been initialized.
  subroutine init_material(this, matl_model, material_ids, stat, errmsg, boussinesq, nfluid)
    class(t2d_flow_model), intent(inout) :: this
    type(material_model), intent(in) :: matl_model
    integer, intent(in) :: material_ids(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    logical, optional, intent(in) :: boussinesq
    integer, optional, intent(in) :: nfluid

    logical :: use_boussinesq
    real(r8) :: fluid_fraction_cutoff, min_face_fraction

    stat = 0
    use_boussinesq = .false.
    if (present(boussinesq)) use_boussinesq = boussinesq
    fluid_fraction_cutoff = this%matl_props%cutoff
    min_face_fraction = this%matl_props%min_face_fraction
    if (present(nfluid)) then
      call this%matl_props%init_material(this%mesh, matl_model, material_ids, this%inviscid, use_boussinesq, &
          stat, errmsg, nfluid)
    else
      call this%matl_props%init_material(this%mesh, matl_model, material_ids, this%inviscid, use_boussinesq, &
          stat, errmsg)
    end if
    this%matl_props%cutoff = fluid_fraction_cutoff
    this%matl_props%min_face_fraction = min_face_fraction
    if (stat == 0) call check_initial_properties(this, this%mesh, stat, errmsg)
  end subroutine


  !! Initialize the mesh-associated operators and boundary conditions.  The
  !! material properties are initialized separately by INIT_MATERIAL.
  subroutine init_core(this, env, mesh, bc_params, stat, errmsg, body_acceleration, inviscid)
    class(t2d_flow_model), intent(inout) :: this
    type(simulation_environment), intent(in) :: env
    type(t2d_unstr_mesh), target, intent(inout) :: mesh
    type(parameter_list), target, intent(inout) :: bc_params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    real(r8), optional, intent(in) :: body_acceleration(:)
    logical, optional, intent(in) :: inviscid

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
    call this%operators%init(mesh)
    call this%bc%init(env, mesh, bc_params, stat, errmsg)
    if (stat /= 0) return
    call this%momentum%init(mesh, this%operators, this%inviscid)
    call this%projection%init(mesh, this%operators)
  end subroutine init_core


  subroutine check_initial_properties(this, mesh, stat, errmsg)
    class(t2d_flow_model), intent(in) :: this
    type(t2d_unstr_mesh), intent(in) :: mesh
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    stat = 0
    if (allocated(this%matl_props%viscosity_c)) then
      if (any(this%matl_props%viscosity_c(:mesh%ncell_onP) <= 0.0_r8)) then
        stat = 1
        errmsg = 'viscosity must be positive at the initial temperature'
      end if
    end if
  end subroutine check_initial_properties


  !! Set the local material volume fractions used to form mixture density and
  !! the face inverse densities used by the pressure projection.  The leading
  !! rows correspond to the real fluid materials in DENSITY order.
  subroutine set_volume_fractions(this, vfrac)
    class(t2d_flow_model), intent(inout) :: this
    real(r8), intent(in) :: vfrac(:,:)

    call this%matl_props%set_volume_fractions(vfrac)
  end subroutine


  !! Set the mobile-fluid fraction below which a cell receives a dummy flow
  !! equation.  This affects only the flow solve; material transport retains
  !! the corresponding volume fractions.
  subroutine set_fluid_fraction_cutoff(this, cutoff, stat, errmsg)
    class(t2d_flow_model), intent(inout) :: this
    real(r8), intent(in) :: cutoff
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    stat = 0
    if (cutoff <= 0.0_r8 .or. cutoff >= 1.0_r8) then
      stat = 1
      errmsg = '"fluid-fraction-cutoff" must be in (0,1)'
      return
    end if
    this%matl_props%cutoff = cutoff
  end subroutine


  !! Set the minimum nonzero interior face density as a fraction of the
  !! lightest real-fluid density.  This limits the projection coefficient at
  !! fluid/VOID interfaces without altering the material distribution.
  subroutine set_min_face_fraction(this, fraction, stat, errmsg)
    class(t2d_flow_model), intent(inout) :: this
    real(r8), intent(in) :: fraction
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    stat = 0
    if (fraction <= 0.0_r8 .or. fraction > 1.0_r8) then
      stat = 1
      errmsg = '"min-face-fraction" must be in (0,1]'
      return
    end if
    this%matl_props%min_face_fraction = fraction
  end subroutine


  subroutine set_initial_material_state(this, vfrac, temperature)
    class(t2d_flow_model), intent(inout) :: this
    real(r8), intent(in) :: vfrac(:,:), temperature(:)

    call this%matl_props%set_initial_state(vfrac, temperature)
  end subroutine


  !! Save the mobile-fluid mass density before thermal phase change.
  subroutine set_pre_solidification_state(this)
    class(t2d_flow_model), intent(inout) :: this

    call this%matl_props%set_pre_solidification_state()
  end subroutine set_pre_solidification_state


  subroutine accept_material_state(this)
    class(t2d_flow_model), intent(inout) :: this

    call this%matl_props%accept()
  end subroutine


  !! Set the cell temperature used to evaluate the material properties and
  !! the Boussinesq buoyancy term.
  subroutine set_buoyancy_temperature(this, temperature)
    class(t2d_flow_model), intent(inout) :: this
    real(r8), intent(in) :: temperature(:)

    call this%matl_props%set_temperature(temperature)
  end subroutine


  !! Assemble the momentum predictor operator for the current flow model.
  !! Inviscid flow has only the cell mass blocks; viscous flow also includes
  !! diffusion and velocity-boundary contributions.
  subroutine assemble_momentum(this, dt, rhs)
    class(t2d_flow_model), intent(inout) :: this
    real(r8), intent(in) :: dt
    real(r8), intent(out) :: rhs(:,:)

    if (this%inviscid) then
      call this%momentum%assemble_inviscid(this%matl_props%density_c, this%matl_props%cell_t, rhs, &
          this%matl_props%solidified_density, this%matl_props%vof)
    else
      call this%momentum%assemble(dt, this%matl_props%density_c, this%matl_props%viscosity_f, &
          this%matl_props%cell_t, this%matl_props%face_t, this%bc, rhs, &
          this%matl_props%solidified_density, this%matl_props%vof)
    end if
  end subroutine


  subroutine compute_bc(this, time, dt, stat, errmsg)
    class(t2d_flow_model), intent(inout) :: this
    real(r8), intent(in) :: time, dt
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    stat = 0
    call this%bc%compute(time, dt)
    call this%bc%check_velocity_flux(stat, errmsg)
  end subroutine


end module t2d_flow_model_type
