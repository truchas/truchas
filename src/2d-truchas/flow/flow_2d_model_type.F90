!!
!! FLOW_2D_MODEL_TYPE
!!
!! This module defines FLOW_2D_MODEL, the mesh-associated, isothermal
!! single-fluid model for a two-dimensional incompressible-flow calculation.
!! It owns the spatial operators, constant fluid properties, boundary
!! conditions, and discrete momentum and pressure-correction operators.  It
!! does not own a solution state or choose a time-integration algorithm.
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
  implicit none
  private

  type, public :: flow_2d_model
    private
    type(unstr_2d_mesh), pointer, public :: mesh => null()  ! unowned reference
    type(flow_2d_operators), pointer, public :: operators => null()
    type(flow_2d_bc), pointer, public :: bc => null()
    type(flow_2d_momentum), pointer, public :: momentum => null()
    type(flow_2d_projection), pointer, public :: projection => null()
    real(r8), allocatable, public :: density_c(:), inv_density_c(:), inv_density_f(:), viscosity_f(:)
  contains
    procedure :: init
    procedure :: compute_bc
    procedure :: pressure_gradient
  end type

contains

  subroutine init(this, mesh, bc_params, density, viscosity, stat, errmsg)
    class(flow_2d_model), intent(out) :: this
    type(unstr_2d_mesh), target, intent(inout) :: mesh
    type(parameter_list), target, intent(inout) :: bc_params
    real(r8), intent(in) :: density, viscosity
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    this%mesh => mesh
    allocate(this%operators, this%bc, this%momentum, this%projection)
    allocate(this%density_c(mesh%ncell), this%inv_density_c(mesh%ncell), &
        this%inv_density_f(mesh%nface), this%viscosity_f(mesh%nface))
    this%density_c = density
    this%inv_density_c = 1.0_r8/density
    this%inv_density_f = 1.0_r8/density
    this%viscosity_f = viscosity

    call this%operators%init(mesh)
    call this%bc%init(mesh, bc_params, stat, errmsg)
    if (stat /= 0) return
    call this%momentum%init(mesh, this%operators)
    call this%projection%init(mesh, this%operators)
  end subroutine


  subroutine compute_bc(this, time, dt)
    class(flow_2d_model), intent(inout) :: this
    real(r8), intent(in) :: time, dt

    call this%bc%compute(time, dt)
  end subroutine


  !! Compute the cell pressure gradient using the physical pressure boundary
  !! conditions evaluated by COMPUTE_BC.
  subroutine pressure_gradient(this, pressure, gradient)
    class(flow_2d_model), intent(in) :: this
    real(r8), intent(in) :: pressure(:)
    real(r8), intent(out) :: gradient(:,:)

    call this%operators%gradient_cc(pressure, gradient, this%bc%pressure_neumann, this%bc%pressure_dirichlet)
  end subroutine

end module flow_2d_model_type
