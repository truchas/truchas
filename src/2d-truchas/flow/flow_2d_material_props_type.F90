!!
!! FLOW_2D_MATERIAL_PROPS_TYPE
!!
!! This module defines FLOW_2D_MATERIAL_PROPS, the material-aware property
!! layer for the two-dimensional incompressible-flow model.  It constructs
!! flow properties from the material model and evaluates their cell- and
!! face-centered values from the reduced flow volume fractions and cell
!! temperatures.
!!
!! The current and committed cell densities are kept separately.  Updating
!! the trial material composition changes DENSITY_C but does not change
!! DENSITY_C_OLD; ACCEPT marks the current material state as committed.
!!
!! Neil Carlson <neil.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module flow_2d_material_props_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use scalar_func_class
  use scalar_func_containers, only: scalar_func_box
  use scalar_func_factories, only: alloc_const_scalar_func, alloc_poly_scalar_func
  use scalar_func_tools, only: is_const
  use unstr_2d_mesh_type
  use material_model_type
  use parallel_communication
  implicit none
  private

  type, public :: flow_2d_material_props
    private
    type(unstr_2d_mesh), pointer, public :: mesh => null()
    real(r8), allocatable, public :: density(:)
    type(scalar_func_box), allocatable :: viscosity(:), density_delta(:)
    real(r8), allocatable, public :: vfrac(:,:), density_c(:), density_c_old(:), &
        density_delta_c(:), inv_density_c(:), inv_density_f(:), viscosity_c(:), viscosity_f(:)
  contains
    procedure :: init
    procedure :: init_material
    procedure :: set_initial_state
    procedure :: set_volume_fractions
    procedure :: set_temperature
    procedure :: accept
  end type flow_2d_material_props

contains

  !! Initialize properties from already-constructed raw flow functions.
  !! This is retained for the standalone single-fluid flow driver.
  subroutine init(this, mesh, density, inviscid, stat, errmsg, viscosity, viscosity_func, density_delta_func)
    class(flow_2d_material_props), intent(out) :: this
    type(unstr_2d_mesh), target, intent(inout) :: mesh
    real(r8), intent(in) :: density(:)
    logical, intent(in) :: inviscid
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    real(r8), intent(in), optional :: viscosity
    class(scalar_func), allocatable, optional, intent(inout) :: viscosity_func, density_delta_func

    this%mesh => mesh
    if (size(density) == 0 .or. any(density <= 0.0_r8)) then
      stat = 1
      errmsg = 'fluid density must be positive'
      return
    end if
    this%density = density
    allocate(this%vfrac(size(density), mesh%ncell), this%density_c(mesh%ncell), &
        this%density_c_old(mesh%ncell), this%density_delta_c(mesh%ncell), &
        this%inv_density_c(mesh%ncell), this%inv_density_f(mesh%nface), &
        this%density_delta(size(density)))
    if (.not.inviscid) then
      allocate(this%viscosity(size(density)))
      if (size(density) /= 1) then
        stat = 1
        errmsg = 'raw flow properties support only one fluid material'
        return
      end if
      if (present(viscosity_func)) then
        call move_alloc(viscosity_func, this%viscosity(1)%f)
      else if (present(viscosity)) then
        call alloc_const_scalar_func(this%viscosity(1)%f, viscosity)
      else
        stat = 1
        errmsg = 'viscosity is required for a viscous flow model'
        return
      end if
      allocate(this%viscosity_c(mesh%ncell), this%viscosity_f(mesh%nface))
    end if
    if (present(density_delta_func)) then
      if (size(density) /= 1) then
        stat = 1
        errmsg = 'raw density-delta properties support only one fluid material'
        return
      end if
      call move_alloc(density_delta_func, this%density_delta(1)%f)
    else
      call alloc_const_scalar_func(this%density_delta(1)%f, 0.0_r8)
    end if
    this%vfrac = 0.0_r8
    this%vfrac(1,:) = 1.0_r8
    call initialize_values(this)
    stat = 0
    errmsg = ''
  end subroutine init


  !! Initialize properties by querying the material model for each flow
  !! material. MATERIAL_IDS are phase indices in DENSITY order.
  subroutine init_material(this, mesh, matl_model, material_ids, inviscid, boussinesq, stat, errmsg)
    class(flow_2d_material_props), intent(out) :: this
    type(unstr_2d_mesh), target, intent(inout) :: mesh
    type(material_model), intent(in) :: matl_model
    integer, intent(in) :: material_ids(:)
    logical, intent(in) :: inviscid
    logical, intent(in) :: boussinesq
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: m
    real(r8) :: rho, alpha, tref
    class(scalar_func), allocatable :: f, alpha_func, tref_func

    this%mesh => mesh
    if (size(material_ids) == 0) then
      stat = 1
      errmsg = 'flow requires at least one fluid material'
      return
    end if
    if (any(material_ids < 1) .or. any(material_ids > matl_model%nphase_real)) then
      stat = 1
      errmsg = 'invalid flow material phase index'
      return
    end if
    allocate(this%density(size(material_ids)), this%vfrac(size(material_ids), mesh%ncell), &
        this%density_c(mesh%ncell), this%density_c_old(mesh%ncell), this%density_delta_c(mesh%ncell), &
        this%inv_density_c(mesh%ncell), this%inv_density_f(mesh%nface), &
        this%density_delta(size(material_ids)))
    if (.not.inviscid) then
      allocate(this%viscosity(size(material_ids)), this%viscosity_c(mesh%ncell), this%viscosity_f(mesh%nface))
    end if

    do m = 1, size(material_ids)
      call matl_model%get_phase_prop(material_ids(m), 'density', f)
      if (.not.allocated(f) .or. .not.is_const(f)) then
        stat = 1
        errmsg = 'fluid density must be a constant property for ' // matl_model%phase_name(material_ids(m))
        return
      end if
      rho = f%eval([real(r8)::])
      if (rho <= 0.0_r8) then
        stat = 1
        errmsg = 'fluid density must be positive for ' // matl_model%phase_name(material_ids(m))
        return
      end if
      this%density(m) = rho

      if (.not.inviscid) then
        call matl_model%get_phase_prop(material_ids(m), 'viscosity', this%viscosity(m)%f)
        if (.not.allocated(this%viscosity(m)%f)) then
          stat = 1
          errmsg = 'material viscosity property is missing for ' // matl_model%phase_name(material_ids(m))
          return
        end if
      end if

      if (.not.boussinesq) then
        call alloc_const_scalar_func(this%density_delta(m)%f, 0.0_r8)
      else
        call matl_model%get_phase_prop(material_ids(m), 'thermal-expan-coef', alpha_func)
        call matl_model%get_phase_prop(material_ids(m), 'expan-ref-temp', tref_func)
        if (.not.allocated(alpha_func) .and. .not.allocated(tref_func)) then
          call alloc_const_scalar_func(this%density_delta(m)%f, 0.0_r8)
        else if (.not.allocated(alpha_func) .or. .not.allocated(tref_func) .or. &
            .not.is_const(alpha_func) .or. .not.is_const(tref_func)) then
          stat = 1
          errmsg = 'thermal expansion properties must be constant and provided together for ' // &
              matl_model%phase_name(material_ids(m))
          return
        else
          alpha = alpha_func%eval([real(r8)::])
          tref = tref_func%eval([real(r8)::])
          if (alpha < 0.0_r8) then
            stat = 1
            errmsg = 'thermal-expan-coef must be nonnegative for ' // matl_model%phase_name(material_ids(m))
            return
          end if
          call alloc_poly_scalar_func(this%density_delta(m)%f, [-rho*alpha], [1], tref)
        end if
      end if
    end do
    this%vfrac = 0.0_r8
    this%vfrac(1,:) = 1.0_r8
    call initialize_values(this)
    stat = 0
    errmsg = ''
  end subroutine init_material


  !! Initialize the current and committed material-property state.
  subroutine set_initial_state(this, vfrac, temperature)
    class(flow_2d_material_props), intent(inout) :: this
    real(r8), intent(in) :: vfrac(:,:), temperature(:)

    call this%set_volume_fractions(vfrac)
    call this%set_temperature(temperature)
    call this%accept()
  end subroutine set_initial_state


  !! Update the current mixture density from the reduced flow composition.
  !! The committed density is deliberately unchanged.
  subroutine set_volume_fractions(this, vfrac)
    class(flow_2d_material_props), intent(inout) :: this
    real(r8), intent(in) :: vfrac(:,:)

    integer :: c1, c2, f
    real(r8) :: weight(2), rho_f

    ASSERT(size(vfrac,1) >= size(this%density))
    ASSERT(size(vfrac,2) >= this%mesh%ncell)
    this%vfrac = vfrac(:size(this%density),:this%mesh%ncell)
    this%density_c = matmul(this%density, this%vfrac)
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
  end subroutine set_volume_fractions


  !! Evaluate temperature-dependent material properties for the current
  !! reduced composition.
  subroutine set_temperature(this, temperature)
    class(flow_2d_material_props), intent(inout) :: this

    real(r8), intent(in) :: temperature(:)
    integer :: c, c1, c2, f, m
    real(r8) :: state(1)

    ASSERT(size(temperature) == this%mesh%ncell_onP)
    do c = 1, this%mesh%ncell_onP
      state(1) = temperature(c)
      this%density_delta_c(c) = 0.0_r8
      if (allocated(this%viscosity_c)) this%viscosity_c(c) = 0.0_r8
      do m = 1, size(this%density)
        this%density_delta_c(c) = this%density_delta_c(c) + &
            this%vfrac(m,c)*this%density_delta(m)%f%eval(state)
        if (allocated(this%viscosity_c)) this%viscosity_c(c) = this%viscosity_c(c) + &
            this%vfrac(m,c)*this%viscosity(m)%f%eval(state)
      end do
    end do
    call this%mesh%cell_imap%gather_offp(this%density_delta_c)
    if (.not.allocated(this%viscosity_c)) return
    call this%mesh%cell_imap%gather_offp(this%viscosity_c)
    do f = 1, this%mesh%nface_onP
      c1 = this%mesh%fcell(1,f)
      c2 = this%mesh%fcell(2,f)
      if (c2 == 0) then
        this%viscosity_f(f) = this%viscosity_c(c1)
      else if (product(this%viscosity_c([c1,c2])) > epsilon(1.0_r8)) then
        this%viscosity_f(f) = 2.0_r8*product(this%viscosity_c([c1,c2])) / &
            sum(this%viscosity_c([c1,c2]))
      else
        this%viscosity_f(f) = 0.0_r8
      end if
    end do
    call this%mesh%face_imap%gather_offp(this%viscosity_f)
  end subroutine set_temperature


  !! Mark current cell density as the committed density for the next step.
  subroutine accept(this)
    class(flow_2d_material_props), intent(inout) :: this

    this%density_c_old = this%density_c
  end subroutine accept


  subroutine initialize_values(this)
    class(flow_2d_material_props), intent(inout) :: this

    integer :: c
    this%density_c = this%density(1)
    this%density_c_old = this%density_c
    this%inv_density_c = 1.0_r8/this%density(1)
    this%inv_density_f = 1.0_r8/this%density(1)
    call this%set_temperature([(0.0_r8, c=1,this%mesh%ncell_onP)])
  end subroutine initialize_values

end module flow_2d_material_props_type
