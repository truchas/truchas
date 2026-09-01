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
!! the trial material distribution changes DENSITY_C but does not change
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
  use flow_domain_types
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
        density_delta_c(:), inv_density_c(:), inv_density_f(:), viscosity_c(:), viscosity_f(:), &
        vof(:), vof_novoid(:)
    integer, allocatable, public :: cell_t(:), face_t(:)
    integer, public :: nfluid = 0
    real(r8), public :: cutoff = 0.01_r8
    logical, public :: any_void = .false., any_real_fluid = .false., any_real_fluid_onP = .false.
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
  subroutine init(this, mesh, density, inviscid, stat, errmsg, viscosity, viscosity_func, density_delta_func, nfluid)
    class(flow_2d_material_props), intent(out) :: this
    type(unstr_2d_mesh), target, intent(inout) :: mesh
    real(r8), intent(in) :: density(:)
    logical, intent(in) :: inviscid
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    real(r8), intent(in), optional :: viscosity
    class(scalar_func), allocatable, optional, intent(inout) :: viscosity_func, density_delta_func
    integer, intent(in), optional :: nfluid

    this%mesh => mesh
    if (size(density) == 0 .or. any(density <= 0.0_r8)) then
      stat = 1
      errmsg = 'fluid density must be positive'
      return
    end if
    this%density = density
    this%nfluid = size(density)
    if (present(nfluid)) this%nfluid = nfluid
    if (this%nfluid < size(density)) then
      stat = 1
      errmsg = 'number of mobile flow materials must include every real fluid'
      return
    end if
    allocate(this%vfrac(size(density), mesh%ncell), this%density_c(mesh%ncell), &
        this%density_c_old(mesh%ncell), this%density_delta_c(mesh%ncell), &
        this%inv_density_c(mesh%ncell), this%inv_density_f(mesh%nface), &
        this%density_delta(size(density)), this%vof(mesh%ncell), this%vof_novoid(mesh%ncell), &
        this%cell_t(mesh%ncell), this%face_t(mesh%nface))
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
  !! material. PHASE_IDS are phase indices in DENSITY order.
  subroutine init_material(this, mesh, matl_model, phase_ids, inviscid, boussinesq, stat, errmsg, nfluid)
    class(flow_2d_material_props), intent(out) :: this
    type(unstr_2d_mesh), target, intent(inout) :: mesh
    type(material_model), intent(in) :: matl_model
    integer, intent(in) :: phase_ids(:)
    logical, intent(in) :: inviscid
    logical, intent(in) :: boussinesq
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, intent(in), optional :: nfluid

    integer :: m
    real(r8) :: rho, alpha, tref
    class(scalar_func), allocatable :: f, alpha_func, tref_func

    this%mesh => mesh
    if (size(phase_ids) == 0) then
      stat = 1
      errmsg = 'flow requires at least one fluid phase'
      return
    end if
    if (any(phase_ids < 1) .or. any(phase_ids > matl_model%nphase_real)) then
      stat = 1
      errmsg = 'invalid flow phase index'
      return
    end if
    allocate(this%density(size(phase_ids)), this%vfrac(size(phase_ids), mesh%ncell), &
        this%density_c(mesh%ncell), this%density_c_old(mesh%ncell), this%density_delta_c(mesh%ncell), &
        this%inv_density_c(mesh%ncell), this%inv_density_f(mesh%nface), &
        this%density_delta(size(phase_ids)), this%vof(mesh%ncell), this%vof_novoid(mesh%ncell), &
        this%cell_t(mesh%ncell), this%face_t(mesh%nface))
    this%nfluid = size(phase_ids)
    if (present(nfluid)) this%nfluid = nfluid
    if (this%nfluid < size(phase_ids)) then
      stat = 1
      errmsg = 'number of mobile flow materials must include every real fluid'
      return
    end if
    if (.not.inviscid) then
      allocate(this%viscosity(size(phase_ids)), this%viscosity_c(mesh%ncell), this%viscosity_f(mesh%nface))
    end if

    do m = 1, size(phase_ids)
      call matl_model%get_phase_prop(phase_ids(m), 'density', f)
      if (.not.allocated(f) .or. .not.is_const(f)) then
        stat = 1
        errmsg = 'fluid density must be a constant property for ' // matl_model%phase_name(phase_ids(m))
        return
      end if
      rho = f%eval([real(r8)::])
      if (rho <= 0.0_r8) then
        stat = 1
        errmsg = 'fluid density must be positive for ' // matl_model%phase_name(phase_ids(m))
        return
      end if
      this%density(m) = rho

      if (.not.inviscid) then
        call matl_model%get_phase_prop(phase_ids(m), 'viscosity', this%viscosity(m)%f)
        if (.not.allocated(this%viscosity(m)%f)) then
          stat = 1
          errmsg = 'material viscosity property is missing for ' // matl_model%phase_name(phase_ids(m))
          return
        end if
      end if

      if (.not.boussinesq) then
        call alloc_const_scalar_func(this%density_delta(m)%f, 0.0_r8)
      else
        call matl_model%get_phase_prop(phase_ids(m), 'thermal-expan-coef', alpha_func)
        call matl_model%get_phase_prop(phase_ids(m), 'expan-ref-temp', tref_func)
        if (.not.allocated(alpha_func) .and. .not.allocated(tref_func)) then
          call alloc_const_scalar_func(this%density_delta(m)%f, 0.0_r8)
        else if (.not.allocated(alpha_func) .or. .not.allocated(tref_func) .or. &
            .not.is_const(alpha_func) .or. .not.is_const(tref_func)) then
          stat = 1
          errmsg = 'thermal expansion properties must be constant and provided together for ' // &
              matl_model%phase_name(phase_ids(m))
          return
        else
          alpha = alpha_func%eval([real(r8)::])
          tref = tref_func%eval([real(r8)::])
          if (alpha < 0.0_r8) then
            stat = 1
            errmsg = 'thermal-expan-coef must be nonnegative for ' // matl_model%phase_name(phase_ids(m))
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


  !! Update the current mixture density from the reduced flow distribution.
  !! The committed density is deliberately unchanged.
  subroutine set_volume_fractions(this, vfrac)
    class(flow_2d_material_props), intent(inout) :: this
    real(r8), intent(in) :: vfrac(:,:)

    integer :: c, c1, c2, f, j
    real(r8) :: weight(2), rho_f, weight_sum

    ASSERT(size(vfrac,1) >= this%nfluid)
    ASSERT(size(vfrac,2) >= this%mesh%ncell)
    this%vfrac = vfrac(:size(this%density),:this%mesh%ncell)
    this%vof = sum(vfrac(:this%nfluid,:this%mesh%ncell), dim=1)
    this%vof_novoid = sum(vfrac(:size(this%density),:this%mesh%ncell), dim=1)
    this%density_c = matmul(this%density, this%vfrac)
    where (this%vof > 0.0_r8) this%density_c = this%density_c/this%vof
    this%inv_density_c = 0.0_r8
    where (this%density_c > 0.0_r8) this%inv_density_c = 1.0_r8/this%density_c

    do c = 1, this%mesh%ncell
      if (this%vof(c) < this%cutoff) then
        this%cell_t(c) = solid_t
      else if (this%vof_novoid(c) == 0.0_r8) then
        this%cell_t(c) = void_t
      else if (this%vof(c) > this%vof_novoid(c)) then
        this%cell_t(c) = regular_void_t
      else
        this%cell_t(c) = regular_t
      end if
    end do

    do c = 1, this%mesh%ncell_onP
      if (this%cell_t(c) /= regular_void_t) cycle
      do j = this%mesh%cstart(c), this%mesh%cstart(c+1)-1
        if (this%mesh%cnhbr(j) > 0) then
          if (this%cell_t(this%mesh%cnhbr(j)) == void_t) then
            this%cell_t(c) = regular_t
            exit
          end if
        end if
      end do
    end do
    call this%mesh%cell_imap%gather_offp(this%cell_t)
    this%any_void = global_any(this%cell_t == void_t)
    this%any_real_fluid_onP = any(this%vof_novoid(:this%mesh%ncell_onP) > this%cutoff)
    this%any_real_fluid = global_any(this%any_real_fluid_onP)

    this%face_t = regular_t
    do f = 1, this%mesh%nface_onP
      c1 = this%mesh%fcell(1,f)
      c2 = this%mesh%fcell(2,f)
      if (c2 > 0) then
        if (any(this%cell_t([c1,c2]) == void_t) .and. &
            any(this%cell_t([c1,c2]) <= regular_t)) then
          this%face_t(f) = regular_void_t
        else
          this%face_t(f) = max(maxval(this%cell_t([c1,c2])), regular_t)
        end if
      else
        this%face_t(f) = max(this%cell_t(c1), regular_t)
      end if
    end do
    call this%mesh%face_imap%gather_offp(this%face_t)

    do f = 1, this%mesh%nface_onP
      c1 = this%mesh%fcell(1,f)
      c2 = this%mesh%fcell(2,f)
      if (c2 == 0) then
        if (this%cell_t(c1) <= regular_t) then
          rho_f = this%density_c(c1)
        else
          rho_f = 0.0_r8
        end if
      else if (any(this%cell_t([c1,c2]) <= regular_t)) then
        weight = this%mesh%volume([c1,c2])
        weight = weight*this%vof([c1,c2])
        weight_sum = sum(weight)
        if (weight_sum > 0.0_r8) then
          rho_f = dot_product(this%density_c([c1,c2]), weight)/weight_sum
        else
          rho_f = 0.0_r8
        end if
      else
        rho_f = 0.0_r8
      end if
      this%inv_density_f(f) = 0.0_r8
      if (rho_f > 0.0_r8) this%inv_density_f(f) = 1.0_r8/rho_f
    end do
    call this%mesh%face_imap%gather_offp(this%inv_density_f)
  end subroutine set_volume_fractions


  !! Evaluate temperature-dependent material properties for the current
  !! reduced distribution.
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
      if (this%vof(c) > 0.0_r8) then
        this%density_delta_c(c) = this%density_delta_c(c)/this%vof(c)
        if (allocated(this%viscosity_c)) this%viscosity_c(c) = this%viscosity_c(c)/this%vof(c)
      end if
    end do
    call this%mesh%cell_imap%gather_offp(this%density_delta_c)
    if (.not.allocated(this%viscosity_c)) return
    call this%mesh%cell_imap%gather_offp(this%viscosity_c)
    do f = 1, this%mesh%nface_onP
      c1 = this%mesh%fcell(1,f)
      c2 = this%mesh%fcell(2,f)
      if (c2 == 0) then
        this%viscosity_f(f) = this%viscosity_c(c1)
      else if (this%face_t(f) == solid_t) then
        this%viscosity_f(f) = maxval(this%viscosity_c([c1,c2]))
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
    this%vof = 1.0_r8
    this%vof_novoid = 1.0_r8
    this%cell_t = regular_t
    this%face_t = regular_t
    this%any_real_fluid = .true.
    this%any_real_fluid_onP = .true.
    call this%set_temperature([(0.0_r8, c=1,this%mesh%ncell_onP)])
  end subroutine initialize_values

end module flow_2d_material_props_type
