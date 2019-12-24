!!
!! SINGLE_PHASE_MATL_TYPE
!!
!! This module provides the concrete implementation SINGLE_PHASE_MATL of
!! the MATERIAL base class that defines a simple single-phase material.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! November 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module single_phase_matl_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use material_class
  use matl_prop_class
  use scalar_func_class
  implicit none
  private

  type, extends(material), public :: single_phase_matl
  contains
    ! Deferred procedures from MATERIAL
    procedure :: alloc_matl_prop
    procedure :: num_phase
    procedure :: phase_name
    procedure :: phase_ref
    procedure :: get_phase_frac
    procedure :: create_enthalpy
  end type single_phase_matl

  type, extends(matl_prop) :: simple_prop
    class(scalar_func), allocatable :: func
  contains
    ! Deferred procedure from MATL_PROP
    procedure :: compute_value
  end type

contains

  integer function num_phase(this)
    class(single_phase_matl), intent(in) :: this
    num_phase = 1
  end function

  function phase_name(this, n) result(name)
    class(single_phase_matl), intent(in) :: this
    integer, intent(in) :: n  ! ignored, should be 1
    character(:), allocatable :: name
    name = this%name
  end function

  function phase_ref(this, n) result(phi)
    class(single_phase_matl), intent(in), target :: this
    integer, intent(in) :: n  ! ignored, should be 1
    type(phase), pointer :: phi
    phi => this%phase ! cast to the parent class
  end function

  !! This method really should not be called for this type.
  subroutine get_phase_frac(this, temp, beta)
    class(single_phase_matl), intent(in) :: this
    real(r8), intent(in) :: temp
    real(r8), intent(out) :: beta(:)
    beta(1) = 1
  end subroutine get_phase_frac

  !! If the material does not already have the specific enthalpy property,
  !! this subroutine attempts to create it from the specific heat property
  !! by integration, if possible. STAT returns a non-zero value if the
  !! specific enthalpy property is not ultimately defined, and ERRMSG
  !! returns an explanatory message.
  !! TODO: Is there some way to make this external to the class?

  subroutine create_specific_enthalpy(this, stat, errmsg)

    use scalar_func_tools

    class(single_phase_matl), intent(inout) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    class(scalar_func), allocatable :: cp, h, rho, hd
    real(r8) :: t0, h0

    if (this%has_prop('specific-enthalpy')) return

    if (this%phase%has_prop('ref-temp')) then
      INSIST(this%phase%has_const_prop('ref-temp'))
      call this%phase%get_prop('ref-temp', t0)
    else
      t0 = 0.0_r8
    end if

    if (this%phase%has_prop('ref-enthalpy')) then
      INSIST(this%phase%has_const_prop('ref-enthalpy'))
      call this%phase%get_prop('ref-enthalpy', h0)
    else
      h0 = 0.0_r8
    end if

    call this%get_prop('specific-heat', cp)
    if (.not.allocated(cp)) then
      stat = 1
      errmsg = 'missing material specific heat'
      return
    end if

    call alloc_scalar_func_antideriv(cp, t0, h0, h, stat, errmsg)
    if (stat /= 0) then
      errmsg = 'unable to integrate specific heat: ' // errmsg
      return
    end if
    call this%add_prop('specific-enthalpy', h)

  end subroutine create_specific_enthalpy

  !! If the material does not already have the enthalpy property, this
  !! subroutine attempts to create it from the specific enthalpy property
  !! by forming its product with the density, if possible. If the specific
  !! enthalpy property does not exist, it will try to create it too, if
  !! possible, using the specific heat property. STAT returns a non-zero
  !! value if the enthalpy property is not ultimately defined, and ERRMSG
  !! returns an explanatory message.
  !! TODO: Is there some way to make this external to the class?

  subroutine create_enthalpy(this, stat, errmsg)

    use scalar_func_tools

    class(single_phase_matl), intent(inout) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    class(scalar_func), allocatable :: h, rho, hd

    if (this%has_prop('enthalpy')) return

    if (.not.this%has_prop('specific-enthalpy')) then
      call create_specific_enthalpy(this, stat, errmsg)
      if (stat /= 0) return
    end if

    call this%get_prop('specific-enthalpy', h)
    call this%get_prop('density', rho)
    if (.not.allocated(rho)) then
      stat = 1
      errmsg = 'missing material density'
      return
    end if
    call alloc_scalar_func_product(rho, h, hd, stat, errmsg)
    INSIST(stat == 0)
    call this%add_prop('enthalpy', hd)

  end subroutine create_enthalpy

  subroutine alloc_matl_prop(this, name, prop, errmsg)

    class(single_phase_matl), intent(in), target :: this
    character(*), intent(in) :: name
    class(matl_prop), allocatable, intent(out) :: prop
    character(:), allocatable, intent(out) :: errmsg

    type(simple_prop), allocatable :: p
    allocate(p)
    call this%get_prop(name, p%func)
    if (.not.allocated(p%func)) then
      errmsg = name // ' not defined for phase ' // this%name
      return
    end if
    call move_alloc(p, prop)

  end subroutine alloc_matl_prop

  subroutine compute_value(this, state, value)
    class(simple_prop), intent(in) :: this
    real(r8), intent(in) :: state(:)
    real(r8), intent(out) :: value
    value = this%func%eval(state)
  end subroutine compute_value

end module single_phase_matl_type
