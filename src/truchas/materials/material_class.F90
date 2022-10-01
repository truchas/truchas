!!
!! MATERIAL_CLASS
!!
!! This module provides the abstract base class MATERIAL that defines the
!! interface to properties and other info associated with a specific material.
!! It also provides a lower-level container PHASE for property functions and
!! attributes associated with a material or material phase.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! November 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module material_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use scalar_func_class
  use scalar_func_map_type
  use string_set_type
  implicit none
  private

  !! PHASE is the basic container for property functions and attributes that
  !! are associated with a material or a particular phase of a material. It
  !! anticipates the possibility of a hierarchical structure where a material
  !! may comprise a PHASE object with properties that span all phases of the
  !! material, and additional PHASE objects specific to each of its phases.
  !! If a property is not found in the container for a specific phase, it may
  !! be found in the container for the parent material of the phase. This
  !! chaining of PHASE containers is implemented via the MATL pointer component.

  type, public :: phase
    private
    character(:), allocatable, public :: name ! READ-ONLY
    type(scalar_func_map) :: prop_map
    type(string_set) :: attr_set
    type(phase), pointer, public :: matl => null() ! parent material as a phase
  contains
    procedure :: add_attr
    procedure :: has_attr
    procedure :: add_prop
    procedure :: has_prop
    procedure :: has_const_prop
    generic   :: get_prop => get_prop_func, get_prop_const
    procedure, private :: get_prop_func, get_prop_const
  end type

  !! MATERIAL is the primary abstract class which defines the interface used
  !! to access information about a material. It extends the PHASE type and
  !! thus is itself a container for properties and attributes that apply to
  !! a material as a whole. Concrete extensions of this class will implement
  !! different categories of a material; e.g. single-phase and multi-phase.

  type, extends(phase), abstract, public :: material
  contains
    procedure(alloc_matl_prop), deferred :: alloc_matl_prop
    procedure(num_phase), deferred :: num_phase
    procedure(phase_name), deferred :: phase_name
    procedure(phase_ref), deferred :: phase_ref
    procedure(get_phase_frac), deferred :: get_phase_frac
    procedure(add_enthalpy_prop), deferred :: add_enthalpy_prop
  end type

  abstract interface
    !! Allocate a MATL_PROP class object for a given property
    subroutine alloc_matl_prop(this, name, prop, errmsg)
      use matl_prop_class
      import material
      class(material), intent(in), target :: this ! possible internal ref to this subobject
      character(*), intent(in) :: name
      class(matl_prop), allocatable, intent(out) :: prop
      character(:), allocatable, intent(out) :: errmsg
    end subroutine

    !! Return the number of phases the material comprises
    integer function num_phase(this)
      import material
      class(material), intent(in) :: this
    end function

    !! Return the name of the nth phase of the material
    function phase_name(this, n)
      import material
      class(material), intent(in) :: this
      integer, intent(in) :: n
      character(:), allocatable :: phase_name
    end function

    !! Return a reference to the nth phase of the material
    function phase_ref(this, n) result(phi)
      import material, phase
      class(material), intent(in), target :: this
      integer, intent(in) :: n
      type(phase), pointer :: phi
    end function

    !! Return the phase fractions of the material at a given temperature
    subroutine get_phase_frac(this, temp, beta)
      import material, r8
      class(material), intent(in) :: this
      real(r8), intent(in) :: temp
      real(r8), intent(out) :: beta(:)
    end subroutine

    !! Build and add the 'enthalpy' property for the material
    !TODO: Is there any way to make this external to the class?
    subroutine add_enthalpy_prop(this, stat, errmsg)
      import material
      class(material), intent(inout) :: this
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
    end subroutine
  end interface

contains

  !! Adds the attribute NAME to the phase.
  subroutine add_attr(this, name)
    class(phase), intent(inout) :: this
    character(*), intent(in) :: name
    call this%attr_set%add(name)
  end subroutine add_attr

  !! Return true if the phase or its parent material has the attribute NAME;
  !! otherwise return false.
  recursive logical function has_attr(this, name)
    class(phase), intent(in) :: this
    character(*), intent(in) :: name
    has_attr = this%attr_set%has(name)
    if (has_attr) return  ! else check the parent material phase
    if (associated(this%matl)) has_attr = this%matl%has_attr(name)
  end function has_attr

  !! Add property NAME with SCALAR_FUNC class function FUNC to the phase.
  !! If the property already exists, its function is replaced by FUNC.
  !! The allocatable FUNC is moved into the object, not copied.
  subroutine add_prop(this, name, func)
    class(phase), intent(inout) :: this
    character(*), intent(in) :: name
    class(scalar_func), allocatable, intent(inout) :: func
    call this%prop_map%insert(name, func)
  end subroutine add_prop

  !! Return true if the phase or its parent material has property NAME;
  !! otherwise false. If the optional argument ONLY is specified, the query
  !! is limited to the case when the phase has the attribute ONLY; if not,
  !! true is always returned. If the optional argument STRICT is present
  !! with value true the parent material is excluded from the query.
  recursive logical function has_prop(this, name, only, strict) result(has)
    class(phase), intent(in) :: this
    character(*), intent(in) :: name
    character(*), intent(in), optional :: only
    logical, intent(in), optional :: strict
    if (present(only)) then
      has = has_prop(this, name, strict=strict) .or. .not.this%has_attr(only)
    else
      has = this%prop_map%mapped(name)
      if (has) return  ! else check the parent material phase
      if (present(strict)) then
        if (strict) return
      end if
      if (associated(this%matl)) has = this%matl%has_prop(name)
    end if
  end function has_prop

  !! Return true if the phase or its parent material has property NAME and
  !! it is constant; otherwise return false.
  recursive logical function has_const_prop(this, name)
    use scalar_func_tools, only: is_const
    class(phase), intent(in) :: this
    character(*), intent(in) :: name
    class(scalar_func), allocatable :: func
    if (this%prop_map%mapped(name)) then
      call this%prop_map%lookup(name, func)
      has_const_prop = is_const(func)
    else  ! check the parent material phase
      has_const_prop = .false.
      if (associated(this%matl)) has_const_prop = this%matl%has_const_prop(name)
    end if
  end function has_const_prop

  !! Return a polymorphic copy FUNC of the SCALAR_FUNC class function for
  !! phase property NAME. The property function associated with the phase
  !! takes precedence over that associated with the parent material. If
  !! neither the phase or its parent material has the property, FUNC is
  !! returned unallocated.
  recursive subroutine get_prop_func(this, name, func)
    class(phase), intent(in) :: this
    character(*), intent(in) :: name
    class(scalar_func), allocatable, intent(out) :: func
    call this%prop_map%lookup(name, func)
    if (allocated(func)) return ! else check the parent material phase
    if (associated(this%matl)) call this%matl%get_prop(name, func)
  end subroutine get_prop_func

  !! Return the value CONST of the constant phase property NAME. It is an
  !! error (unchecked) if the property is not constant.
  subroutine get_prop_const(this, name, const)
    class(phase), intent(in) :: this
    character(*), intent(in) :: name
    real(r8), intent(out) :: const
    class(scalar_func), allocatable :: func
    call get_prop_func(this, name, func)
    const = func%eval([real(r8)::])
  end subroutine get_prop_const

end module material_class
