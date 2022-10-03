!!
!! MATERIAL_UTILITIES
!!
!! This provides a collection of specialized utilities for
!! Neil N. Carlson <nnc@lanl.gov>
!! November 2019
!!

module material_utilities

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use material_class
  use material_model_type
  implicit none
  private

  public :: required_property_check, required_fluid_property_check
  public :: optional_property_check, constant_property_check
  public :: define_property_default, define_fluid_property_default
  public :: add_enthalpy_prop

contains

  !! Check that property NAME is defined for all materials. If not, STAT
  !! returns a nonzero value and an error message is returned in ERRMSG.

  subroutine required_property_check(model, name, stat, errmsg)
    class(material_model), intent(in) :: model
    character(*), intent(in) :: name
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call req_prop_check_aux(model, name, stat, errmsg)
  end subroutine required_property_check

  !! Check that property NAME is defined for all fluid material phases. If not,
  !! STAT returns a nonzero value and an error message is returned in ERRMSG.

  subroutine required_fluid_property_check(model, name, stat, errmsg)
    class(material_model), intent(in) :: model
    character(*), intent(in) :: name
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call req_prop_check_aux(model, name, stat, errmsg, only='fluid')
    if (stat /= 0) errmsg = 'fluid ' // errmsg
  end subroutine required_fluid_property_check

  subroutine req_prop_check_aux(model, name, stat, errmsg, only)
    class(material_model), intent(in) :: model
    character(*), intent(in) :: name
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    character(*), intent(in), optional :: only
    integer :: n
    class(material), pointer :: matl
    do n = 1, model%nmatl_real
      call model%get_matl_ref(n, matl)
      if (matl%has_prop(name, only)) cycle
      stat = 1
      errmsg = name // ' undefined for material ' // matl%name
      return
    end do
    stat = 0
  end subroutine req_prop_check_aux

  !! Ensure that property NAME is defined for all phases by assigning the
  !! constant value DEFAULT to those phases without the property.

  subroutine define_property_default(model, name, default)
    class(material_model), intent(in) :: model
    character(*), intent(in) :: name
    real(r8), intent(in) :: default
    call def_prop_default(model, name, default)
  end subroutine define_property_default

  !! Ensure that property NAME is defined for all fluid phases by assigning
  !! the constant value DEFAULT to those fluid phases without the property.
  subroutine define_fluid_property_default(model, name, default)
    class(material_model), intent(in) :: model
    character(*), intent(in) :: name
    real(r8), intent(in) :: default
    call def_prop_default(model, name, default, only='fluid')
  end subroutine define_fluid_property_default

  subroutine def_prop_default(model, name, default, only)
    use scalar_func_factories
    class(material_model), intent(in) :: model
    character(*), intent(in) :: name
    real(r8), intent(in) :: default
    character(*), intent(in), optional :: only
    integer :: n
    type(phase), pointer :: phi
    class(scalar_func), allocatable :: f_default
    do n = 1, model%nphase_real
      call model%get_phase_ref(n, phi)
      if (phi%has_prop(name, only)) cycle
      call alloc_const_scalar_func(f_default, default)
      call phi%add_prop(name, f_default)
    end do
  end subroutine def_prop_default

  subroutine add_enthalpy_prop(model, stat, errmsg)

    class(material_model), intent(inout) :: model
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    class(material), pointer :: matl

    do n = 1, model%nmatl_real
      call model%get_matl_ref(n, matl)
      call matl%add_enthalpy_prop(stat, errmsg)
      if (stat /= 0) then
        errmsg = 'cannot define enthalpy density for material ' // matl%name // ': ' // errmsg
        return
      end if
    end do

  end subroutine add_enthalpy_prop

  !TODO: Rethink. This replicates the original procedure and seems clumsy.
  !! Examine the material_model for property NAME. If all phases have the
  !! property, STAT returns 0. If none of the phases have the property
  !! STAT returns 1. Otherwise some phases have the property while others
  !! do not, and STAT returns -1 and an explanatory error message is
  !! returned in ERRMSG.

  subroutine optional_property_check(model, name, stat, errmsg)

    class(material_model), intent(inout) :: model
    character(*), intent(in) :: name
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    logical :: some_have, some_lack
    type(phase), pointer :: phi

    some_have = .false.; some_lack = .false.
    do n = 1, model%nphase_real
      call model%get_phase_ref(n, phi)
      if (phi%has_prop(name)) then
        some_have = .true.
      else
        some_lack = .true.
      end if
      if (some_have .and. some_lack) then
        stat = -1
        errmsg = 'incomplete specification of property ' // name
        return
      end if
    end do

    if (some_have) then ! all phases have the property
      stat = 0
    else  ! none have the property
      stat = 1
    end if

  end subroutine optional_property_check

  !! Check that each real material has the constant property NAME. The value
  !! of the constant may differ between materials. If the check fails, STAT
  !! returns a nonzero value and an explanatory message is returned in ERRMSG.

  subroutine constant_property_check(model, name, stat, errmsg)
    class(material_model), intent(in) :: model
    character(*), intent(in) :: name
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer :: n
    class(material), pointer :: matl
    do n = 1, model%nmatl_real
      call model%get_matl_ref(n, matl)
      if (matl%has_const_prop(name)) cycle
      if (allocated(errmsg)) then
        errmsg = errmsg // ', ' // matl%name
      else
        errmsg = matl%name
      end if
    end do
    stat = merge(1, 0, allocated(errmsg))
  end subroutine constant_property_check

end module material_utilities
