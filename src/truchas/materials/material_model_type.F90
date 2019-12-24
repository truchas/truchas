!!
!! MATERIAL_MODEL_TYPE
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

module material_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use material_class
  use matl_prop_class
  implicit none
  private

  type :: material_box
    class(material), pointer :: matl
  end type

  type :: phase_box
    type(phase), pointer :: phi
  end type

  type, public :: material_model
    private
    !TODO: *read-only* public components, or function getters?
    integer, public :: nmatl       ! number of materials, including void
    integer, public :: nmatl_real  ! number of materials, excluding void
    integer, public :: nphase      ! total number of phases, including void
    integer, public :: nphase_real ! total number of phases, excluding void
    logical, public :: have_void
    integer, public :: void_index     ! either 0 or NUM_PHASE
    logical, allocatable, public :: is_fluid(:) ! read-only
    type(phase_box), allocatable :: plist(:)
    type(material_box), allocatable :: mlist(:)
    integer, allocatable :: p2m(:)  ! phase index to parent material index
    integer, allocatable :: m2p(:)  ! material index to first material phase index
  contains
    procedure :: init
    procedure :: phase_name
    procedure :: phase_index
    procedure :: has_phase
    procedure :: matl_name
    procedure :: matl_index
    procedure :: has_matl
    procedure :: num_matl_phase
    procedure :: get_matl_phase_index_range
    procedure :: get_matl_phase_frac
    procedure :: get_phase_ref
    procedure :: get_matl_ref
    procedure :: const_phase_prop
    procedure :: alloc_phase_prop
    procedure :: add_phase_prop
    procedure :: alloc_matl_prop
    generic   :: alloc_avg_matl_prop => alloc_avg_matl_prop_list, alloc_avg_matl_prop_all
    procedure, private :: alloc_avg_matl_prop_list
    procedure, private :: alloc_avg_matl_prop_all
    procedure :: alloc_avg_phase_prop => alloc_avg_phase_prop_list
    !TODO: really would like these to be external to this class
    procedure :: required_property_check
    procedure :: required_fluid_property_check
    procedure :: define_property_default
    procedure :: define_fluid_property_default
    procedure :: optional_property_check
    procedure :: constant_property_check
    procedure :: create_enthalpy
  end type

  character(4), parameter :: VOID = 'VOID'

contains


  subroutine init(this, matl_names, matl_db, stat, errmsg)

    use material_database_type

    class(material_model), intent(out) :: this
    character(*), intent(in) :: matl_names(:)
    type(material_database), intent(in), target :: matl_db
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: j, k, m, p
    class(material), pointer :: matl

    !! Check for any 0-length or repeated material names.
    do j = 1, size(matl_names)
      if (matl_names(j) == '') then
        stat = 1
        errmsg = 'invalid 0-length material name'
        return
      end if
      do k = j+1, size(matl_names)
        if (matl_names(k) == matl_names(j)) then
          stat = 1
          errmsg = 'repeated material name: ' // trim(matl_names(j))
          return
        end if
      end do
    end do

    !! Verify that the given materials are defined in the database.
    do j = 1, size(matl_names)
      if (matl_names(j) == VOID) cycle
      if (matl_db%has_matl(matl_names(j))) cycle
      stat = 1
      errmsg = 'material ' // trim(matl_names(j)) // ' not defined'
      return
    end do

    !! Number of materials, possibly including VOID.
    this%nmatl = size(matl_names)

    !! Number of phases, possibly including VOID.
    this%nphase = 0
    do j = 1, size(matl_names)
      if (matl_names(j) == VOID) then
        this%nphase = this%nphase + 1
      else
        matl => matl_db%matl_ref(matl_names(j))
        this%nphase = this%nphase + matl%num_phase()
      end if
    end do

    !! Account for the possible presence of VOID.
    this%have_void = any(matl_names == VOID)
    this%nmatl_real  = this%nmatl  - merge(1,0,this%have_void)
    this%nphase_real = this%nphase - merge(1,0,this%have_void)
    this%void_index = merge(this%nphase, 0, this%have_void) !TODO: phase index; matl?

    !! Create the arrays of MATERIAL and PHASE object references. Phases are
    !! grouped by their parent material in order. At the same time we generate
    !! the mapping P2M from phase index to parent material index, and the
    !! mapping M2P from material index to the index of its first phase. If VOID
    !! is listed it is skipped and will be associated with the last material/
    !! phase index, and will have no corresponding MATERIAL or PHASE object
    !! references.
    allocate(this%plist(this%nphase_real), this%mlist(this%nmatl_real))
    allocate(this%p2m(this%nphase), this%m2p(this%nmatl+1))
    m = 0; p = 0
    do j = 1, size(matl_names)
      if (matl_names(j) == VOID) cycle
      m = m + 1
      this%m2p(m) = p + 1
      this%mlist(m)%matl => matl_db%matl_ref(matl_names(j))
      associate (matl => this%mlist(m)%matl)
        do k = 1, matl%num_phase()
          p = p + 1
          this%p2m(p) = m
          this%plist(p)%phi => matl%phase_ref(k)
        end do
      end associate
    end do
    this%m2p(this%nmatl+1) = this%nphase + 1
    INSIST(m == this%nmatl_real)
    INSIST(p == this%nphase_real)

    !! Extend the M2P and P2M mapping to account for VOID if present.
    if (this%have_void) then
      this%p2m(this%nphase) = this%nmatl
      this%m2p(this%nmatl) = this%nphase
    end if

    !! Define the IS_FLUID phase mask as a convenience. VOID is regarded
    !! as a fluid.
    allocate(this%is_fluid(this%nphase))
    do p = 1, this%nphase_real
      this%is_fluid(p) = this%plist(p)%phi%has_attr('is-fluid')
    end do
    if (this%have_void) this%is_fluid(this%nphase) = .true.

    stat = 0

  end subroutine init

  !! Return the name of phase N. N must be in [1, NUM_PHASE] (unchecked)
  !! and will return VOID for the pseudo void phase, if it exists.
  function phase_name(this, n) result(name)
    class(material_model), intent(in) :: this
    integer, intent(in) :: n
    character(:), allocatable :: name
    if (n == this%void_index) then
      name = VOID
    else
      name = this%plist(n)%phi%name
    end if
  end function phase_name

  !! Return the index of the phase with name NAME. If no such phase exists,
  !! 0 is returned. NAME may be the reserved name 'VOID'.
  elemental integer function phase_index(this, name) result(n)
    class(material_model), intent(in) :: this
    character(*), intent(in) :: name
    do n = this%nphase_real, 1, -1
      if (this%plist(n)%phi%name == name) return
    end do
    if (this%have_void .and. name == VOID) n = this%void_index
  end function phase_index

  !! Return true if a phase with name NAME exists; otherwise false.
  !! NAME may be the reserved name 'VOID'.
  elemental logical function has_phase(this, name)
    class(material_model), intent(in) :: this
    character(*), intent(in) :: name
    has_phase = (phase_index(this, name) > 0)
  end function has_phase

  !! Return the name of material N. N must be in [1, NUM_MATL] (unchecked)
  !! and will return VOID for the pseudo void material, if it exists.
  function matl_name(this, n) result(name)
    class(material_model), intent(in) :: this
    integer, intent(in) :: n
    character(:), allocatable :: name
    if (this%have_void .and. n == this%nmatl) then
      name = VOID
    else
      name = this%mlist(n)%matl%name
    end if
  end function matl_name

  !! Return the index of the material with name NAME. if no such material
  !! exists, 0 is returned. NAME may be the reserved name 'VOID'.
  elemental integer function matl_index(this, name) result(n)
    class(material_model), intent(in) :: this
    character(*), intent(in) :: name
    do n = this%nmatl_real, 1, -1
      if (this%mlist(n)%matl%name == name) return
    end do
    if (this%have_void .and. name == VOID) n = this%nmatl
  end function matl_index

  !! Return true if a material with name NAME exists; otherwise false.
  !! NAME may be the reserved name 'VOID'.
  elemental logical function has_matl(this, name)
    class(material_model), intent(in) :: this
    character(*), intent(in) :: name
    has_matl = (matl_index(this, name) > 0)
  end function has_matl

  !! Return the number of phases that comprise the Nth material.
  !! N must be in [1, NUM_REAL_MATL) (unchecked).
  integer function num_matl_phase(this, n)
    class(material_model), intent(in) :: this
    integer, intent(in) :: n
    num_matl_phase = this%mlist(n)%matl%num_phase()
  end function

  !! Phases are grouped by the material to which they belong. This subroutine
  !! returns the phase index range [FIRST, LAST] of the phases that comprise
  !! material N. Currently the phases are ordered from low to high-temperature.
  !! N must be in [1, NUM_MATL] and may be the VOID material index.
  subroutine get_matl_phase_index_range(this, n, first, last)
    class(material_model), intent(in) :: this
    integer, intent(in) :: n
    integer, intent(out) :: first, last
    first = this%m2p(n)
    last  = this%m2p(n+1) - 1
  end subroutine get_matl_phase_index_range

  !! Return the phase fractions BETA for material N at the temperature TEMP.
  !! The order of values in the BETA array correspond to the order of the phases
  !! of the material. The size of BETA may be larger than the number of phases.
  !! Unused elements at the end of the array are left unchanged. N must be in
  !! [1, NUM_REAL_MATL] (unchecked).
  subroutine get_matl_phase_frac(this, n, temp, beta)
    class(material_model), intent(in) :: this
    integer,  intent(in)  :: n
    real(r8), intent(in)  :: temp
    real(r8), intent(inout) :: beta(:)  ! INOUT to not alter unused elements
    call this%mlist(n)%matl%get_phase_frac(temp, beta)
  end subroutine

  !! Return a MATERIAL class pointer to material N. N must be in
  !! [1, NUM_REAL_MATL] (unchecked).
  subroutine get_matl_ref(this, n, matl)
    class(material_model), intent(in), target :: this
    integer, intent(in) :: n
    class(material), pointer :: matl
    matl => this%mlist(n)%matl
  end subroutine

  !! Return a PHASE class pointer to phase N. N must be in [1, NUM_REAL_PHASE]
  !! (unchecked).
  !TODO: should this be a TYPE(PHASE) pointer?
  subroutine get_phase_ref(this, n, phi)
    class(material_model), intent(in), target :: this
    integer, intent(in) :: n
    class(phase), pointer :: phi
    phi => this%plist(n)%phi
  end subroutine

  !! Return the value of the constant property NAME of phase N. N must be in
  !! [1, NUM_REAL_PHASE] (unchecked). The property must exist and be constant
  !! valued (unchecked).
  function const_phase_prop(this, n, name) result(const)
    class(material_model), intent(in) :: this
    integer, intent(in) :: n
    character(*), intent(in) :: name
    real(r8) :: const
    call this%plist(n)%phi%get_prop(name, const)
  end function const_phase_prop

  !! Allocate a SCALAR_FUNC class copy FUNC of property NAME of phase N.
  !! N must be in [1, NUM_REAL_PHASE] (unchecked). If the property does not
  !! exist, FUNC is returned unallocated.
  subroutine alloc_phase_prop(this, n, name, func)
    use scalar_func_class
    class(material_model), intent(in) :: this
    integer, intent(in) :: n
    character(*), intent(in) :: name
    class(scalar_func), allocatable, intent(out) :: func
    call this%plist(n)%phi%get_prop(name, func)
  end subroutine alloc_phase_prop

  !! Adds property NAME with SCALAR_FUNC class function FUNC to phase N.
  !! N must be in [1, NUM_REAL_PHASE] (unchecked). FUNC is moved into the
  !! phase object, not copied.
  subroutine add_phase_prop(this, n, name, func)
    use scalar_func_class
    class(material_model), intent(in) :: this
    integer, intent(in) :: n
    character(*), intent(in) :: name
    class(scalar_func), allocatable, intent(inout) :: func
    call this%plist(n)%phi%add_prop(name, func)
  end subroutine

  !! Allocate a MATL_PROP class object PROP for property NAME of material N.
  !! N must be in [1, NUM_REAL_MATL] (unchecked). If an error occurs, PROP is
  !! returned unallocated and an explanatory message is returned in ERRMSG.
  subroutine alloc_matl_prop(this, n, name, prop, errmsg)
    class(material_model), intent(in) :: this
    integer, intent(in) :: n
    character(*), intent(in) :: name
    class(matl_prop), allocatable, intent(out) :: prop
    character(:), allocatable, intent(out) :: errmsg
    call this%mlist(n)%matl%alloc_matl_prop(name, prop, errmsg)
  end subroutine alloc_matl_prop

  !! Allocate an AVG_MATL_PROP type object FUNC for property NAME. MIDS
  !! is a list of material indexes to include in the object. All must be in
  !! the range [1, NUM_REAL_MATL]. Later evaluation of FUNC requires an array
  !! of material weights corresponding to MIDS. If an error occurs, FUNC is
  !! returned unallocated and an explanatory message is returned in ERRMSG.
  subroutine alloc_avg_matl_prop_list(this, name, mids, func, errmsg)
    use avg_matl_prop_type
    class(material_model), intent(in) :: this
    character(*), intent(in) :: name
    integer, intent(in) :: mids(:)
    type(avg_matl_prop), allocatable, intent(out) :: func
    character(:), allocatable, intent(out) :: errmsg
    integer :: n
    allocate(func)
    allocate(func%matl(size(mids)))
    do n = 1, size(mids)
      associate (matl => this%mlist(mids(n))%matl)
        call matl%alloc_matl_prop(name, func%matl(n)%prop, errmsg)
        if (.not.allocated(func%matl(n)%prop)) return
      end associate
    end do
  end subroutine alloc_avg_matl_prop_list

  !! Allocate an AVG_MATL_PROP type object FUNC for property NAME. All real
  !! materials are included in the object. Later evaluation of FUNC requires
  !! an array of material weights corresponding to the real materials in order.
  !! If an error occurs, FUNC is returned unallocated and an explanatory message
  !! is returned in ERRMSG.
  subroutine alloc_avg_matl_prop_all(this, name, func, errmsg)
    use avg_matl_prop_type
    class(material_model), intent(in) :: this
    character(*), intent(in) :: name
    type(avg_matl_prop), allocatable, intent(out) :: func
    character(:), allocatable, intent(out) :: errmsg
    integer :: n, all_mids(this%nmatl_real)
    all_mids = [(n, n=1,this%nmatl_real)]
    call alloc_avg_matl_prop_list(this, name, all_mids, func, errmsg)
  end subroutine alloc_avg_matl_prop_all

  !! Allocate an AVG_MATL_PROP type object FUNC for property NAME. PIDS is a
  !! list of phase indexes to include in the object. All must be in the range
  !! [1, NUM_REAL_PHASE]. Later evaluation of FUNC requires an array of phase
  !! weights corresponding to PIDS. If an error occurs, FUNC is returned
  !! unallocated and an explanatory message is returned in ERRMSG.

  subroutine alloc_avg_phase_prop_list(this, name, pids, func, errmsg)
    use avg_phase_prop_type
    class(material_model), intent(in) :: this
    character(*), intent(in) :: name
    integer, intent(in) :: pids(:)
    type(avg_phase_prop), allocatable, intent(out) :: func
    character(:), allocatable, intent(out) :: errmsg
    integer :: n
    allocate(func)
    allocate(func%phase(size(pids)))
    do n = 1, size(pids)
      associate (phi => this%plist(pids(n))%phi)
        call phi%get_prop(name, func%phase(n)%func)
        if (.not.allocated(func%phase(n)%func)) then
          errmsg = name // ' undefined for phase ' // phi%name
          return
        end if
      end associate
    end do
  end subroutine alloc_avg_phase_prop_list

  !! Check that property NAME is defined for all materials. If not, STAT
  !! returns a nonzero value and and explanatory message is returned in ERRMSG.

  subroutine required_property_check(this, name, stat, errmsg)
    class(material_model), intent(in) :: this
    character(*), intent(in) :: name
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call req_prop_check(this, name, stat, errmsg)
  end subroutine required_property_check

  !! Check that property NAME is defined for all fluid material phases. If not,
  !! STAT returns a nonzero value and and explanatory message is returned in
  !! ERRMSG.

  subroutine required_fluid_property_check(this, name, stat, errmsg)
    class(material_model), intent(in) :: this
    character(*), intent(in) :: name
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call req_prop_check(this, name, stat, errmsg, only='is-fluid')
    if (stat /= 0) errmsg = 'fluid ' // errmsg
  end subroutine required_fluid_property_check

  subroutine req_prop_check(this, name, stat, errmsg, only)
    class(material_model), intent(in) :: this
    character(*), intent(in) :: name
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    character(*), intent(in), optional :: only
    integer :: n
    do n = 1, size(this%mlist)
      associate (matl => this%mlist(n)%matl)
        if (matl%has_prop(name, only)) cycle
        stat = 1
        errmsg = name // ' undefined for material ' // matl%name
        return
      end associate
    end do
    stat = 0
  end subroutine req_prop_check

  !! Ensure that property NAME is defined for all phases by assigning the
  !! constant value DEFAULT to those phases without the property.

  subroutine define_property_default(this, name, default)
    class(material_model), intent(in) :: this
    character(*), intent(in) :: name
    real(r8), intent(in) :: default
    call def_prop_default(this, name, default)
  end subroutine define_property_default

  !! Ensure that property NAME is defined for all fluid phases by assigning
  !! the constant value DEFAULT to those fluid phases without the property.
  subroutine define_fluid_property_default(this, name, default)
    class(material_model), intent(in) :: this
    character(*), intent(in) :: name
    real(r8), intent(in) :: default
    call def_prop_default(this, name, default, only='is-fluid')
  end subroutine define_fluid_property_default

  subroutine def_prop_default(this, name, default, only)
    use scalar_func_factories
    class(material_model), intent(in) :: this
    character(*), intent(in) :: name
    real(r8), intent(in) :: default
    character(*), intent(in), optional :: only
    integer :: n
    class(scalar_func), allocatable :: f_default
    do n = 1, size(this%plist)
      associate (phi => this%plist(n)%phi)
        if (phi%has_prop(name, only)) cycle
        call alloc_const_scalar_func(f_default, default)
        call phi%add_prop(name, f_default)
      end associate
    end do
  end subroutine def_prop_default

  !TODO: Rethink. This replicates the original procedure and seems clumsy.
  !! Examine the material_model for property NAME. If all phases have the
  !! property, STAT returns 0. If none of the phases have the property
  !! STAT returns 1. Otherwise some phases have the property while others
  !! do not, and STAT returns -1 and an explanatory error message is
  !! returned in ERRMSG.

  subroutine optional_property_check(this, name, stat, errmsg)

    class(material_model), intent(inout) :: this
    character(*), intent(in) :: name
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    logical :: some_have, some_lack

    some_have = .false.; some_lack = .false.
    do n = 1, size(this%plist)
      associate (phi => this%plist(n)%phi)
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
      end associate
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

  subroutine constant_property_check(this, name, stat, errmsg)
    class(material_model), intent(in) :: this
    character(*), intent(in) :: name
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer :: n
    do n = 1, size(this%mlist)
      associate (matl => this%mlist(n)%matl)
        if (matl%has_const_prop(name)) cycle
        if (allocated(errmsg)) then
          errmsg = errmsg // ', ' // matl%name
        else
          errmsg = matl%name
        end if
      end associate
    end do
    stat = merge(1, 0, allocated(errmsg))
  end subroutine constant_property_check

  !! Create an enthalpy function for each material from density, specific heat,
  !! and latent heats. If an error occurs, STAT returns a nonzero value and
  !! and explanatory message is returned in ERRMSG.

  subroutine create_enthalpy(this, stat, errmsg)

    class(material_model), intent(inout) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n

    do n = 1, size(this%mlist)
      associate (matl => this%mlist(n)%matl)
        call matl%create_enthalpy(stat, errmsg)
        if (stat /= 0) then
          errmsg = 'cannot define enthalpy density for material ' // matl%name // ': ' // errmsg
          return
        end if
      end associate
    end do

  end subroutine create_enthalpy

end module material_model_type
