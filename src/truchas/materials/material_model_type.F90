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
    private ! public components are READ ONLY
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
    procedure :: phase_matl_index
    procedure :: get_matl_phase_index_range
    procedure :: get_matl_phase_frac
    procedure :: get_phase_ref
    procedure :: get_matl_ref
    procedure :: has_const_phase_prop
    procedure :: const_phase_prop
    procedure :: get_phase_prop
    procedure :: get_matl_prop
  end type material_model

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
      this%is_fluid(p) = this%plist(p)%phi%has_attr('fluid')
    end do
    if (this%have_void) this%is_fluid(this%nphase) = .true.

    stat = 0

  end subroutine init

  !! Return the name of phase PID. PID must be in [1, NUM_PHASE] (unchecked)
  !! and will return VOID for the pseudo void phase, if it exists.
  function phase_name(this, pid) result(name)
    class(material_model), intent(in) :: this
    integer, intent(in) :: pid
    character(:), allocatable :: name
    if (pid == this%void_index) then
      name = VOID
    else
      name = this%plist(pid)%phi%name
    end if
  end function

  !! Return the index of the phase with name NAME. If no such phase exists,
  !! 0 is returned. NAME may be the reserved name 'VOID'.
  elemental integer function phase_index(this, name) result(n)
    class(material_model), intent(in) :: this
    character(*), intent(in) :: name
    do n = this%nphase_real, 1, -1
      if (this%plist(n)%phi%name == name) return
    end do
    if (this%have_void .and. name == VOID) n = this%void_index
  end function

  !! Return the index of the parent material for the phase with index PID.
  !! PID must be in [1, NUM_PHASE].
  elemental integer function phase_matl_index(this, pid) result(mid)
    class(material_model), intent(in) :: this
    integer, intent(in) :: pid
    mid = this%p2m(pid)
  end function

  !! Return true if a phase with name NAME exists; otherwise false.
  !! NAME may be the reserved name 'VOID'.
  elemental logical function has_phase(this, name)
    class(material_model), intent(in) :: this
    character(*), intent(in) :: name
    has_phase = (phase_index(this, name) > 0)
  end function

  !! Return the name of material MID. MID must be in [1, NUM_MATL] (unchecked)
  !! and will return VOID for the pseudo void material, if it exists.
  function matl_name(this, mid) result(name)
    class(material_model), intent(in) :: this
    integer, intent(in) :: mid
    character(:), allocatable :: name
    if (this%have_void .and. mid == this%nmatl) then
      name = VOID
    else
      name = this%mlist(mid)%matl%name
    end if
  end function

  !! Return the index of the material with name NAME. if no such material
  !! exists, 0 is returned. NAME may be the reserved name 'VOID'.
  elemental integer function matl_index(this, name) result(mid)
    class(material_model), intent(in) :: this
    character(*), intent(in) :: name
    do mid = this%nmatl_real, 1, -1
      if (this%mlist(mid)%matl%name == name) return
    end do
    if (this%have_void .and. name == VOID) mid = this%nmatl
  end function

  !! Return true if a material with name NAME exists; otherwise false.
  !! NAME may be the reserved name 'VOID'.
  elemental logical function has_matl(this, name)
    class(material_model), intent(in) :: this
    character(*), intent(in) :: name
    has_matl = (matl_index(this, name) > 0)
  end function

  !! Return the number of phases that comprise the material with index MID.
  !! MID must be in [1, NUM_REAL_MATL) (unchecked).
  integer function num_matl_phase(this, mid)
    class(material_model), intent(in) :: this
    integer, intent(in) :: mid
    num_matl_phase = this%mlist(mid)%matl%num_phase()
  end function

  !! Phases are grouped by the material to which they belong. This subroutine
  !! returns the phase index range [FIRST, LAST] of the phases that comprise
  !! material MID. Currently the phases are ordered from low to high-temperature.
  !! MID must be in [1, NUM_MATL] and may be the VOID material index.
  subroutine get_matl_phase_index_range(this, mid, first, last)
    class(material_model), intent(in) :: this
    integer, intent(in) :: mid
    integer, intent(out) :: first, last
    first = this%m2p(mid)
    last  = this%m2p(mid+1) - 1
  end subroutine

  !! Return the phase fractions BETA for material MID at the temperature TEMP.
  !! The order of values in the BETA array correspond to the order of the phases
  !! of the material. The size of BETA may be larger than the number of phases.
  !! Unused elements at the end of the array are left unchanged. MID must be in
  !! [1, NUM_REAL_MATL] (unchecked).
  subroutine get_matl_phase_frac(this, mid, temp, beta)
    class(material_model), intent(in) :: this
    integer,  intent(in)  :: mid
    real(r8), intent(in)  :: temp
    real(r8), intent(inout) :: beta(:)  ! INOUT to not alter unused elements
    call this%mlist(mid)%matl%get_phase_frac(temp, beta)
  end subroutine

  !! Return a MATERIAL class pointer to material MID. MID must be in
  !! [1, NUM_REAL_MATL] (unchecked).
  subroutine get_matl_ref(this, mid, matl)
    class(material_model), intent(in) :: this
    integer, intent(in) :: mid
    class(material), pointer :: matl
    matl => this%mlist(mid)%matl
  end subroutine

  !! Return a PHASE class pointer to phase PID. PID must be in [1, NUM_REAL_PHASE]
  !! (unchecked).
  !TODO: should this be a TYPE(PHASE) pointer?
  subroutine get_phase_ref(this, pid, phi)
    class(material_model), intent(in) :: this
    integer, intent(in) :: pid
    type(phase), pointer :: phi
    phi => this%plist(pid)%phi
  end subroutine

  !! Return true if phase PID has the constant property NAME. PID must be in
  !! [1, NUM_REAL_PHASE] (unchecked).
  logical function has_const_phase_prop(this, pid, name)
    class(material_model), intent(in) :: this
    integer, intent(in) :: pid
    character(*), intent(in) :: name
    has_const_phase_prop = this%plist(pid)%phi%has_const_prop(name)
  end function

  !! Return the value of the constant property NAME of phase PID. PID must be in
  !! [1, NUM_REAL_PHASE] (unchecked). The property must exist and be constant
  !! valued (unchecked).
  function const_phase_prop(this, pid, name) result(const)
    class(material_model), intent(in) :: this
    integer, intent(in) :: pid
    character(*), intent(in) :: name
    real(r8) :: const
    call this%plist(pid)%phi%get_prop(name, const)
  end function

  !! Allocate a SCALAR_FUNC class copy PROP of property NAME of phase PID.
  !! PID must be in [1, NUM_REAL_PHASE] (unchecked). If the property does not
  !! exist, PROP is returned unallocated.
  subroutine get_phase_prop(this, pid, name, prop)
    use scalar_func_class
    class(material_model), intent(in) :: this
    integer, intent(in) :: pid
    character(*), intent(in) :: name
    class(scalar_func), allocatable, intent(out) :: prop
    call this%plist(pid)%phi%get_prop(name, prop)
  end subroutine

  !! Allocate a MATL_PROP class object PROP for property NAME of material MID.
  !! MID must be in [1, NUM_REAL_MATL] (unchecked). If an error occurs, PROP is
  !! returned unallocated and an explanatory message is returned in ERRMSG.
  subroutine get_matl_prop(this, mid, name, prop, errmsg)
    class(material_model), intent(in) :: this
    integer, intent(in) :: mid
    character(*), intent(in) :: name
    class(matl_prop), allocatable, intent(out) :: prop
    character(:), allocatable, intent(out) :: errmsg
    call this%mlist(mid)%matl%alloc_matl_prop(name, prop, errmsg)
  end subroutine

end module material_model_type
