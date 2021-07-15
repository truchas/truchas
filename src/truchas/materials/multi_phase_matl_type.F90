!! MULTI_PHASE_MATL_TYPE
!!
!! This module provides the concrete implementation MULTI_PHASE_MATL of
!! the MATERIAL base class that defines a multi-phase material.
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

module multi_phase_matl_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use material_class
  use phase_change_class
  use matl_prop_class
  use scalar_func_class
  implicit none
  private

  type :: phase_change_box
    class(phase_change), allocatable :: pc
  end type

  type, extends(material), public :: multi_phase_matl
    !private
    type(phase), allocatable :: phi(:)
    type(phase_change_box), allocatable :: pc_seq(:)
  contains
    ! Overloaded procedures from PHASE
    procedure :: has_attr
    procedure :: has_prop
    procedure :: has_const_prop
    ! Deferred procedures from MATERIAL
    procedure :: alloc_matl_prop
    procedure :: num_phase
    procedure :: phase_name
    procedure :: phase_ref
    procedure :: get_phase_frac
    procedure :: add_enthalpy_prop
    procedure :: write_solid_frac_plotfile
  end type

  type :: scalar_func_box
    class(scalar_func), allocatable :: func
  end type

  type, extends(matl_prop) :: multiphase_prop
    type(multi_phase_matl), pointer :: matl => null()  ! reference only -- not owned
    type(scalar_func_box), allocatable :: farray(:)
  contains
    ! Deferred procedure from MATL_PROP
    procedure :: compute_value
  end type

contains

  !! Overload the base class HAS_ATTR procedure. Returns true if every phase of
  !! the material has attribue NAME.
  logical function has_attr(this, name)
    class(multi_phase_matl), intent(in) :: this
    character(*), intent(in) :: name
    has_attr = all(has_prop_mask(this, name))
  end function

  !! Overload the base class HAS_PROP procedure. Returns true if every phase of
  !! the material has property NAME. If the optional argument ONLY is present
  !! the query is limited to phases having the attribute ONLY; if no phase has
  !! that attribute, true is returned. The optional argument STRICT is currently
  !! ignored for this overloading procedure.
  !! NB: If STRICT were passed through HAS_PROP_MASK to PHASE%HAS_PROP, the
  !! behavior would be to ignore a property defined at the material level and
  !! restrict the query to the property defined per material phase. Useful?
  logical function has_prop(this, name, only, strict)
    class(multi_phase_matl), intent(in) :: this
    character(*), intent(in) :: name
    character(*), intent(in), optional :: only
    logical, intent(in), optional :: strict ! ignored
    if (present(only)) then
      has_prop = all(has_prop_mask(this, name) .or. .not.has_attr_mask(this, only))
    else
      has_prop = all(has_prop_mask(this, name))
    end if
  end function has_prop

  function has_prop_mask(this, name) result(mask)
    class(multi_phase_matl), intent(in) :: this
    character(*), intent(in) :: name
    logical :: mask(size(this%phi))
    integer :: n
    do n = 1, size(mask)
      mask(n) = this%phi(n)%has_prop(name)
    end do
  end function has_prop_mask

  function has_attr_mask(this, name) result(mask)
    class(multi_phase_matl), intent(in) :: this
    character(*), intent(in) :: name
    logical :: mask(size(this%phi))
    integer :: n
    do n = 1, size(mask)
      mask(n) = this%phi(n)%has_attr(name)
    end do
  end function has_attr_mask

  !! Overload the base class HAS_CONST_PROP procedure. Return true if every
  !! phase of the material has the constant property NAME and its value is
  !! the same for all phases.
  logical function has_const_prop(this, name)
    class(multi_phase_matl), intent(in) :: this
    character(*), intent(in) :: name
    integer :: n
    real(r8) :: c, cref
    has_const_prop = .false.
    do n = 1, size(this%phi)
      if (.not.this%phi(n)%has_const_prop(name)) return
      call this%phi(n)%get_prop(name, c)
      if (n == 1) then
        cref = c
      else
        if (c /= cref) return
      end if
    end do
    has_const_prop = .true.
  end function has_const_prop

  integer function num_phase(this)
    class(multi_phase_matl), intent(in) :: this
    num_phase = size(this%phi)
  end function

  function phase_name(this, n) result(name)
    class(multi_phase_matl), intent(in) :: this
    integer, intent(in) :: n
    character(:), allocatable :: name
    name = this%phi(n)%name
    ! name = this%name // ':' // this%phi(n)%name !TODO?
  end function

  function phase_ref(this, n) result(phi)
    class(multi_phase_matl), intent(in), target :: this
    integer, intent(in) :: n
    type(phase), pointer :: phi
    phi => this%phi(n)
  end function

  subroutine get_phase_frac(this, temp, beta)

    class(multi_phase_matl), intent(in) :: this
    real(r8), intent(in) :: temp
    real(r8), intent(out) :: beta(:)

    integer :: n

    ASSERT(size(beta) == size(this%pc_seq)+1)

    beta = 0
    do n = 1, size(this%pc_seq)
      associate (pc => this%pc_seq(n)%pc)
        if (temp <= pc%solidus_temp()) then
          beta(n) = 1
          return
        else if (temp < pc%liquidus_temp()) then
          beta(n) = pc%solid_frac(temp)
          beta(n+1) = 1 - beta(n)
          return
        end if
      end associate
    end do
    beta(n) = 1

  end subroutine get_phase_frac

  !! If the material does not already have the specific enthalpy property,
  !! this subroutine attempts to build it. For each phase that is missing
  !! it, the property is built from the phase's specific heat property by
  !! integration, if possible. This requires the latent heat with respect
  !! to the adjacent lower temperature phase. STAT returns a non-zero value
  !! if the enthalpy property is not ultimately defined, and ERRMSG returns
  !! an explanatory message.
  !! TODO: Is there some way to make this external to the class?

  subroutine add_specific_enthalpy_prop(this, stat, errmsg)

    use scalar_func_tools

    class(multi_phase_matl), intent(inout) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    class(scalar_func), allocatable :: cp, h, hd
    real(r8) :: t0, h0
    integer :: n

    if (this%has_prop('specific-enthalpy')) return

    !! Set the initial reference temperature and specific enthalpy.
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

    !! While it is possible to have a single material-wide specific heat, it
    !! would not be understood to incorporate the latent heat of phase change.
    !! As such, the integrated specific enthalpy would still be defined per phase.

    do n = 1, size(this%phi)

      if (.not.this%phi(n)%has_prop('specific-enthalpy')) then
        !! Set the reference temperature and specific enthalpy.
        if (n > 1) then
          if (.not.allocated(this%pc_seq(n-1)%pc%latent_heat)) then
            stat = 1
            errmsg = 'missing latent heat for phase change ' // &
                     this%phi(n-1)%name // ':' // this%phi(n)%name
            return
          end if
          h0 = h0 + this%pc_seq(n-1)%pc%latent_heat
        end if

        !! Integrate the specific heat to get the specific enthalpy.
        call this%phi(n)%get_prop('specific-heat', cp)
        if (.not.allocated(cp)) then
          stat = 1
          errmsg = 'missing specific-heat for phase ' // this%phi(n)%name
          return
        end if
        call alloc_scalar_func_antideriv(cp, t0, h0, h, stat, errmsg)
        if (stat /= 0) then
          errmsg = 'unable to integrate specific heat for phase ' // &
                   this%phi(n)%name // ': ' // errmsg
          return
        end if
        call this%phi(n)%add_prop('specific-enthalpy', h)
      end if

      !! Set the reference temperature and specific enthalpy for the next phase
      if (n < size(this%phi)) then
        call this%phi(n)%get_prop('specific-enthalpy', h)
        t0 = this%pc_seq(n)%pc%ref_liquidus_temp()
        h0 = h%eval([t0])
      end if

    end do

  end subroutine add_specific_enthalpy_prop

  !! This subroutine attempts to build the enthalpy property for the material
  !! by forming the product of its specific enthalpy and density properties,
  !! if possible. Where the specific enthalpy property does not exist, it will
  !! attempt to build it too, if possible, using specific heat properties and
  !! and phase change latent heats (see ADD_SPECIFIC_ENTHALPY_PROP). STAT
  !! returns a non-zero value if the enthalpy property is not ultimately
  !! defined, and ERRMSG returns an explanatory message.
  !!
  !! NB: If the enthalpy property is already defined, the subroutine simply
  !! returns with no action taken. It is an error, however, for some phases
  !! to have the enthalpy property but others not.
  !! TODO: Is there some way to make this external to the class?

  subroutine add_enthalpy_prop(this, stat, errmsg)

    use scalar_func_tools

    class(multi_phase_matl), intent(inout) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    class(scalar_func), allocatable :: cp, h, rho, hd
    real(r8) :: t0, h0
    integer :: n

    if (this%has_prop('enthalpy')) return

    if (.not.this%has_prop('specific-enthalpy')) then
      call add_specific_enthalpy_prop(this, stat, errmsg)
      if (stat /= 0) return
    end if

    !! Add per-phase enthalpy property where indicated.
    do n = 1, size(this%phi)
      if (this%phi(n)%has_prop('enthalpy')) then
        stat = 1
        errmsg = 'attempting to redefine enthalpy for phase ' // this%phi(n)%name
        return
      end if
      if (.not.this%phi(n)%has_prop('density')) then
        stat = 1
        errmsg = 'missing density for phase ' // this%phi(n)%name
        return
      end if
      if (.not.this%phi(n)%has_prop('density', strict=.true.) .and. &
          .not.this%phi(n)%has_prop('specific-enthalpy', strict=.true.)) cycle
      call this%phi(n)%get_prop('specific-enthalpy', h)
      call this%phi(n)%get_prop('density', rho)
      call alloc_scalar_func_product(rho, h, hd, stat, errmsg)
      INSIST(stat == 0)
      call this%phi(n)%add_prop('enthalpy', hd)
    end do

    !! Add material-wide enthalpy property if indicated.
    if (this%phase%has_prop('specific-enthalpy') .and. this%phase%has_prop('density')) then
      call this%phase%get_prop('specific-enthalpy', h)
      call this%phase%get_prop('density', rho)
      call alloc_scalar_func_product(rho, h, hd, stat, errmsg)
      INSIST(stat == 0)
      call this%add_prop('enthalpy', hd)
    end if

    ASSERT(this%has_prop('enthalpy'))

  end subroutine add_enthalpy_prop

  subroutine alloc_matl_prop(this, name, prop, errmsg)

    class(multi_phase_matl), intent(in), target :: this
    character(*), intent(in) :: name
    class(matl_prop), allocatable, intent(out) :: prop
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    type(multiphase_prop), allocatable :: p

    allocate(p)
    n = size(this%phi)
    allocate(p%farray(n))
    p%matl => this

    do n = 1, size(this%phi)
      call this%phi(n)%get_prop(name, p%farray(n)%func)
      if (.not.allocated(p%farray(n)%func)) then
        errmsg = name // ' not defined for phase ' // this%phi(n)%name
        return
      end if
    end do
    call move_alloc(p, prop)

  end subroutine alloc_matl_prop

  subroutine compute_value(this, state, value)
    class(multiphase_prop), intent(in) :: this
    real(r8), intent(in) :: state(:)
    real(r8), intent(out) :: value
    integer :: n
    real(r8) :: fs
    do n = 1, size(this%matl%pc_seq)
      associate (pc => this%matl%pc_seq(n)%pc)
        if (state(1) <= pc%solidus_temp()) then
          value = this%farray(n)%func%eval(state)
          return
        else if (state(1) < pc%liquidus_temp()) then
          fs = pc%solid_frac(state(1))
          value = fs*this%farray(n)%func%eval(state) + &
                  (1-fs)*this%farray(n+1)%func%eval(state)
          return
        end if
      end associate
    end do
    value = this%farray(n)%func%eval(state)
  end subroutine compute_value

  subroutine write_solid_frac_plotfile(this, n, filename, digits, npoints, iostat)
    class(multi_phase_matl), intent(in) :: this
    character(*), intent(in) :: filename
    integer, intent(in) :: n, digits, npoints
    integer, intent(out) :: iostat
    ASSERT(n > 0 .and. n <= size(this%pc_seq))
    call this%pc_seq(n)%pc%write_solid_frac_plotfile(filename, digits, npoints, iostat)
  end subroutine

end module multi_phase_matl_type
