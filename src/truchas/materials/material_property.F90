!!
!! MATERIAL_PROPERTY
!!
!! The derived type MAT_PROP is an opaque structure that encapsulates the
!! information required to evaluate a specific thermo-physical property
!! for a specific material as a function of the material's state variables.
!!
!! No defined assignment is provided; do not use instances in an assignment
!! statement unless you really know what the default assignment is doing.
!!
!!  CALL MP_CREATE (THIS, MATERIAL_ID, PROPERTY_ID, STAT, ERRMSG) configures
!!    the MAT_PROP object THIS to describe property PROPERTY_ID for material
!!    MATERIAL_ID.  The integer STAT returns a nonzero value if an error
!!    condition occurred, and an explanatory message is returned in the
!!    character string ERRMSG.  Errors are an invalid material or property ID,
!!    or an unassigned phase property for one or more of the material phases.
!!
!!  CALL MP_EVAL (THIS, STATE, VALUE) evaluates the material property described
!!    by the MAT_PROP object THIS as a function of the state variable values
!!    given in the real rank-1 array STATE and returns the property value in the
!!    real scalar VALUE.
!!
!!  CALL DESTROY (THIS) deallocates any storage associated with the MAT_PROP
!!    object THIS, returning it to its default initialization state.
!!

#include "f90_assert.fpp"

module material_property

  use kinds, only: r8
  use scalar_func_class
  use phase_property_table
  use material_system
  use material_table
  implicit none
  private

  public :: mp_create, mp_eval, mp_eval_deriv, destroy

  type :: box
    integer :: id = 0
    class(scalar_func), allocatable :: prop
  end type box

  type, public :: mat_prop
    private
    integer :: material_id = -1
    integer :: property_id = -1
    integer :: eval_type = -1
    type(mat_system), pointer :: ms => null()
    type(box), allocatable :: phase(:)
  end type mat_prop

  !! Values for EVAL_TYPE.
  integer, parameter :: SINGLE_PHASE = 1
  integer, parameter :: MULTI_PHASE_SINGLE_COMPONENT = 2
  integer, parameter :: MULTI_PHASE_MULTI_COMPONENT = 3

  interface destroy
    module procedure destroy_mat_prop
  end interface

contains

  subroutine mp_create (this, material_id, property_id, stat, errmsg)

    type(mat_prop), intent(out) :: this
    integer, intent(in) :: material_id, property_id
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    integer :: j
    integer, pointer :: phase_id(:) => null()

    if (.not.mt_valid_material(material_id)) then
      stat = -1
      write(errmsg,'(i0,a,i0)') 'invalid material ID: ', material_id
      return
    end if

    if (.not.ppt_valid_property(property_id)) then
      stat = -2
      write(errmsg,'(a,i0)') 'invalid property ID: ', property_id
      return
    end if

    this%material_id = material_id
    this%property_id = property_id

    !! Initialize the list of phase properties.
    this%ms => mt_get_material(material_id)
    ASSERT(associated(this%ms))
    call ms_get_phase_id(this%ms, phase_id)
    allocate(this%phase(size(phase_id)))
    do j = 1, size(phase_id)
      this%phase(j)%id = phase_id(j)
      call ppt_get_phase_property (phase_id(j), property_id, this%phase(j)%prop)
      if (.not.allocated(this%phase(j)%prop)) then
        stat = 1
        errmsg = 'no property for phase: ' // trim(ppt_phase_name(phase_id(j)))
        return
      end if
    end do
    deallocate(phase_id)

    !! Determine the type of property evaluation required by this material system.
    if (ms_num_phase(this%ms) == 1) then
      this%eval_type = SINGLE_PHASE
    else if (ms_num_component(this%ms) == 1) then
      this%eval_type = MULTI_PHASE_SINGLE_COMPONENT
    else
      this%eval_type = MULTI_PHASE_MULTI_COMPONENT
    end if

    stat = 0
    errmsg = ''

  end subroutine mp_create

  subroutine mp_eval (this, state, value)

    type(mat_prop), intent(in) :: this
    real(r8), intent(in)  :: state(:)
    real(r8), intent(out) :: value

    integer :: i
    real(r8), allocatable :: beta(:)

    select case (this%eval_type)
    case (SINGLE_PHASE)

      !! Straight pass-through to the phase property.
      value = this%phase(1)%prop%eval(state)

    case (MULTI_PHASE_SINGLE_COMPONENT)

      !! Temperature-dependent phase fraction average of phase properties.
      allocate(beta(size(this%phase)))
      call ms_phase_mixture (this%ms, state, beta)
      value = 0.0_r8
      do i = 1, size(beta)
        if (beta(i) > 0.0_r8) value = value + beta(i) * this%phase(i)%prop%eval(state)
      end do
      deallocate(beta)

    case (MULTI_PHASE_MULTI_COMPONENT)

      !! State-dependent phase fraction average of phase properties evaluated
      !! at the state-dependent composition of the phase.
      INSIST( .false. ) !!! NOT YET IMPLEMENTED !!!

    end select

  end subroutine mp_eval
  
  subroutine mp_eval_deriv (this, state, n, value)
  
    type(mat_prop), intent(in) :: this
    real(r8), intent(in) :: state(:)
    integer, intent(in) :: n
    real(r8), intent(out) :: value
    
    real(r8) :: fdinc, v1, v2, pstate(size(state))
    
    ASSERT(n >= 1 .and. n <= size(state))
    
    fdinc = max(1.0_r8, abs(state(n))) * sqrt(epsilon(1.0_r8))
    pstate = state
    pstate(n) = state(n) + fdinc
    call mp_eval (this, pstate, v2)
    pstate(n) = state(n) - fdinc
    call mp_eval (this, pstate, v1)
    value = (v2 - v1) / (2*fdinc)
    
  end subroutine mp_eval_deriv

  elemental subroutine destroy_mat_prop (this)
    type(mat_prop), intent(inout) :: this
    type(mat_prop) :: default
    !! N.B.  The MAT_PROP structure does not own the targets of the
    !! MAT_SYSTEM and SCAFUN pointer components and so this routine
    !! must not deallocate/destroy them.
    if (allocated(this%phase)) deallocate(this%phase)
    this = default  ! assign default initialization values
  end subroutine destroy_mat_prop

end module material_property
