!!
!! USTRUC_TIME_TYPE
!!
!! A concrete implementation of USTRUC_PLUGIN/USTRUC_COMP that adds the
!! computation of the solidification time.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! August 2014
!!
!! PROGRAMMING INTERFACE
!!
!!  See the comments that accompany the source for the USTRUC_COMP and
!!  USTRUC_PLUGIN classes which define the interface to this object.  The
!!  only public entity provided by the module is the following function.
!!
!!  NEW_USTRUC_TIME(COMP, PARAMS) returns a pointer to a new USTRUC_COMP class
!!    object whose dynamic type is USTRUC_TIME.  The USTRUC_TIME type is itself
!!    private.  COMP is a USTRUC_COMP pointer whose target is being wrapped by
!!    this new analysis component specified by PARAMS.  The new object takes
!!    ownership of the target.  PARAMS is a PARAMETER_LIST object; the relevant
!!    parameters in PARAMS are:
!!
!!      'theta1' -- low solid fraction threshold; solidification is deemed
!!          to have started when the solid fraction crosses this threshold.
!!      'theta2' -- high solid fraction threshold; solidification is deemed
!!          to have finished when the solid fraction crosses this threshold.
!!      'theta1p' -- if the solid fraction drops below this threshold while
!!          solidifying, it is deemed to have returned to a liquid state;
!!          theta1p <= theta1, optional, default theta1p = theta1.
!!      'theta2p' -- if the solid fraction drops below this threshold while
!!          considered solid, it is considered to have remelted and the
!!          previously computed solidification time erased.
!!
!!  Objects of this type respond to the following data names in the generic
!!  GET subroutine: 'solid-time', 'invalid-solid-time'
!!

#include "f90_assert.fpp"

module ustruc_time_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use ustruc_plugin_class
  implicit none
  private

  public :: new_ustruc_time

  type, extends(ustruc_plugin) :: ustruc_time
    real(r8) :: f1, f1p, f2, f2p
    integer,  allocatable :: state(:)
    real(r8), allocatable :: dt(:)
  contains
    procedure :: set_state
    procedure :: update_state
    procedure :: getl1
    procedure :: getr1
  end type ustruc_time

  integer, parameter :: STATE_INVALID   = 0
  integer, parameter :: STATE_UNDEFINED = 1
  integer, parameter :: STATE_LIQUID    = 2
  integer, parameter :: STATE_MUSHY     = 3
  integer, parameter :: STATE_SOLID     = 4

contains

  function new_ustruc_time (comp, params) result (this)

    use parameter_list_type

    class(ustruc_comp), pointer :: comp
    type(parameter_list) :: params
    type(ustruc_time), pointer :: this

    allocate(this)
    call this%init (comp)

    call params%get ('theta1',  this%f1)
    INSIST(this%f1 >= 0.0_r8)
    call params%get ('theta2',  this%f2)
    INSIST(this%f2 <= 1.0_r8)
    INSIST(this%f1 <= this%f2)
    call params%get ('theta1p', this%f1p, default=this%f1)
    INSIST(this%f1p >= 0.0_r8 .and. this%f1p <= this%f1)
    call params%get ('theta2p', this%f2p, default=this%f2)
    INSIST(this%f2p >= this%f1 .and. this%f2p <= this%f2)

    allocate(this%state(this%n), this%dt(this%n))

  end function new_ustruc_time

  subroutine set_state (this, t, temp, temp_grad, frac, invalid)

    class(ustruc_time), intent(inout) :: this
    real(r8), intent(in) :: t, temp(:), temp_grad(:,:), frac(:)
    logical,  intent(in) :: invalid(:)

    integer :: j

    !! Set the initial state for the next analysis component in the chain.
    !! This will define the CORE component state arrays.
    call this%ustruc_plugin%set_state (t, temp, temp_grad, frac, invalid)

    do j = 1, this%n
      if (invalid(j)) then
        this%state(j) = STATE_INVALID
      else if (frac(j) <= this%f1) then
        this%state(j) = STATE_LIQUID
      else
        this%state(j) = STATE_UNDEFINED
      end if
    end do

    !! Assign dummy values to the remaining arrays.
    this%dt = 0.0_r8

  end subroutine set_state

  subroutine update_state (this, t, temp, temp_grad, frac, invalid)

    class(ustruc_time), intent(inout) :: this
    real(r8), intent(in) :: t, temp(:), temp_grad(:,:), frac(:)
    logical,  intent(in) :: invalid(:)

    integer  :: j

    !! Update the state of this analysis component.  We do this before the
    !! the core analysis component is updated (see comments below) so that
    !! the core component state arrays give the state at the previous time.
    do j = 1, this%n
      if (invalid(j)) then
        this%state(j) = STATE_INVALID
      else
        associate (prev_t => this%core%t, prev_frac => this%core%frac)
          select case (this%state(j))
          case (STATE_LIQUID)
            if (frac(j) >= this%f1) then
              this%dt(j) = interp_frac(prev_t, prev_frac(j), t, frac(j), this%f1)
              this%state(j) = STATE_MUSHY
            end if
            if (frac(j) > this%f2) then
              this%dt(j) = interp_frac(prev_t, prev_frac(j), t, frac(j), this%f2) - this%dt(j)
              this%state(j) = STATE_SOLID
            end if
          case (STATE_MUSHY)
            if (frac(j) < this%f1p) then
              this%state(j) = STATE_LIQUID
            else if (frac(j) > this%f2) then
              this%dt(j) = interp_frac(prev_t, prev_frac(j), t, frac(j), this%f2) - this%dt(j)
              this%state(j) = STATE_SOLID
            end if
          case (STATE_SOLID)
            if (frac(j) < this%f1) then
              this%state(j) = STATE_LIQUID
            else if (frac(j) <= this%f2p) then
              this%state(j) = STATE_UNDEFINED
            end if
          case (STATE_UNDEFINED)
            if (frac(j) < this%f1) this%state(j) = STATE_LIQUID
          case (STATE_INVALID)
            if (frac(j) < this%f1) then
              this%state(j) = STATE_LIQUID
            else
              this%dt(j) = interp_frac(prev_t, 0.0_r8, t, frac(j), this%f1)
              this%state(j) = STATE_MUSHY
              if (frac(j) > this%f2) then
                this%dt(j) = interp_frac(prev_t, 0.0_r8, t, frac(j), this%f2) - this%dt(j)
                this%state(j) = STATE_SOLID
              end if
            end if
          case default
            INSIST(.false.)
          end select
        end associate
      end if
    end do

    !! Update the next analysis component in the chain.
    !! Note that the core analysis component always ends the chain.
    call this%ustruc_plugin%update_state (t, temp, temp_grad, frac, invalid)

  contains

    pure function interp_frac (t1, f1, t2, f2, f) result (t)
      real(r8), intent(in) :: t1, f1, t2, f2, f
      real(r8) :: t
      t = t1*((f2-f)/(f2-f1)) + t2*((f-f1)/(f2-f1))
    end function

  end subroutine update_state

  subroutine getl1 (this, name, array)
    class(ustruc_time), intent(in) :: this
    character(*), intent(in) :: name
    logical, intent(out) :: array(:)
    ASSERT(size(array) == this%n)
    select case (name)
    case ('invalid-solid-time')
      array = (this%state /= STATE_SOLID)
    case default
      call this%ustruc_plugin%get (name, array)
    end select
  end subroutine getl1

  subroutine getr1 (this, name, array, invalid)
    class(ustruc_time), intent(in) :: this
    character(*), intent(in) :: name
    real(r8), intent(out) :: array(:)
    logical, intent(out), optional :: invalid(:)
    select case (name)
    case ('solid-time')
      ASSERT(size(array) == this%n)
      where (this%state == STATE_SOLID)
        array = this%dt
      else where
        array = 0.0_r8
      end where
      if (present(invalid)) then
        ASSERT(size(invalid) == this%n)
        invalid = this%state /= STATE_SOLID
      end if
    case default
      call this%ustruc_plugin%get (name, array, invalid)
    end select
  end subroutine getr1

end module ustruc_time_type
