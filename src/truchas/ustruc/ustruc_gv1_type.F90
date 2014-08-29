!!
!! USTRUC_GV1
!!
!! A concrete implementation of USTRUC_PLUGIN/USTRUC_COMP that adds the
!! identification of microstructure type and characteristics (first model).
!!
!! NB: THIS IS AN INCOMPLETE STUB.  ONLY THE THERMAL GRADIENT MAGNITUDE (G)
!! AND SOLIDIFICATION FRONT SPEED (V) ARE COMPUTED.  THESE WILL BE THE INPUTS
!! TO A TABLE LOOKUP FOR THE TYPE OF MICROSTRUCTURE (DENDRITIC, PLANAR, ...)
!! AND THE COMPUTATION OF CHARACTERISTICS (PRIMARY AND SECONDARY ARM SPACING)
!! AWAITING THIS INFO FROM SETH IMHOFF AND PAUL GIBBS (MST-6)
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
!!  NEW_USTRUC_GV1(COMP, PARAMS) returns a pointer to a new USTRUC_COMP class
!!    object whose dynamic type is USTRUC_GV1.  The USTRUC_GV1 type is itself
!!    private.  COMP is a USTRUC_COMP pointer whose target is being wrapped
!!    by this new analysis component specified by PARAMS.  The new object
!!    takes ownership of the target.  PARAMS is a PARAMETER_LIST object;
!!    the relevant parameters in PARAMS are:
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
!!  GET subroutine: 'solid-time', 'g', and 'v'.
!!

#include "f90_assert.fpp"

module ustruc_gv1_type

  use kinds, only: r8
  use ustruc_plugin_class
  implicit none
  private

  public :: new_ustruc_gv1

  type, extends(ustruc_plugin) :: ustruc_gv1
    real(r8) :: f1, f1p, f2, f2p
    integer,  allocatable :: state(:)
    real(r8), allocatable :: dt(:), g(:), v(:)
    integer,  allocatable :: count(:)
  contains
    procedure :: set_state
    procedure :: update_state
    procedure :: getl1
    procedure :: getr1
  end type ustruc_gv1

  integer, parameter :: STATE_INVALID   = 0
  integer, parameter :: STATE_UNDEFINED = 1
  integer, parameter :: STATE_LIQUID    = 2
  integer, parameter :: STATE_MUSHY     = 3
  integer, parameter :: STATE_SOLID     = 4

contains

  function new_ustruc_gv1 (comp, params) result (this)

    use parameter_list_type

    class(ustruc_comp), pointer, intent(in) :: comp
    type(parameter_list) :: params
    type(ustruc_gv1), pointer :: this

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

    allocate(this%state(this%n), this%dt(this%n), this%g(this%n), this%v(this%n))
    allocate(this%count(this%n))

  end function new_ustruc_gv1

  subroutine set_state (this, t, temp, temp_grad, frac, frac_grad, invalid)

    class(ustruc_gv1), intent(inout) :: this
    real(r8), intent(in) :: t, temp(:), temp_grad(:,:), frac(:), frac_grad(:,:)
    logical,  intent(in) :: invalid(:)

    integer :: j

    !! Set the initial state for the next analysis component in the chain.
    !! This will define the CORE component state arrays.
    call this%ustruc_plugin%set_state (t, temp, temp_grad, frac, frac_grad, invalid)

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
    this%g = 0.0_r8
    this%v = 0.0_r8
    this%count = 0

  end subroutine set_state

  subroutine update_state (this, t, temp, temp_grad, frac, frac_grad, invalid)

    class(ustruc_gv1), intent(inout) :: this
    real(r8), intent(in) :: t, temp(:), temp_grad(:,:), frac(:), frac_grad(:,:)
    logical,  intent(in) :: invalid(:)

    integer  :: j
    real(r8) :: prev_t, prev_frac(this%n)

    prev_t = this%core%t
    prev_frac = this%core%frac

    !! Update the next analysis component in the chain.
    !! Note that the core analysis component always ends the chain.
    call this%ustruc_plugin%update_state (t, temp, temp_grad, frac, frac_grad, invalid)

    !! Update the state of this analysis component.
    do j = 1, this%n
      if (invalid(j)) then
        this%state(j) = STATE_INVALID
      else
        associate (temp_grad => this%core%temp_grad, velocity => this%core%velocity)
          select case (this%state(j))
          case (STATE_LIQUID)
            if (frac(j) >= this%f1) then
              this%dt(j) = interp_frac(prev_t, prev_frac(j), t, frac(j), this%f1)
              this%state(j) = STATE_MUSHY
              this%g(j) = this%vector_magnitude(temp_grad(:,j))
              this%v(j) = this%vector_magnitude(velocity(:,j))
            else if (frac(j) <= this%f2) then
              this%count(j) = 1
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
            else
              this%count(j) = this%count(j) + 1
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
              this%g(j) = this%vector_magnitude(temp_grad(:,j))
              this%v(j) = this%vector_magnitude(velocity(:,j))
              if (frac(j) > this%f2) then
                this%dt(j) = interp_frac(prev_t, 0.0_r8, t, frac(j), this%f2) - this%dt(j)
                this%state(j) = STATE_SOLID
              else
                this%count(j) = 1
              end if
            end if
          case default
            INSIST(.false.)
          end select
        end associate
      end if
    end do

  contains

    pure function interp_frac (t1, f1, t2, f2, f) result (t)
      real(r8), intent(in) :: t1, f1, t2, f2, f
      real(r8) :: t
      t = t1*((f-f1)/(f2-f1)) + t2*((f2-f)/(f2-f1))
    end function

  end subroutine update_state

  subroutine getl1 (this, name, array)
    class(ustruc_gv1), intent(in) :: this
    character(*), intent(in) :: name
    logical, intent(out) :: array(:)
    select case (name)
    case ('invalid-gv')
      ASSERT(size(array) == this%n)
      array = (this%state /= STATE_SOLID)
    case default
      call this%ustruc_plugin%get (name, array)
    end select
  end subroutine getl1

  subroutine getr1 (this, name, array, invalid)
    class(ustruc_gv1), intent(in) :: this
    character(*), intent(in) :: name
    real(r8), intent(out) :: array(:)
    logical, intent(out), optional :: invalid(:)
    select case (name)
    case ('count')
      ASSERT(size(array) == this%n)
      where (this%state == STATE_SOLID)
        array = this%count
      else where
        array = 0
      end where
      if (present(invalid)) then
        ASSERT(size(invalid) == this%n)
        invalid = (this%state /= STATE_SOLID)
      end if
    case ('solid-time')
      ASSERT(size(array) == this%n)
      where (this%state == STATE_SOLID)
        array = this%dt
      else where
        array = 0.0_r8
      end where
      if (present(invalid)) then
        ASSERT(size(invalid) == this%n)
        invalid = (this%state /= STATE_SOLID)
      end if
    case ('g')
      ASSERT(size(array) == this%n)
      where (this%state == STATE_SOLID)
        array = this%g
      else where
        array = 0.0_r8
      end where
      if (present(invalid)) then
        ASSERT(size(invalid) == this%n)
        invalid = (this%state /= STATE_SOLID)
      end if
    case ('v')
      ASSERT(size(array) == this%n)
      where (this%state == STATE_SOLID)
        array = this%v
      else where
        array = 0.0_r8
      end where
      if (present(invalid)) then
        ASSERT(size(invalid) == this%n)
        invalid = (this%state /= STATE_SOLID)
      end if
    case default
      call this%ustruc_plugin%get (name, array, invalid)
    end select
  end subroutine getr1

end module ustruc_gv1_type
