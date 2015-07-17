!!
!! USTRUC_CORE_TYPE
!!
!! A concrete implementation the abstract base class USTRUC_CORE.  This
!! implementation defines the non-optional core microstructure analysis
!! component that is referenced by the optional analysis components.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! August 2014
!!
!! PROGRAMMING INTERFACE
!!
!!  See the comments accompanying the source of the abstract base class
!!  USTRUC_COMP.  The design of low-level microstructure model is implemented
!!  using the decorator programming pattern, whose "concrete component" class
!!  corresponds to this derived type.  All of the optional components will
!!  contain a reference to an instance of this type through which they can
!!  access the common state variables which are stored as array components.
!!
!!  An object of this type should be instantiated using the function
!!  NEW_USTRUC_CORE(N) which returns a pointer to a new TYPE(USTRUC_CORE)
!!  object with state array components allocated for N independent points.
!!  The content of the arrays is initialized by a subsequent call to
!!  SET_STATE.
!!
!!  Objects of this type respond to the following data names in the generic
!!  GET subroutine: 'temp', 'temp-grad', 'frac', 'frac-grad', and 'invalid'.
!!  These correspond exactly to the data passed in the SET_STATE and
!!  UPDATE_STATE calls.  They also respond to 'frac-rate' which returns the
!!  time rate of change of solid fraction, which is computed using a simple
!!  backward time difference between successive states.  Before the first
!!  state update, it returns a dummy value of 0.
!!
!! NOTES
!!
!! The design of low-level microstructure model is implemented using the
!! decorator programming pattern, whose "concrete component" class corresponds
!! to this derived type.
!!
!! All of the specific procedures for the GET method are defined.  If at this
!! point (the end of the chain of analysis components) the specified data name
!! is not recognized a fatal error message is written and execution is halted.
!!
!! The internal state includes some arrays that are managed by one of the
!! optional analysis components that computes solidification velocities.
!! This is rather funky, and in retrospect perhaps not such a great idea.
!! The idea is that we want to be able to have different ways to compute
!! this velocity, and treating them as "decorators" was a really simple way
!! to do this.  But because that data is needed by other analysis components
!! it needed to be stored with the core component rather than private to the
!! decorator as would be normal.  It might be better to incorporate an
!! option-driven velocity computation directly into this core component.
!! As it stands now, we must be careful to include the velocity component at
!! the appropriate position in the chain of components (next to last) whenever
!! velocity data is required.
!!

#include "f90_assert.fpp"

module ustruc_core_type

  use kinds, only: r8
  use ustruc_comp_class
  use truchas_logging_services
  implicit none
  private

  public :: new_ustruc_core

  type, extends(ustruc_comp), public :: ustruc_core
    real(r8) :: t = 0.0_r8
    real(r8), allocatable :: temp(:)
    real(r8), allocatable :: temp_grad(:,:)
    real(r8), allocatable :: frac(:)
    real(r8), allocatable :: frac_grad(:,:)
    real(r8), allocatable :: frac_rate(:)
    logical,  allocatable :: invalid(:)
    !! These are managed by the USTRUC_VELOCITY type
    real(r8), allocatable :: speed(:)
    real(r8), allocatable :: velocity(:,:)
    logical,  allocatable :: invalid_velocity(:)
  contains
    procedure, private :: init
    procedure :: set_state
    procedure :: update_state
    procedure :: has
    procedure :: getl1
    procedure :: geti1
    procedure :: getr1
    procedure :: getr2
  end type ustruc_core

contains

  function new_ustruc_core (n) result (this)
    integer, intent(in) :: n
    type(ustruc_core), pointer :: this
    allocate(this)
    call this%init (n)
  end function

  subroutine init (this, n)
    integer, intent(in) :: n
    class(ustruc_core), intent(out) :: this
    ASSERT(n >= 0)
    this%n = n
    allocate(this%temp(n), this%temp_grad(3,n))
    allocate(this%frac(n), this%frac_grad(3,n), this%frac_rate(n))
    allocate(this%invalid(n))
  end subroutine init

  subroutine set_state (this, t, temp, temp_grad, frac, frac_grad, invalid)

    class(ustruc_core), intent(inout) :: this
    real(r8), intent(in) :: t, temp(:), temp_grad(:,:), frac(:), frac_grad(:,:)
    logical,  intent(in) :: invalid(:)

    ASSERT(size(temp) == this%n)
    ASSERT(size(temp_grad,1) == 3)
    ASSERT(size(temp_grad,2) == size(temp))
    ASSERT(size(frac) == this%n)
    ASSERT(size(frac_grad,1) == 3)
    ASSERT(size(frac_grad,2) == size(frac))
    ASSERT(size(invalid) == this%n)

    this%t = t
    this%temp = temp
    this%temp_grad = temp_grad
    this%frac = frac
    this%frac_grad = frac_grad
    this%frac_rate = 0.0_r8 ! no valid data here
    this%invalid = invalid

  end subroutine set_state

  subroutine update_state (this, t, temp, temp_grad, frac, frac_grad, invalid)

    class(ustruc_core), intent(inout) :: this
    real(r8), intent(in) :: t, temp(:), temp_grad(:,:), frac(:), frac_grad(:,:)
    logical,  intent(in) :: invalid(:)

    real(r8) :: dt

    ASSERT(size(temp) == this%n)
    ASSERT(size(temp_grad,1) == 3)
    ASSERT(size(temp_grad,2) == size(temp))
    ASSERT(size(frac) == this%n)
    ASSERT(size(frac_grad,1) == 3)
    ASSERT(size(frac_grad,2) == size(frac))
    ASSERT(size(invalid) == this%n)

    dt = t - this%t
    this%t = t

    this%temp = temp
    this%temp_grad = temp_grad
    this%frac_rate = (frac - this%frac)/dt
    where (invalid) this%frac_rate = 0.0_r8
    this%frac = frac
    this%frac_grad = frac_grad

    this%invalid = invalid

  end subroutine update_state

  logical function has (this, name)
    class(ustruc_core), intent(in) :: this
    character(*), intent(in) :: name
    select case (name)
    case ('invalid', 'temp', 'frac-rate', 'temp-grad', 'frac-grad')
      has = .true.
    case default
      has = .false.
    end select
  end function has

  subroutine getl1 (this, name, array)
    class(ustruc_core), intent(in) :: this
    character(*), intent(in) :: name
    logical, intent(out) :: array(:)
    select case (name)
    case ('invalid')
      ASSERT(size(array) == size(this%invalid))
      array = this%invalid
    case default
      call TLS_fatal ('USTRUCT_COMP%GET: unknown name: ' // name)
    end select
  end subroutine getl1

  subroutine geti1 (this, name, array, invalid)
    class(ustruc_core), intent(in) :: this
    character(*), intent(in) :: name
    integer, intent(out) :: array(:)
    logical, intent(out), optional :: invalid(:)
    call TLS_fatal ('USTRUCT_COMP%GET: unknown name: ' // name)
  end subroutine geti1

  subroutine getr1 (this, name, array, invalid)
    class(ustruc_core), intent(in) :: this
    character(*), intent(in) :: name
    real(r8), intent(out) :: array(:)
    logical, intent(out), optional :: invalid(:)
    select case (name)
    case ('temp')
      ASSERT(size(array) == this%n)
      array = this%temp
      if (present(invalid)) then
        ASSERT(size(invalid) == this%n)
        invalid = this%invalid
      else
        where (this%invalid)  array = 0.0_r8
      end if
    case ('frac')
      ASSERT(size(array) == this%n)
      array = this%frac
      if (present(invalid)) then
        ASSERT(size(invalid) == this%n)
        invalid = this%invalid
      else
        where (this%invalid)  array = 0.0_r8
      end if
    case ('frac-rate')
      ASSERT(size(array) == this%n)
      array = this%frac_rate
      if (present(invalid)) then
        ASSERT(size(invalid) == this%n)
        invalid = this%invalid
      else
        where (this%invalid)  array = 0.0_r8
      end if
    case default
      call TLS_fatal ('USTRUCT_COMP%GET: unknown name: ' // name)
    end select
  end subroutine getr1

  subroutine getr2 (this, name, array, invalid)
    class(ustruc_core), intent(in) :: this
    character(*), intent(in) :: name
    real(r8), intent(out) :: array(:,:)
    logical, intent(out), optional :: invalid(:)
    integer :: j
    select case (name)
    case ('temp-grad')
      ASSERT(all(shape(array) == shape(this%temp_grad)))
      array = this%temp_grad
      if (present(invalid)) then
        ASSERT(size(invalid) == this%n)
        invalid = this%invalid
      else
        do j = 1, this%n
          if (this%invalid(j)) array(:,j) = 0.0_r8
        end do
      end if
    case ('frac-grad')
      ASSERT(all(shape(array) == shape(this%frac_grad)))
      array = this%frac_grad
      if (present(invalid)) then
        ASSERT(size(invalid) == this%n)
        invalid = this%invalid
      else
        do j = 1, this%n
          if (this%invalid(j)) array(:,j) = 0.0_r8
        end do
      end if
    case default
      call TLS_fatal ('USTRUCT_COMP%GET: unknown name: ' // name)
    end select
  end subroutine getr2

end module ustruc_core_type
