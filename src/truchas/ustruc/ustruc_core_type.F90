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
!!  GET subroutine: 'temp', 'temp-grad', 'frac', and 'invalid'.
!!  These correspond exactly to the data passed in the SET_STATE and
!!  UPDATE_STATE calls.  They also respond to 'temp-rate' which returns the
!!  time rate of change of temperature, which is computed using a simple
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

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use ustruc_comp_class
  use truchas_logging_services
  implicit none
  private

  public :: new_ustruc_core

  type, extends(ustruc_comp), public :: ustruc_core
    real(r8) :: t = 0.0_r8
    real(r8), allocatable :: temp(:)
    real(r8), allocatable :: temp_rate(:)   ! R
    real(r8), allocatable :: temp_grad(:,:) ! magnitude is G
    real(r8), allocatable :: frac(:)
    logical,  allocatable :: invalid(:)
    real(r8), allocatable :: speed(:)       ! V
    real(r8), allocatable :: velocity(:,:)
    logical,  allocatable :: invalid_velocity(:)
    real(r8), private :: vmax, theta1, theta2
  contains
    procedure, private :: init
    procedure :: set_state
    procedure :: update_state
    procedure :: get_comp_list
    procedure :: has
    procedure :: getl1
    procedure :: geti1
    procedure :: getr1
    procedure :: getr2
    procedure :: serialize
    procedure :: deserialize
  end type ustruc_core

#ifndef INTEL_BUG20200721
  !! Number of bytes (per cell) of internal state for serialization/deserialization
  type(ustruc_core), allocatable :: dummy  ! only use is in the following parameter declaration
  integer, parameter :: NBYTES = storage_size(dummy%speed)/8 + &
                               3*storage_size(dummy%velocity)/8 + &
                                 storage_size(dummy%invalid_velocity)/8
#endif

contains

  function new_ustruc_core (n, params) result (this)
    use parameter_list_type
    integer, intent(in) :: n
    type(parameter_list) :: params
    type(ustruc_core), pointer :: this
    allocate(this)
    call this%init (n, params)
  end function

  subroutine init (this, n, params)

    use parameter_list_type

    class(ustruc_core), intent(out) :: this
    integer, intent(in) :: n
    type(parameter_list) :: params

    ASSERT(n >= 0)

    this%n = n
    allocate(this%temp(n), this%temp_grad(3,n), this%temp_rate(n))
    allocate(this%frac(n))
    allocate(this%invalid(n))
    allocate(this%invalid_velocity(n), this%velocity(3,n), this%speed(n))

    call params%get ('vel-lo-solid-frac', this%theta1)
    INSIST(this%theta1 > 0.0_r8)
    call params%get ('vel-hi-solid-frac', this%theta2)
    INSIST(this%theta2 < 1.0_r8)
    INSIST(this%theta1 <= this%theta2)
    call params%get ('vel-max', this%vmax)
    INSIST(this%vmax >= 0.0_r8)

  end subroutine init

  subroutine set_state (this, t, temp, temp_grad, frac, invalid)

    class(ustruc_core), intent(inout) :: this
    real(r8), intent(in) :: t, temp(:), temp_grad(:,:), frac(:)
    logical,  intent(in) :: invalid(:)

    ASSERT(size(temp) == this%n)
    ASSERT(size(temp_grad,1) == 3)
    ASSERT(size(temp_grad,2) == size(temp))
    ASSERT(size(frac) == this%n)
    ASSERT(size(invalid) == this%n)

    this%t = t
    this%temp = temp
    this%temp_grad = temp_grad
    this%temp_rate = 0.0_r8 ! no valid data here
    this%frac = frac
    this%invalid = invalid
    this%invalid_velocity = .true.
    this%velocity = 0.0_r8
    this%speed = 0.0_r8

  end subroutine set_state

  subroutine update_state (this, t, temp, temp_grad, frac, invalid)

    class(ustruc_core), intent(inout) :: this
    real(r8), intent(in) :: t, temp(:), temp_grad(:,:), frac(:)
    logical,  intent(in) :: invalid(:)

    integer :: j
    real(r8) :: dt, grad_norm

    ASSERT(size(temp) == this%n)
    ASSERT(size(temp_grad,1) == 3)
    ASSERT(size(temp_grad,2) == size(temp))
    ASSERT(size(frac) == this%n)
    ASSERT(size(invalid) == this%n)

    dt = t - this%t
    this%t = t

    this%temp_rate = (temp - this%temp)/dt
    where (invalid) this%temp_rate = 0.0_r8
    this%temp = temp
    this%temp_grad = temp_grad
    this%frac = frac
    this%invalid = invalid

    do j = 1, this%n
      this%invalid_velocity(j) = .true.
      this%velocity(:,j) = 0.0_r8
      this%speed(j) = 0.0_r8
      if (invalid(j)) cycle
      if (frac(j) >= this%theta1 .and. frac(j) <= this%theta2) then
        grad_norm = this%vector_magnitude(temp_grad(:,j))
        if (abs(this%temp_rate(j)) < this%vmax * grad_norm) then
          this%speed(j) = -this%temp_rate(j) / grad_norm
          this%velocity(:,j) = this%speed(j) * (temp_grad(:,j)/grad_norm)
          this%invalid_velocity(j) = .false.
        end if
      end if
    end do

  end subroutine update_state

  subroutine get_comp_list (this, list)
    class(ustruc_core), intent(in) :: this
    integer, allocatable, intent(out) :: list(:)
    list = [0]
  end subroutine get_comp_list

  logical function has (this, name)
    class(ustruc_core), intent(in) :: this
    character(*), intent(in) :: name
    select case (name)
    case ('invalid', 'temp', 'temp-rate', 'temp-grad', 'speed', 'velocity')
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
      call TLS_fatal ('USTRUC_COMP%GET: unknown name: ' // name)
    end select
  end subroutine getl1

  subroutine geti1 (this, name, array, invalid)
    class(ustruc_core), intent(in) :: this
    character(*), intent(in) :: name
    integer, intent(out) :: array(:)
    logical, intent(out), optional :: invalid(:)
    call TLS_fatal ('USTRUC_COMP%GET: unknown name: ' // name)
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
    case ('temp-rate')
      ASSERT(size(array) == this%n)
      array = this%temp_rate
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
    case ('speed')
      ASSERT(size(array) == this%n)
      array = this%speed
      if (present(invalid)) then
        ASSERT(size(invalid) == this%n)
        invalid = this%invalid_velocity
      end if
    case default
      call TLS_fatal ('USTRUC_COMP%GET: unknown name: ' // name)
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
    case ('velocity')
      ASSERT(all(shape(array) == shape(this%velocity)))
      array = this%velocity
      if (present(invalid)) then
        ASSERT(size(invalid) == this%n)
        invalid = this%invalid_velocity
      end if
    case default
      call TLS_fatal ('USTRUC_COMP%GET: unknown name: ' // name)
    end select
  end subroutine getr2


  subroutine serialize (this, cid, array)

    use,intrinsic :: iso_fortran_env, only: int8
    use serialization_tools, only: copy_to_bytes

    class(ustruc_core), intent(in) :: this
    integer, intent(in) :: cid
    integer(int8), allocatable, intent(out) :: array(:,:)

    integer :: j, offset
#ifdef INTEL_BUG20200721
    integer :: NBYTES
    NBYTES = storage_size(this%speed)/8 + &
           3*storage_size(this%velocity)/8 + &
             storage_size(this%invalid_velocity)/8
#endif

    if (cid == 1) then
      allocate(array(NBYTES,this%n))
      do j = 1, this%n
        offset = 0
        call copy_to_bytes (this%speed(j), array(:,j), offset)
        call copy_to_bytes (this%velocity(:,j), array(:,j), offset)
        call copy_to_bytes (this%invalid_velocity(j), array(:,j), offset)
      end do
    end if

  end subroutine serialize

  subroutine deserialize (this, cid, array)

    use,intrinsic :: iso_fortran_env, only: int8
    use serialization_tools, only: copy_from_bytes

    class(ustruc_core), intent(inout) :: this
    integer, intent(in) :: cid
    integer(int8), intent(in) :: array(:,:)

    integer :: j, offset
#ifdef INTEL_BUG20200721
    integer :: NBYTES
    NBYTES = storage_size(this%speed)/8 + &
           3*storage_size(this%velocity)/8 + &
             storage_size(this%invalid_velocity)/8
#endif

    if (cid == 1) then
      INSIST(size(array,1) == NBYTES)
      INSIST(size(array,2) == this%n)
      do j = 1, this%n
        offset = 0
        call copy_from_bytes (array(:,j), offset, this%speed(j))
        call copy_from_bytes (array(:,j), offset, this%velocity(:,j))
        call copy_from_bytes (array(:,j), offset, this%invalid_velocity(j))
      end do
    end if

  end subroutine deserialize

end module ustruc_core_type
