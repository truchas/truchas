!!
!! USTRUC_GV0
!!
!! A concrete implementation of USTRUC_PLUGIN/USTRUC_COMP that adds the
!! identification of the thermal gradient (G) and solidification front
!! velocity (V) at the onset of solidification.  View this as a starting
!! point for implementing a microstructure prediction model.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! August 2014; updated July 2015
!!
!! PROGRAMMING INTERFACE
!!
!!  See the comments that accompany the source for the USTRUC_COMP and
!!  USTRUC_PLUGIN classes which define the interface to this object.  The
!!  only public entity provided by the module is the following function.
!!
!!  NEW_USTRUC_GV0(COMP, PARAMS) returns a pointer to a new USTRUC_COMP class
!!    object whose dynamic type is USTRUC_GV0.  The USTRUC_GV0 type is itself
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
!!      'theta-gv' -- the solid fraction at which G and V are taken
!!
!!  Objects of this type respond to the following data names in the generic
!!  GET subroutine: 'solid-time', 'g', and 'v'.
!!

#include "f90_assert.fpp"

module ustruc_gv0_type

  use kinds, only: r8
  use,intrinsic :: iso_fortran_env, only: int8, int16
  use ustruc_plugin_class
  implicit none
  private

  public :: new_ustruc_gv0

  type, extends(ustruc_plugin) :: ustruc_gv0
    real(r8) :: f1, f1p, f2, f2p, theta
    real(r8), allocatable :: dt(:), g(:), v(:)
    integer(int8), allocatable :: state(:), gv_state(:)
    integer(int16), allocatable :: count(:)
  contains
    procedure :: set_state
    procedure :: update_state
    procedure :: get_comp_list
    procedure :: has
    procedure :: getl1
    procedure :: getr1
    procedure :: serialize
    procedure :: deserialize
  end type ustruc_gv0

  integer, parameter :: STATE_INVALID   = 0
  integer, parameter :: STATE_UNDEFINED = 1
  integer, parameter :: STATE_LIQUID    = 2
  integer, parameter :: STATE_MUSHY     = 3
  integer, parameter :: STATE_SOLID     = 4

  integer, parameter :: GV_INVALID   = 0
  integer, parameter :: GV_UNDEFINED = 1
  integer, parameter :: GV_DEFINED   = 2
  
  !! Number of bytes (per cell) of internal state for serialization/deserialization
  type(ustruc_gv0), allocatable :: dummy  ! only use is in the following parameter declaration
  integer, parameter :: NBYTES = storage_size(dummy%dt)/8 + storage_size(dummy%g)/8 + &
                                 storage_size(dummy%v)/8 +  storage_size(dummy%state)/8 + &
                                 storage_size(dummy%gv_state)/8 + storage_size(dummy%count)/8

contains

  function new_ustruc_gv0 (comp, params) result (this)

    use parameter_list_type
    use truchas_logging_services

    class(ustruc_comp), pointer, intent(in) :: comp
    type(parameter_list) :: params
    type(ustruc_gv0), pointer :: this

    integer :: stat
    character(:), allocatable :: errmsg

    allocate(this)
    call this%init (comp)

    call params%get ('theta1', this%f1, stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      call TLS_fatal (errmsg)
    else if (this%f1 <= 0.0 .or. this%f1 >= 1.0) then
      call TLS_fatal ('theta1 must be > 0.0 and < 1.0')
    end if

    call params%get ('theta2', this%f2, stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      call TLS_fatal (errmsg)
    else if (this%f2 <= 0.0 .or. this%f2 >= 1.0) then
      call TLS_fatal ('theta2 must be > 0.0 and < 1.0')
    end if
    if (this%f2 <= this%f1) call TLS_fatal ('theta2 <= theta1')

    call params%get ('theta1p', this%f1p, default=this%f1, stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      call TLS_fatal (errmsg)
    else if (this%f1p > this%f1 .or. this%f1p < 0.0) then
      call TLS_fatal ('theta1p must be >= 0.0 and <= theta1')
    end if

    call params%get ('theta2p', this%f2p, default=this%f2, stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      call TLS_fatal (errmsg)
    else if (this%f1p > this%f1 .or. this%f1p < 0.0) then
      call TLS_fatal ('theta2p must be >= theta1 and <= theta2')
    end if

    call params%get ('theta-gv', this%theta, default=this%f1, stat=stat, errmsg=errmsg)
    if (this%theta < this%f1 .or. this%theta > this%f2) then
      call TLS_fatal ('theta-gv must be >= theta1 and <= theta2')
    end if

    allocate(this%state(this%n), this%count(this%n), this%dt(this%n))
    allocate(this%g(this%n), this%v(this%n), this%gv_state(this%n))

    this%gv_state = GV_INVALID

  end function new_ustruc_gv0

  subroutine set_state (this, t, temp, temp_grad, frac, frac_grad, invalid)

    class(ustruc_gv0), intent(inout) :: this
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

    class(ustruc_gv0), intent(inout) :: this
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
        this%gv_state(j) = GV_INVALID
      else
        associate (temp_grad => this%core%temp_grad, velocity => this%core%velocity, &
                   invalid_velocity => this%core%invalid_velocity)
          select case (this%state(j))
          case (STATE_LIQUID)
            if (frac(j) >= this%f1) then
              this%dt(j) = interp_frac(prev_t, prev_frac(j), t, frac(j), this%f1)
              this%state(j) = STATE_MUSHY
            else if (frac(j) <= this%f2) then
              this%count(j) = 1
            end if
            if (frac(j) > this%f2) then
              this%dt(j) = interp_frac(prev_t, prev_frac(j), t, frac(j), this%f2) - this%dt(j)
              this%state(j) = STATE_SOLID
            end if
            if (this%state(j) == STATE_MUSHY) then
              if (frac(j) >= this%theta) then
                if (invalid_velocity(j)) then
                  this%gv_state(j) = GV_INVALID
                else
                  this%g(j) = this%vector_magnitude(temp_grad(:,j))
                  this%v(j) = this%vector_magnitude(velocity(:,j))
                  this%gv_state(j) = GV_DEFINED
                end if
              else
                this%gv_state(j) = GV_UNDEFINED
              end if
            end if
          case (STATE_MUSHY)
            if (frac(j) < this%f1p) then
              this%state(j) = STATE_LIQUID
              this%state(j) = GV_INVALID
            else if (frac(j) > this%f2) then
              this%dt(j) = interp_frac(prev_t, prev_frac(j), t, frac(j), this%f2) - this%dt(j)
              this%state(j) = STATE_SOLID
              if (this%gv_state(j) == GV_UNDEFINED) this%gv_state(j) = GV_INVALID
            else
              this%count(j) = this%count(j) + 1
              if (this%gv_state(j) == GV_UNDEFINED) then
                if (frac(j) >= this%theta) then
                  if (invalid_velocity(j)) then
                    this%gv_state(j) = GV_INVALID
                  else
                    this%g(j) = this%vector_magnitude(temp_grad(:,j))
                    this%v(j) = this%vector_magnitude(velocity(:,j))
                    this%gv_state(j) = GV_DEFINED
                  end if
                end if
              end if
            end if
          case (STATE_SOLID)
            if (frac(j) < this%f1) then
              this%state(j) = STATE_LIQUID
            else if (frac(j) <= this%f2p) then
              this%state(j) = STATE_UNDEFINED
            end if
            if (this%state(j) /= STATE_SOLID) this%gv_state(j) = GV_INVALID
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

  subroutine get_comp_list (this, list)
    class(ustruc_gv0), intent(in) :: this
    integer, allocatable, intent(out) :: list(:)
    integer, allocatable :: rest(:)
    call this%ustruc_plugin%get_comp_list (rest)
    allocate(list(size(rest)+1))
    list(1) = USTRUC_GV0_ID
    list(2:) = rest
  end subroutine get_comp_list

  logical function has (this, name)
    class(ustruc_gv0), intent(in) :: this
    character(*), intent(in) :: name
    select case (name)
    case ('invalid-gv', 'count', 'solid-time', 'g', 'v')
      has = .true.
    case default
      has = this%ustruc_plugin%has (name)
    end select
  end function has

  subroutine getl1 (this, name, array)
    class(ustruc_gv0), intent(in) :: this
    character(*), intent(in) :: name
    logical, intent(out) :: array(:)
    select case (name)
    case ('invalid-gv')
      ASSERT(size(array) == this%n)
      array = (this%gv_state /= GV_DEFINED)
    case default
      call this%ustruc_plugin%get (name, array)
    end select
  end subroutine getl1

  subroutine getr1 (this, name, array, invalid)
    class(ustruc_gv0), intent(in) :: this
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
      where (this%gv_state == GV_DEFINED)
        array = this%g
      else where
        array = 0.0_r8
      end where
      if (present(invalid)) then
        ASSERT(size(invalid) == this%n)
        invalid = (this%gv_state /= GV_DEFINED)
      end if
    case ('v')
      ASSERT(size(array) == this%n)
      where (this%gv_state == GV_DEFINED)
        array = this%v
      else where
        array = 0.0_r8
      end where
      if (present(invalid)) then
        ASSERT(size(invalid) == this%n)
        invalid = (this%gv_state /= GV_DEFINED)
      end if
    case default
      call this%ustruc_plugin%get (name, array, invalid)
    end select
  end subroutine getr1

  subroutine serialize (this, cid, array)

    use,intrinsic :: iso_fortran_env, only: int8
    use serialization_tools, only: copy_to_bytes

    class(ustruc_gv0), intent(in) :: this
    integer, intent(in) :: cid
    integer(int8), allocatable, intent(out) :: array(:,:)

    integer :: j, offset

    if (cid == USTRUC_GV0_ID) then
      allocate(array(NBYTES,this%n))
      do j = 1, this%n
        offset = 0
        call copy_to_bytes (this%dt(j), array(:,j), offset)
        call copy_to_bytes (this%g(j), array(:,j), offset)
        call copy_to_bytes (this%v(j), array(:,j), offset)
        call copy_to_bytes (this%state(j), array(:,j), offset)
        call copy_to_bytes (this%gv_state(j), array(:,j), offset)
        call copy_to_bytes (this%count(j), array(:,j), offset)
      end do
    else
      call this%ustruc_plugin%serialize (cid, array)
    end if

  end subroutine serialize

  subroutine deserialize (this, cid, array)

    use,intrinsic :: iso_fortran_env, only: int8
    use serialization_tools, only: copy_from_bytes

    class(ustruc_gv0), intent(inout) :: this
    integer, intent(in) :: cid
    integer(int8), intent(in) :: array(:,:)

    integer :: j, offset

    if (cid == USTRUC_GV0_ID) then
      INSIST(size(array,1) == NBYTES)
      INSIST(size(array,2) == this%n)
      do j = 1, this%n
        offset = 0
        call copy_from_bytes (array(:,j), offset, this%dt(j))
        call copy_from_bytes (array(:,j), offset, this%g(j))
        call copy_from_bytes (array(:,j), offset, this%v(j))
        call copy_from_bytes (array(:,j), offset, this%state(j))
        call copy_from_bytes (array(:,j), offset, this%gv_state(j))
        call copy_from_bytes (array(:,j), offset, this%count(j))
      end do
    else
      call this%ustruc_plugin%deserialize (cid, array)
    end if

  end subroutine deserialize

end module ustruc_gv0_type
