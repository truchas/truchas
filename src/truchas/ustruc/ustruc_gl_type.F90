!!
!! USTRUC_GL
!!
!! A concrete implementation of USTRUC_PLUGIN/USTRUC_COMP that adds the
!! identification of the thermal gradient (G) and cooling rate (L) at the
!! onset of solidification, and the duration of solidification. 
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! November 2020; a modification of GV0 from 2015.
!!
!! PROGRAMMING INTERFACE
!!
!!  See the comments that accompany the source for the USTRUC_COMP and
!!  USTRUC_PLUGIN classes which define the interface to this object.  The
!!  only public entity provided by the module is the following function.
!!
!!  NEW_USTRUC_GL(COMP, PARAMS) returns a pointer to a new USTRUC_COMP class
!!    object whose dynamic type is USTRUC_gl.  The USTRUC_gl type is itself
!!    private.  COMP is a USTRUC_COMP pointer whose target is being wrapped
!!    by this new analysis component specified by PARAMS.  The new object
!!    takes ownership of the target.  PARAMS is a PARAMETER_LIST object;
!!    the relevant parameters in PARAMS are:
!!
!!      'liquidus-temp' -- solidification is deemed to have started when the
!!          temperature drops below this threshold; it need not be exactly
!!          the actual liquidus temperature. Required.
!!      'solidus-temp' -- solidification is deemed to have finished when the
!!          temperature drops below this threshold; it need not be exactly
!!          the actual solidus temperature. Required.
!!      'gl-temp' -- when the temperature drops below this threshold, G and L
!!          are set to the current value of the thermal gradient and cooling
!!          rate. Optional; its default value is the value of liquidus-temp.
!!      'liquidus-temp-reset' -- if the temperature rises back above this
!!          threshold while solidifying, it is deemed to have returned to a
!!          liquid state. Optional; its default value is liquidus-temp, and
!!          liquidus-temp-reset must be >= liquidus-temp.
!!      'solidus-temp-reset' -- if the temperature rises back above this
!!          threshold after becoming solid, it is deemed to have (partially)
!!          re-melted, and the previously computed G, L, and solidification
!!          time are erased. Optional; its default value is solidus-temp, and
!!          solidus-temp-reset must be <= solidus-temp.
!!
!!  Objects of this type respond to the following data names in the generic
!!  GET subroutine: 'solid-time', 'G', and 'L'.
!!

#include "f90_assert.fpp"

module ustruc_gl_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use,intrinsic :: iso_fortran_env, only: int8, int16
  use ustruc_plugin_class
  implicit none
  private

  public :: new_ustruc_gl

  type, extends(ustruc_plugin) :: ustruc_gl
    real(r8) :: Tliq, Tsol, Tgl, Tliq_reset, Tsol_reset
    real(r8), allocatable :: dt(:), g(:,:), l(:)
    integer(int8), allocatable :: state(:), gl_state(:)
  contains
    procedure :: set_state
    procedure :: update_state
    procedure :: get_comp_list
    procedure :: has
    procedure :: getl1
    procedure :: getr1
    procedure :: getr2
    procedure :: serialize
    procedure :: deserialize
  end type ustruc_gl

  integer, parameter :: STATE_INVALID   = 0
  integer, parameter :: STATE_UNDEFINED = 1
  integer, parameter :: STATE_LIQUID    = 2
  integer, parameter :: STATE_MUSHY     = 3
  integer, parameter :: STATE_SOLID     = 4

  integer, parameter :: GL_INVALID   = 0
  integer, parameter :: GL_UNDEFINED = 1
  integer, parameter :: GL_DEFINED   = 2

#ifndef INTEL_COMPILER_WORKAROUND
  !! Number of bytes (per cell) of internal state for serialization/deserialization
  type(ustruc_gl), allocatable :: dummy  ! only use is in the following parameter declaration
  integer, parameter :: NBYTES = storage_size(dummy%dt)/8 + 3*storage_size(dummy%g)/8 + &
      storage_size(dummy%l)/8 + storage_size(dummy%state)/8 + storage_size(dummy%gl_state)/8
#endif

contains

  function new_ustruc_gl(comp, params) result(this)

    use parameter_list_type
    use truchas_logging_services

    class(ustruc_comp), pointer, intent(in) :: comp
    type(parameter_list) :: params
    type(ustruc_gl), pointer :: this

    integer :: stat
    character(:), allocatable :: errmsg

    allocate(this)
    call this%init(comp)

    call params%get('liquidus-temp', this%Tliq, stat=stat, errmsg=errmsg)
    if (stat /= 0) call TLS_fatal(errmsg)

    call params%get('solidus-temp', this%Tsol, stat=stat, errmsg=errmsg)
    if (stat /= 0) call TLS_fatal(errmsg)

    if (this%Tliq <= this%Tsol) call TLS_fatal('liquidus-temp <= solidus-temp')

    call params%get('liquidus-temp-reset', this%Tliq_reset, default=this%Tliq, stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      call TLS_fatal(errmsg)
    else if (this%Tliq_reset < this%Tliq) then
      call TLS_fatal('liquidus-temp-reset < liquidus-temp')
    end if

    call params%get('solidus-temp-reset', this%Tsol_reset, default=this%Tsol, stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      call TLS_fatal(errmsg)
    else if (this%Tsol_reset > this%Tsol) then
      call TLS_fatal('solidus-temp-reset > solidus-temp')
    end if

    call params%get('gl-temp', this%Tgl, default=this%Tliq, stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      call TLS_fatal(errmsg)
    else if (this%Tgl > this%Tliq) then
      call TLS_fatal('gl-temp > liquidus-temp')
    else if (this%Tgl <= this%Tsol) then
      call TLS_fatal('gl-temp <= solidus-temp')
    end if

    allocate(this%state(this%n), this%gl_state(this%n))
    allocate(this%dt(this%n), this%g(3,this%n), this%l(this%n))

    this%gl_state = GL_INVALID

  end function new_ustruc_gl

  subroutine set_state(this, t, temp, temp_grad, frac, invalid)

    class(ustruc_gl), intent(inout) :: this
    real(r8), intent(in) :: t, temp(:), temp_grad(:,:), frac(:)
    logical,  intent(in) :: invalid(:)

    integer :: j

    !! Set the initial state for the next analysis component in the chain.
    !! This will define the CORE component state arrays.
    call this%ustruc_plugin%set_state(t, temp, temp_grad, frac, invalid)

    do j = 1, this%n
      if (invalid(j)) then
        this%state(j) = STATE_INVALID
      else if (temp(j) >= this%Tliq) then
        this%state(j) = STATE_LIQUID
      else
        this%state(j) = STATE_UNDEFINED
      end if
    end do

    !! Assign dummy values to the remaining arrays.
    this%dt = 0.0_r8
    this%g = 0.0_r8
    this%l = 0.0_r8

  end subroutine set_state

  subroutine update_state(this, t, temp, temp_grad, frac, invalid)

    class(ustruc_gl), intent(inout) :: this
    real(r8), intent(in) :: t, temp(:), temp_grad(:,:), frac(:)
    logical,  intent(in) :: invalid(:)

    integer  :: j

    !! Update the state of this analysis component.
    associate (prev_t => this%core%t, prev_temp => this%core%temp)
      do j = 1, this%n
        if (invalid(j)) then
          this%state(j) = STATE_INVALID
          this%gl_state(j) = GL_INVALID
        else
          select case (this%state(j))
          case (STATE_LIQUID)
            if (temp(j) <= this%Tliq) then  ! began solidifying
              this%dt(j) = lin_interp(prev_t, prev_temp(j), t, temp(j), this%Tliq)
              this%state(j) = STATE_MUSHY
            end if
            if (temp(j) < this%Tsol) then   ! completely solidified
              this%dt(j) = lin_interp(prev_t, prev_temp(j), t, temp(j), this%Tsol) - this%dt(j)
              this%state(j) = STATE_SOLID
            end if
            if (this%state(j) == STATE_MUSHY) then
              if (temp(j) <= this%Tgl) then ! crossed threshold for picking off G and L data
                this%g(:,j) = temp_grad(:,j)
                this%l(j) = (prev_temp(j) - temp(j))/(t - prev_t)
                this%gl_state(j) = GL_DEFINED
              else
                this%gl_state(j) = GL_UNDEFINED
              end if
            end if
          case (STATE_MUSHY)
            if (temp(j) > this%Tliq_reset) then ! re-melted
              this%state(j) = STATE_LIQUID
              this%state(j) = GL_INVALID
            else if (temp(j) < this%Tsol) then  ! completely solidified
              this%dt(j) = lin_interp(prev_t, prev_temp(j), t, temp(j), this%Tsol) - this%dt(j)
              this%state(j) = STATE_SOLID
              if (this%gl_state(j) == GL_UNDEFINED) this%gl_state(j) = GL_INVALID
            else  ! still mushy
              if (this%gl_state(j) == GL_UNDEFINED) then
                if (temp(j) <= this%Tgl) then ! crossed threshold for picking off G and L data
                  this%g(:,j) = temp_grad(:,j)
                  this%l(j) = (prev_temp(j) - temp(j))/(t - prev_t)
                  this%gl_state(j) = GL_DEFINED
                end if
              end if
            end if
          case (STATE_SOLID)
            if (temp(j) > this%Tliq) then ! completely re-melted
              this%state(j) = STATE_LIQUID
            else if (temp(j) >= this%Tsol_reset) then ! partially re-melted
              this%state(j) = STATE_UNDEFINED
            end if
            if (this%state(j) /= STATE_SOLID) this%gl_state(j) = GL_INVALID
          case (STATE_UNDEFINED)
            if (temp(j) > this%Tliq) this%state(j) = STATE_LIQUID
          case (STATE_INVALID)
            ! A cell containing an insignificant fraction of material is
            ! marked invalid. We are here because over the last time step
            ! the cell has filled with a significant amount of the material.
            if (temp(j) > this%Tliq) then
              this%state(j) = STATE_LIQUID
            else
              this%state(j) = STATE_UNDEFINED
              !NB: if important, with assumptions we could put in other states.
            end if
          case default
            INSIST(.false.)
          end select
        end if
      end do
    end associate

    !! Update the next analysis component in the chain.
    !! Note that the core analysis component always ends the chain.
    call this%ustruc_plugin%update_state(t, temp, temp_grad, frac, invalid)

  contains

    pure function lin_interp(t1, x1, t2, x2, x) result(t)
      real(r8), intent(in) :: t1, x1, t2, x2, x
      real(r8) :: t
      t = t1*((x2-x)/(x2-x1)) + t2*((x-x1)/(x2-x1))
    end function

  end subroutine update_state

  subroutine get_comp_list(this, list)
    class(ustruc_gl), intent(in) :: this
    integer, allocatable, intent(out) :: list(:)
    integer, allocatable :: rest(:)
    call this%ustruc_plugin%get_comp_list (rest)
    allocate(list(size(rest)+1))
    list(1) = USTRUC_GL_ID
    list(2:) = rest
  end subroutine get_comp_list

  logical function has(this, name)
    class(ustruc_gl), intent(in) :: this
    character(*), intent(in) :: name
    select case (name)
    case ('invalid-gv', 'solid-time', 'G', 'L')
      has = .true.
    case default
      has = this%ustruc_plugin%has(name)
    end select
  end function has

  subroutine getl1(this, name, array)
    class(ustruc_gl), intent(in) :: this
    character(*), intent(in) :: name
    logical, intent(out) :: array(:)
    select case (name)
    case ('invalid-gl')
      ASSERT(size(array) == this%n)
      array = (this%gl_state /= GL_DEFINED)
    case default
      call this%ustruc_plugin%get(name, array)
    end select
  end subroutine getl1

  subroutine getr1(this, name, array, invalid)
    class(ustruc_gl), intent(in) :: this
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
        invalid = (this%state /= STATE_SOLID)
      end if
    case ('L')
      ASSERT(size(array) == this%n)
      where (this%gl_state == GL_DEFINED)
        array = this%l
      else where
        array = 0.0_r8
      end where
      if (present(invalid)) then
        ASSERT(size(invalid) == this%n)
        invalid = (this%gl_state /= GL_DEFINED)
      end if
    case default
      call this%ustruc_plugin%get(name, array, invalid)
    end select
  end subroutine getr1

  subroutine getr2(this, name, array, invalid)
    class(ustruc_gl), intent(in) :: this
    character(*), intent(in) :: name
    real(r8), intent(out) :: array(:,:)
    logical, intent(out), optional :: invalid(:)
    integer :: j
    select case (name)
    case ('G')
      ASSERT(all(shape(array) == shape(this%g)))
      do j = 1, this%n
        if (this%gl_state(j) == GL_DEFINED) then
          array(:,j) = this%g(:,j)
        else
          array(:,j) = 0
        end if
      end do
      if (present(invalid)) then
        ASSERT(size(invalid) == this%n)
        invalid = (this%gl_state /= GL_DEFINED)
      end if
    case default
      call this%ustruc_plugin%get(name, array, invalid)
    end select
  end subroutine getr2

  subroutine serialize(this, cid, array)

    use,intrinsic :: iso_fortran_env, only: int8
    use serialization_tools, only: copy_to_bytes

    class(ustruc_gl), intent(in) :: this
    integer, intent(in) :: cid
    integer(int8), allocatable, intent(out) :: array(:,:)

    integer :: j, offset
#ifdef INTEL_COMPILER_WORKAROUND
    integer :: NBYTES
    NBYTES = storage_size(this%dt)/8 + 3*storage_size(this%g)/8 + &
        storage_size(this%l)/8 + storage_size(this%state)/8 + storage_size(this%gl_state)/8
#endif

    if (cid == USTRUC_GL_ID) then
      allocate(array(NBYTES,this%n))
      do j = 1, this%n
        offset = 0
        call copy_to_bytes(this%dt(j), array(:,j), offset)
        call copy_to_bytes(this%g(:,j), array(:,j), offset)
        call copy_to_bytes(this%l(j), array(:,j), offset)
        call copy_to_bytes(this%state(j), array(:,j), offset)
        call copy_to_bytes(this%gl_state(j), array(:,j), offset)
      end do
    else
      call this%ustruc_plugin%serialize(cid, array)
    end if

  end subroutine serialize

  subroutine deserialize(this, cid, array)

    use,intrinsic :: iso_fortran_env, only: int8
    use serialization_tools, only: copy_from_bytes

    class(ustruc_gl), intent(inout) :: this
    integer, intent(in) :: cid
    integer(int8), intent(in) :: array(:,:)

    integer :: j, offset
#ifdef INTEL_COMPILER_WORKAROUND
    integer :: NBYTES
    NBYTES = storage_size(this%dt)/8 + 3*storage_size(this%g)/8 + &
        storage_size(this%l)/8 +  storage_size(this%state)/8 + storage_size(this%gl_state)/8
#endif

    if (cid == USTRUC_GL_ID) then
      INSIST(size(array,1) == NBYTES)
      INSIST(size(array,2) == this%n)
      do j = 1, this%n
        offset = 0
        call copy_from_bytes(array(:,j), offset, this%dt(j))
        call copy_from_bytes(array(:,j), offset, this%g(:,j))
        call copy_from_bytes(array(:,j), offset, this%l(j))
        call copy_from_bytes(array(:,j), offset, this%state(j))
        call copy_from_bytes(array(:,j), offset, this%gl_state(j))
      end do
    else
      call this%ustruc_plugin%deserialize(cid, array)
    end if

  end subroutine deserialize

end module ustruc_gl_type
