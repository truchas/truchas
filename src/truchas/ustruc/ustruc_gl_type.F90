!!
!! USTRUC_GL
!!
!! A concrete implementation of USTRUC_COMP/USTRUC_ANALYSIS that adds the
!! identification of the thermal gradient (G) and cooling rate (L) at the
!! onset of solidification, and the duration of solidification.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! November 2020; a modification of GV0 from 2015.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module ustruc_gl_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use,intrinsic :: iso_fortran_env, only: int8
  use ustruc_comp_class
  implicit none
  private

  public :: new_ustruc_gl

  type, extends(ustruc_comp), abstract :: ustruc_gl
    real(r8), allocatable :: dt(:), g(:,:), l(:)
    integer(int8), allocatable :: state(:), gl_state(:)
  contains
    procedure :: get_comp_list
    procedure :: has
    procedure :: getl1
    procedure :: getr1
    procedure :: getr2
    procedure :: serialize
    procedure :: deserialize
  end type ustruc_gl

  !! Specific implementation using solid fraction based thresholds
  type, extends(ustruc_gl) :: ustruc_gl_frac
    real(r8) :: begin_frac, end_frac, gl_frac, begin_frac_reset, end_frac_reset
  contains
    procedure :: set_state => set_state_frac
    procedure :: update_state => update_state_frac
  end type

  !! Specific implementation using temperature based thresholds
  type, extends(ustruc_gl) :: ustruc_gl_temp
    real(r8) :: begin_temp, end_temp, gl_temp, begin_temp_reset, end_temp_reset
  contains
    procedure :: set_state => set_state_temp
    procedure :: update_state => update_state_temp
  end type

  integer, parameter :: STATE_INVALID   = 0
  integer, parameter :: STATE_UNDEFINED = 1
  integer, parameter :: STATE_LIQUID    = 2
  integer, parameter :: STATE_MUSHY     = 3
  integer, parameter :: STATE_SOLID     = 4

  integer, parameter :: GL_INVALID   = 0
  integer, parameter :: GL_UNDEFINED = 1
  integer, parameter :: GL_DEFINED   = 2

#ifndef INTEL_BUG20200721
  !! Number of bytes (per cell) of internal state for serialization/deserialization
  type(ustruc_gl_frac), allocatable :: dummy  ! only use is in the following parameter declaration
  integer, parameter :: NBYTES = storage_size(dummy%dt)/8 + 3*storage_size(dummy%g)/8 + &
      storage_size(dummy%l)/8 + storage_size(dummy%state)/8 + storage_size(dummy%gl_state)/8
#endif

contains

  function new_ustruc_gl(comp, params) result(this)

    use parameter_list_type
    use truchas_logging_services

    class(ustruc_analysis), pointer, intent(in) :: comp
    type(parameter_list) :: params
    class(ustruc_gl), pointer :: this

    if (params%is_parameter('begin-frac')) then ! solid fraction based thresholds
      this => new_ustruc_gl_frac(comp, params)
    else if (params%is_parameter('begin-temp')) then ! temperature based thresholds
      this => new_ustruc_gl_temp(comp, params)
    else
      call TLS_fatal('either begin-frac or begin-temp must be specified')
    end if

  end function

  function new_ustruc_gl_frac(comp, params) result(this)

    use parameter_list_type
    use truchas_logging_services

    class(ustruc_analysis), pointer, intent(in) :: comp
    type(parameter_list) :: params
    type(ustruc_gl_frac), pointer :: this

    integer :: stat
    character(:), allocatable :: errmsg

    allocate(this)
    call this%init(comp)

    allocate(this%state(this%n), this%gl_state(this%n))
    allocate(this%dt(this%n), this%g(3,this%n), this%l(this%n))

    this%gl_state = GL_INVALID

    call params%get('begin-frac', this%begin_frac, stat, errmsg)
    if (stat /= 0) then
      call TLS_fatal(errmsg)
    else if (this%begin_frac <= 0 .or. this%begin_frac >= 1) then
      call TLS_fatal('begin-frac must belong to (0, 1)')
    end if

    call params%get('end-frac', this%end_frac, stat, errmsg)
    if (stat /= 0) then
      call TLS_fatal(errmsg)
    else if (this%end_frac <= this%begin_frac .or. this%end_frac >= 1) then
      call TLS_fatal('end-frac must belong to (begin-frac, 1)')
    end if

    call params%get('gl-frac', this%gl_frac, stat, errmsg, default=this%begin_frac)
    if (stat /= 0) then
      call TLS_fatal(errmsg)
    else if (this%gl_frac < this%begin_frac .or. this%gl_frac >= this%end_frac) then
      call TLS_fatal('gl-frac must belong to [begin-frac, end-frac)')
    end if

    call params%get('begin-frac-reset', this%begin_frac_reset, stat, errmsg, default=this%begin_frac)
    if (stat /= 0) then
      call TLS_fatal(errmsg)
    else if (this%begin_frac_reset <= 0 .or. this%begin_frac_reset > this%begin_frac) then
      call TLS_fatal('begin-frac-reset must belong to (0, begin-frac]')
    end if

    call params%get('end-frac-reset', this%end_frac_reset, stat, errmsg, default=this%end_frac)
    if (stat /= 0) then
      call TLS_fatal(errmsg)
    else if (this%end_frac_reset <= this%begin_frac .or. this%end_frac_reset > this%end_frac) then
      call TLS_fatal('end-frac-reset must belong to (begin-frac, end-frac]')
    end if

  end function new_ustruc_gl_frac

  function new_ustruc_gl_temp(comp, params) result(this)

    use parameter_list_type
    use truchas_logging_services

    class(ustruc_analysis), pointer, intent(in) :: comp
    type(parameter_list) :: params
    type(ustruc_gl_temp), pointer :: this

    integer :: stat
    character(:), allocatable :: errmsg

    allocate(this)
    call this%init(comp)

    allocate(this%state(this%n), this%gl_state(this%n))
    allocate(this%dt(this%n), this%g(3,this%n), this%l(this%n))

    this%gl_state = GL_INVALID

    call params%get('begin-temp', this%begin_temp, stat, errmsg)
    if (stat /= 0) call TLS_fatal(errmsg)

    call params%get('end-temp', this%end_temp, stat, errmsg)
    if (stat /= 0) then
      call TLS_fatal(errmsg)
    else if (this%end_temp >= this%begin_temp) then
      call TLS_fatal('end-temp must be < begin-temp')
    end if

    call params%get('gl-temp', this%gl_temp, stat, errmsg, default=this%begin_temp)
    if (stat /= 0) then
      call TLS_fatal(errmsg)
    else if (this%gl_temp <= this%end_temp .or. this%gl_temp > this%begin_temp) then
      call TLS_fatal('gl-temp must belong to (end-temp, begin-temp]')
    end if

    call params%get('begin-temp-reset', this%begin_temp_reset, stat, errmsg, default=this%begin_temp)
    if (stat /= 0) then
      call TLS_fatal(errmsg)
    else if (this%begin_temp_reset < this%begin_temp) then
      call TLS_fatal('begin-temp-reset must be > begin-temp')
    end if

    call params%get('end-temp-reset', this%end_temp_reset, stat, errmsg, default=this%end_temp)
    if (stat /= 0) then
      call TLS_fatal(errmsg)
    else if (this%end_temp_reset < this%end_temp) then
      call TLS_fatal('end-temp-reset must be > end-temp')
    end if

  end function new_ustruc_gl_temp

  subroutine set_state_frac(this, t, temp, temp_grad, frac, invalid)

    class(ustruc_gl_frac), intent(inout) :: this
    real(r8), intent(in) :: t, temp(:), temp_grad(:,:), frac(:)
    logical,  intent(in) :: invalid(:)

    integer :: j

    !! Set the initial state for the next analysis component in the chain.
    !! This will define the CORE component state arrays.
    call this%ustruc_comp%set_state(t, temp, temp_grad, frac, invalid)

    do j = 1, this%n
      if (invalid(j)) then
        this%state(j) = STATE_INVALID
      else if (frac(j) <= this%begin_frac) then
        this%state(j) = STATE_LIQUID
      else
        this%state(j) = STATE_UNDEFINED
      end if
    end do

    !! Assign dummy values to the remaining arrays.
    this%dt = 0.0_r8
    this%g = 0.0_r8
    this%l = 0.0_r8

  end subroutine set_state_frac

  subroutine set_state_temp(this, t, temp, temp_grad, frac, invalid)

    class(ustruc_gl_temp), intent(inout) :: this
    real(r8), intent(in) :: t, temp(:), temp_grad(:,:), frac(:)
    logical,  intent(in) :: invalid(:)

    integer :: j

    !! Set the initial state for the next analysis component in the chain.
    !! This will define the CORE component state arrays.
    call this%ustruc_comp%set_state(t, temp, temp_grad, frac, invalid)

    do j = 1, this%n
      if (invalid(j)) then
        this%state(j) = STATE_INVALID
      else if (temp(j) >= this%begin_temp) then
        this%state(j) = STATE_LIQUID
      else
        this%state(j) = STATE_UNDEFINED
      end if
    end do

    !! Assign dummy values to the remaining arrays.
    this%dt = 0.0_r8
    this%g = 0.0_r8
    this%l = 0.0_r8

  end subroutine set_state_temp

  subroutine update_state_frac(this, t, temp, temp_grad, frac, invalid)

    class(ustruc_gl_frac), intent(inout) :: this
    real(r8), intent(in) :: t, temp(:), temp_grad(:,:), frac(:)
    logical,  intent(in) :: invalid(:)

    integer  :: j

    !! Update the state of this analysis component.
    associate (prev_t => this%core%t, prev_temp => this%core%temp, prev_frac => this%core%frac)
      do j = 1, this%n
        if (invalid(j)) then
          this%state(j) = STATE_INVALID
          this%gl_state(j) = GL_INVALID
        else
          select case (this%state(j))
          case (STATE_LIQUID)
            if (frac(j) >= this%begin_frac) then  ! began solidifying
              this%dt(j) = lin_interp(prev_t, prev_frac(j), t, frac(j), this%begin_frac)
              this%state(j) = STATE_MUSHY
            end if
            if (frac(j) > this%end_frac) then   ! completely solidified
              this%dt(j) = lin_interp(prev_t, prev_frac(j), t, frac(j), this%end_frac) - this%dt(j)
              this%state(j) = STATE_SOLID
            end if
            if (this%state(j) == STATE_MUSHY) then
              if (frac(j) >= this%gl_frac) then ! crossed threshold for picking off G and L data
                this%g(:,j) = temp_grad(:,j)
                this%l(j) = (prev_temp(j) - temp(j))/(t - prev_t)
                this%gl_state(j) = GL_DEFINED
              else
                this%gl_state(j) = GL_UNDEFINED
              end if
            end if
          case (STATE_MUSHY)
            if (frac(j) < this%begin_frac_reset) then ! re-melted
              this%state(j) = STATE_LIQUID
              this%state(j) = GL_INVALID
            else if (frac(j) > this%end_frac) then  ! completely solidified
              this%dt(j) = lin_interp(prev_t, prev_frac(j), t, frac(j), this%end_frac) - this%dt(j)
              this%state(j) = STATE_SOLID
              if (this%gl_state(j) == GL_UNDEFINED) this%gl_state(j) = GL_INVALID
            else  ! still mushy
              if (this%gl_state(j) == GL_UNDEFINED) then
                if (frac(j) >= this%gl_frac) then ! crossed threshold for picking off G and L data
                  this%g(:,j) = temp_grad(:,j)
                  this%l(j) = (prev_temp(j) - temp(j))/(t - prev_t)
                  this%gl_state(j) = GL_DEFINED
                end if
              end if
            end if
          case (STATE_SOLID)
            if (frac(j) < this%begin_frac) then ! completely re-melted
              this%state(j) = STATE_LIQUID
            else if (frac(j) < this%end_frac_reset) then ! partially re-melted
              this%state(j) = STATE_UNDEFINED
            end if
            if (this%state(j) /= STATE_SOLID) this%gl_state(j) = GL_INVALID
          case (STATE_UNDEFINED)
            if (frac(j) < this%begin_frac) this%state(j) = STATE_LIQUID
          case (STATE_INVALID)
            ! A cell containing an insignificant fraction of material is
            ! marked invalid. We are here because over the last time step
            ! the cell has filled with a significant amount of the material.
            if (frac(j) < this%begin_frac) then
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
    call this%ustruc_comp%update_state(t, temp, temp_grad, frac, invalid)

  contains

    pure function lin_interp(t1, x1, t2, x2, x) result(t)
      real(r8), intent(in) :: t1, x1, t2, x2, x
      real(r8) :: t
      t = t1*((x2-x)/(x2-x1)) + t2*((x-x1)/(x2-x1))
    end function

  end subroutine update_state_frac

  subroutine update_state_temp(this, t, temp, temp_grad, frac, invalid)

    class(ustruc_gl_temp), intent(inout) :: this
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
            if (temp(j) <= this%begin_temp) then  ! began solidifying
              this%dt(j) = lin_interp(prev_t, prev_temp(j), t, temp(j), this%begin_temp)
              this%state(j) = STATE_MUSHY
            end if
            if (temp(j) < this%end_temp) then   ! completely solidified
              this%dt(j) = lin_interp(prev_t, prev_temp(j), t, temp(j), this%end_temp) - this%dt(j)
              this%state(j) = STATE_SOLID
            end if
            if (this%state(j) == STATE_MUSHY) then
              if (temp(j) <= this%gl_temp) then ! crossed threshold for picking off G and L data
                this%g(:,j) = temp_grad(:,j)
                this%l(j) = (prev_temp(j) - temp(j))/(t - prev_t)
                this%gl_state(j) = GL_DEFINED
              else
                this%gl_state(j) = GL_UNDEFINED
              end if
            end if
          case (STATE_MUSHY)
            if (temp(j) > this%begin_temp_reset) then ! re-melted
              this%state(j) = STATE_LIQUID
              this%state(j) = GL_INVALID
            else if (temp(j) < this%end_temp) then  ! completely solidified
              this%dt(j) = lin_interp(prev_t, prev_temp(j), t, temp(j), this%end_temp) - this%dt(j)
              this%state(j) = STATE_SOLID
              if (this%gl_state(j) == GL_UNDEFINED) this%gl_state(j) = GL_INVALID
            else  ! still mushy
              if (this%gl_state(j) == GL_UNDEFINED) then
                if (temp(j) <= this%gl_temp) then ! crossed threshold for picking off G and L data
                  this%g(:,j) = temp_grad(:,j)
                  this%l(j) = (prev_temp(j) - temp(j))/(t - prev_t)
                  this%gl_state(j) = GL_DEFINED
                end if
              end if
            end if
          case (STATE_SOLID)
            if (temp(j) > this%begin_temp) then ! completely re-melted
              this%state(j) = STATE_LIQUID
            else if (temp(j) > this%end_temp_reset) then ! partially re-melted
              this%state(j) = STATE_UNDEFINED
            end if
            if (this%state(j) /= STATE_SOLID) this%gl_state(j) = GL_INVALID
          case (STATE_UNDEFINED)
            if (temp(j) > this%begin_temp) this%state(j) = STATE_LIQUID
          case (STATE_INVALID)
            ! A cell containing an insignificant fraction of material is
            ! marked invalid. We are here because over the last time step
            ! the cell has filled with a significant amount of the material.
            if (temp(j) > this%begin_temp) then
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
    call this%ustruc_comp%update_state(t, temp, temp_grad, frac, invalid)

  contains

    pure function lin_interp(t1, x1, t2, x2, x) result(t)
      real(r8), intent(in) :: t1, x1, t2, x2, x
      real(r8) :: t
      t = t1*((x2-x)/(x2-x1)) + t2*((x-x1)/(x2-x1))
    end function

  end subroutine update_state_temp

  subroutine get_comp_list(this, list)
    class(ustruc_gl), intent(in) :: this
    integer, allocatable, intent(out) :: list(:)
    integer, allocatable :: rest(:)
    call this%ustruc_comp%get_comp_list(rest)
    allocate(list(size(rest)+1))
    list(1) = USTRUC_GL_ID
    list(2:) = rest
  end subroutine

  logical function has(this, name)
    class(ustruc_gl), intent(in) :: this
    character(*), intent(in) :: name
    select case (name)
    case ('gl-invalid-GL', 'gl-t_sol', 'gl-G', 'gl-L')
      has = .true.
    case default
      has = this%ustruc_comp%has(name)
    end select
  end function

  subroutine getl1(this, name, array)
    class(ustruc_gl), intent(in) :: this
    character(*), intent(in) :: name
    logical, intent(out) :: array(:)
    select case (name)
    case ('gl-invalid-GL')
      ASSERT(size(array) == this%n)
      array = (this%gl_state /= GL_DEFINED)
    case default
      call this%ustruc_comp%get(name, array)
    end select
  end subroutine

  subroutine getr1(this, name, array, invalid)
    class(ustruc_gl), intent(in) :: this
    character(*), intent(in) :: name
    real(r8), intent(out) :: array(:)
    logical, intent(out), optional :: invalid(:)
    select case (name)
    case ('gl-t_sol')
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
    case ('gl-L')
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
      call this%ustruc_comp%get(name, array, invalid)
    end select
  end subroutine getr1

  subroutine getr2(this, name, array, invalid)
    class(ustruc_gl), intent(in) :: this
    character(*), intent(in) :: name
    real(r8), intent(out) :: array(:,:)
    logical, intent(out), optional :: invalid(:)
    integer :: j
    select case (name)
    case ('gl-G')
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
      call this%ustruc_comp%get(name, array, invalid)
    end select
  end subroutine

  subroutine serialize(this, cid, array)

    use,intrinsic :: iso_fortran_env, only: int8
    use serialization_tools, only: copy_to_bytes

    class(ustruc_gl), intent(in) :: this
    integer, intent(in) :: cid
    integer(int8), allocatable, intent(out) :: array(:,:)

    integer :: j, offset
#ifdef INTEL_BUG20200721
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
      call this%ustruc_comp%serialize(cid, array)
    end if

  end subroutine serialize

  subroutine deserialize(this, cid, array)

    use,intrinsic :: iso_fortran_env, only: int8
    use serialization_tools, only: copy_from_bytes

    class(ustruc_gl), intent(inout) :: this
    integer, intent(in) :: cid
    integer(int8), intent(in) :: array(:,:)

    integer :: j, offset
#ifdef INTEL_BUG20200721
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
      call this%ustruc_comp%deserialize(cid, array)
    end if

  end subroutine deserialize

end module ustruc_gl_type
