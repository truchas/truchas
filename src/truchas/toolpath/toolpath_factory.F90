!!
!! TOOLPATH_FACTORY
!!
!! This module provides a procedure for instantiating a new TOOLPATH object.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! November 2016
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  A TOOLPATH object is intended to describe the motion through space of some
!!  object or part; for example, a laser or piece of the problem domain. It is
!!  basically just a path function giving position coordinates as a function of
!!  time.  The path is described as a sequence of path segments.
!!
!!  CALL ALLOC_TOOLPATH(THIS, PARAMS, STAT, ERRMSG) allocates and initializes
!!    the TYPE(TOOLPATH) argument THIS. The parameter list PARAMS provides the
!!    data used to initialize the object. The following parameters are recognized.
!!
!!      'start-time' -- starting time of the path (default 0)
!!      'start-coord' -- starting position of the path (default (0,0,0))
!!      'time-scale-factor' -- multiplicative time scale factor (default 1)
!!      'coord-scale-factor' -- multiplicative space scale factor (default 1)
!!      'command-string' -- string containing the list of path commands
!!      'command-file' -- file from which to read the list of path commands
!!
!!    Either 'command-string' or 'command-file' are required.  The path
!!    commands are JSON-format text of the following form:
!!
!!      [<command>, <command>, ...]
!!
!!    where <command> is one of the following:
!!
!!      ["dwell", DT]
!!
!!          Remain at the current position for the time interval DT.
!!
!!      ["moverel", [DX,DY,DZ], S, A, D]
!!
!!          Linear displacement from the current position.  Motion accelerates
!!          from rest to a constant speed, and then decelerates to rest at the
!!          position (DX,DY,DZ) relative to the current position.  The linear
!!          speed, acceleration, and deceleration are S, A, and D, respectively.
!!          If D is omitted, its value is taken to be A.  If both A and D are
!!          omitted then instantaneous acceleration/deceleration to speed S is
!!          assumed.
!!
!!      ["setflag", N1, N2, ...]
!!      ["clrflag", N1, N2, ...]
!!
!!          Sets or clears the indicated flags.  32 flags (0 through 31) are
!!          available.  The setting (set or clear) of a flag holds for all
!!          subsequent motions (above) until changed.  All flags start clear.
!!
!!    Except for the integer flag numbers, real numeric values may be integer
!!    or floating point numbers; i.e., 1 is an acceptable alternative to 1.0
!!
!! IMPLEMENTATION NOTES
!!
!!  The choice to use JSON for the path commands was made to leverage existing
!!  JSON parsing capabilities.  One alternative is to parse the JSON text into
!!  an internal data structure (using json.F90) and then parse it for the path
!!  commands.  That was shaping up to be really messy and clumsy code.  The
!!  alternative approach taken here is to customize the YAJL parsing of the
!!  JSON stream by implementing versions of its call-back functions that
!!  directly implement the path command language (a subset of JSON).  This has
!!  the benefit of being able to validate the input as it is being parsed, and
!!  in the event of an error, being able to point (literally) to its location
!!  in the input.
!!
!!  The implementation of the yajl call-backs for the current path command
!!  language is adequate.  But it appears it will quickly become unwieldy
!!  as more commands are added.  In that event, a promising solution is to
!!  push/pop different sets of call-back functions to the parser that handle
!!  specific path commands.  This was briefly investigated but deemed not yet
!!  needed.
!!

module toolpath_factory

  use toolpath_type
  use parameter_list_type
  use yajl_fort
  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  public :: alloc_toolpath, toolpath

  type, extends(fyajl_callbacks) :: path_builder
    type(toolpath), allocatable :: path
    character(:), allocatable :: errmsg
    integer :: state, flags
    real(r8) :: tsf, csf
    real(r8) :: t, r(3)
    real(r8) :: dt, dr(3), speed, accel, decel
  contains
    procedure :: init
    !! Deferred callback functions from the abstract base type
    procedure :: start_map
    procedure :: end_map
    procedure :: map_key
    procedure :: null_value
    procedure :: logical_value
    procedure :: integer_value
    procedure :: double_value
    procedure :: string_value
    procedure :: start_array
    procedure :: end_array
  end type path_builder

  !! Parser state flags
  enum, bind(c)
    enumerator :: STATE_INIT=0, STATE_DONE, STATE_COMMAND, STATE_ACTION, &
        STATE_SETFLAG1, STATE_SETFLAG2, STATE_CLRFLAG1, STATE_CLRFLAG2,  &
        STATE_DR, STATE_DX, STATE_DY, STATE_DZ, STATE_DR_DONE, &
        STATE_SPEED, STATE_ACCEL, STATE_DECEL, STATE_MOVEREL_DONE, &
        STATE_DWELL_DT, STATE_DWELL_DONE
  end enum

contains

  subroutine alloc_toolpath(path, params, stat, errmsg)

    use,intrinsic :: iso_c_binding, only: c_loc, c_f_pointer
    use,intrinsic :: iso_fortran_env, only: iostat_end
    use string_utilities, only: i_to_c

    type(toolpath), allocatable, intent(out) :: path
    type(parameter_list) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(path_builder), target :: builder
    type(fyajl_parser), target :: parser
    type(fyajl_status) :: yajl_stat
    character, pointer :: buffer(:)
    character(:), allocatable, target :: string
    integer :: lun, ios, last_pos, curr_pos, buflen

    stat = 0

    call builder%init(params)
    call parser%init(builder)

    if (params%is_parameter('command-string')) then ! parse string

      call params%get('command-string', string)
      call c_f_pointer(c_loc(string(1:1)), buffer, shape=[len(string)]) ! possibly dicey
      call parser%parse(buffer, yajl_stat)
      if (yajl_stat /= FYAJL_STATUS_OK) then
        errmsg = fyajl_get_error(parser, .true., buffer)
        if (yajl_stat == FYAJL_STATUS_CLIENT_CANCELED) errmsg = errmsg // builder%errmsg
        stat = 1
        return
      end if

      call parser%complete_parse(yajl_stat)
      if (yajl_stat /= FYAJL_STATUS_OK) then
        errmsg = fyajl_get_error(parser, .false., buffer)
        if (yajl_stat == FYAJL_STATUS_CLIENT_CANCELED) errmsg = errmsg // builder%errmsg
        stat = 1
        return
      else
        call move_alloc(builder%path, path)
      end if

    else  ! parse file

      call params%get('command-file', string)
      open(newunit=lun,file=string,action='read',access='stream',form='unformatted',iostat=ios)
      if (ios /= 0) then
        errmsg = 'unable to open file "' // string // '": iostat=' // i_to_c(ios)
        stat = -1
        return
      end if
      call parser%set_option(FYAJL_ALLOW_COMMENTS)
      allocate(buffer(1000))
      inquire(lun,pos=last_pos)  ! starting position in stream
      do
        !! Fill buffer a character at a time
        do buflen = 1, size(buffer)
          read(lun,iostat=ios) buffer(buflen)
          if (ios /= 0) exit
        end do
        if (ios /= 0 .and. ios /= iostat_end) then
          errmsg = 'read error: iostat=' // i_to_c(ios)
          stat = -1
          exit
        end if

        inquire(lun,pos=curr_pos)
        buflen = curr_pos - last_pos
        last_pos = curr_pos
        if (buflen > 0) then
          call parser%parse(buffer(:buflen), yajl_stat)
          if (yajl_stat /= FYAJL_STATUS_OK) then
            errmsg = fyajl_get_error(parser, .true., buffer(:buflen))
            if (yajl_stat == FYAJL_STATUS_CLIENT_CANCELED) errmsg = errmsg // builder%errmsg
            stat = 1
            exit
          end if
        end if

        if (ios == iostat_end) then
          call parser%complete_parse(yajl_stat)
          if (yajl_stat /= FYAJL_STATUS_OK) then
            errmsg = fyajl_get_error(parser, .false., buffer(:buflen))
            if (yajl_stat == FYAJL_STATUS_CLIENT_CANCELED) errmsg = errmsg // builder%errmsg
            stat = 1
          else
            call move_alloc(builder%path, path)
          end if
          exit
        end if
      end do
      deallocate(buffer)

    end if

  end subroutine alloc_toolpath

  !! Initialize the PATH_BUILDER object
  subroutine init(this, params)

    class(path_builder), intent(out) :: this
    type(parameter_list) :: params

    real(r8), allocatable :: r(:)

    call params%get('start-time',  this%t, default=0.0_r8)
    call params%get('start-coord', r, default=[0.0_r8, 0.0_r8, 0.0_r8])
    call params%get('time-scale-factor',  this%tsf, default=1.0_r8)
    call params%get('coord-scale-factor', this%csf, default=1.0_r8)

    this%t = this%tsf * this%t
    this%r = this%csf * r
    this%flags = 0  ! all start cleared

    allocate(this%path)
    call this%path%append_path_segment(new_start_segment(this%t, this%r, this%flags))

    this%state = STATE_INIT

  end subroutine init

!!!! YAJL_FORT PARSER CALLBACK FUNCTIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function start_array(this) result(status)
    class(path_builder) :: this
    select case (this%state)
    case (STATE_INIT)
      this%state = STATE_COMMAND
      status = FYAJL_CONTINUE_PARSING
    case (STATE_COMMAND)
      this%state = STATE_ACTION
      status = FYAJL_CONTINUE_PARSING
    case (STATE_DR)
      this%state = STATE_DX
      status = FYAJL_CONTINUE_PARSING
    case default
      this%errmsg = 'invalid syntax: array not expected here'
      status = FYAJL_TERMINATE_PARSING
    end select
  end function start_array

  integer function end_array(this) result(status)
    use xyz_motion_class
    use dwell_xyz_motion_type
    use linear_xyz_motion_type
    class(path_builder) :: this
    class(xyz_motion), allocatable :: move
    select case (this%state)
    case (STATE_COMMAND)
      call this%path%append_path_segment(new_final_segment(this%t, this%r, this%flags))
      this%state = STATE_DONE
      status = FYAJL_CONTINUE_PARSING
    case (STATE_SETFLAG2, STATE_CLRFLAG2)
      this%state = STATE_COMMAND
      status = FYAJL_CONTINUE_PARSING
    case (STATE_DWELL_DONE)
#ifdef NO_2008_LHS_POLY_REALLOC
      allocate(move, source=dwell_xyz_motion(this%r, t0=this%t, dt=this%dt))
#else
      move = dwell_xyz_motion(this%r, t0=this%t, dt=this%dt)
#endif
      this%t = move%final_time(); this%r = move%final_coord()
      call this%path%append_path_segment(new_path_segment(move, this%flags))
      this%state = STATE_COMMAND
      status = FYAJL_CONTINUE_PARSING
    case (STATE_DR_DONE)
      this%state = STATE_SPEED
      status = FYAJL_CONTINUE_PARSING
    case (STATE_ACCEL)
#ifdef NO_2008_LHS_POLY_REALLOC
      allocate(move, source=linear_xyz_motion(this%t, this%r, this%dr, this%speed))
#else
      move = linear_xyz_motion(this%t, this%r, this%dr, this%speed)
#endif
      this%t = move%final_time(); this%r = move%final_coord()
      call this%path%append_path_segment(new_path_segment(move, this%flags))
      this%state = STATE_COMMAND
      status = FYAJL_CONTINUE_PARSING
    case (STATE_DECEL)
#ifdef NO_2008_LHS_POLY_REALLOC
      allocate(move, source=linear_xyz_motion(this%t, this%r, this%dr, this%speed, this%accel))
#else
      move = linear_xyz_motion(this%t, this%r, this%dr, this%speed, this%accel)
#endif
      this%t = move%final_time(); this%r = move%final_coord()
      call this%path%append_path_segment(new_path_segment(move, this%flags))
      this%state = STATE_COMMAND
      status = FYAJL_CONTINUE_PARSING
    case (STATE_MOVEREL_DONE)
#ifdef NO_2008_LHS_POLY_REALLOC
      allocate(move, source=linear_xyz_motion(this%t, this%r, this%dr, this%speed, this%accel, this%decel))
#else
      move = linear_xyz_motion(this%t, this%r, this%dr, this%speed, this%accel, this%decel)
#endif
      this%t = move%final_time(); this%r = move%final_coord()
      call this%path%append_path_segment(new_path_segment(move, this%flags))
      this%state = STATE_COMMAND
      status = FYAJL_CONTINUE_PARSING
    case default
      this%errmsg = 'invalid syntax: premature end to array'
      status = FYAJL_TERMINATE_PARSING
    end select
  end function end_array

  integer function integer_value(this, value) result(status)
    class(path_builder) :: this
    integer(fyajl_integer_kind), intent(in) :: value
    select case (this%state)
    case (STATE_SETFLAG1,STATE_SETFLAG2)
      if (value >= 0 .and. value < bit_size(this%flags)) then
        this%flags = ibset(this%flags,pos=value)
        this%state = STATE_SETFLAG2
        status = FYAJL_CONTINUE_PARSING
      else
        this%errmsg = 'invalid flag number'
        status = FYAJL_TERMINATE_PARSING
      end if
    case (STATE_CLRFLAG1,STATE_CLRFLAG2)
      if (value >= 0 .and. value < bit_size(this%flags)) then
        this%flags = ibclr(this%flags,pos=value)
        this%state = STATE_CLRFLAG2
        status = FYAJL_CONTINUE_PARSING
      else
        this%errmsg = 'invalid flag number'
        status = FYAJL_TERMINATE_PARSING
      end if
    case (STATE_DWELL_DT)
      if (value > 0) then
        this%dt = this%tsf * value
        this%state = STATE_DWELL_DONE
        status = FYAJL_CONTINUE_PARSING
      else
        this%errmsg = 'positive value expected'
        status = FYAJL_TERMINATE_PARSING
      end if
    case (STATE_DX)
      this%dr(1) = this%csf * value
      this%state = STATE_DY
      status = FYAJL_CONTINUE_PARSING
    case (STATE_DY)
      this%dr(2) = this%csf * value
      this%state = STATE_DZ
      status = FYAJL_CONTINUE_PARSING
    case (STATE_DZ)
      this%dr(3) = this%csf * value
      this%state = STATE_DR_DONE
      status = FYAJL_CONTINUE_PARSING
    case (STATE_SPEED)
      if (value > 0) then
        this%speed = (this%csf/this%tsf) * value
        this%state = STATE_ACCEL
        status = FYAJL_CONTINUE_PARSING
      else
        this%errmsg = 'positive value expected'
        status = FYAJL_TERMINATE_PARSING
      end if
    case (STATE_ACCEL)
      if (value > 0) then
        this%accel = (this%csf/this%tsf**2) * value
        this%state = STATE_DECEL
        status = FYAJL_CONTINUE_PARSING
      else
        this%errmsg = 'positive value expected'
        status = FYAJL_TERMINATE_PARSING
      end if
    case (STATE_DECEL)
      if (value > 0) then
        this%decel = (this%csf/this%tsf**2) * value
        this%state = STATE_MOVEREL_DONE
        status = FYAJL_CONTINUE_PARSING
      else
        this%errmsg = 'positive value expected'
        status = FYAJL_TERMINATE_PARSING
      end if
    case default
      this%errmsg = 'invalid syntax: number not expected here'
      status = FYAJL_TERMINATE_PARSING
    end select
  end function integer_value

  integer function double_value(this, value) result(status)
    class(path_builder) :: this
    real(fyajl_real_kind), intent(in) :: value
    select case (this%state)
    case (STATE_DWELL_DT)
      if (value > 0) then
        this%dt = this%tsf * value
        this%state = STATE_DWELL_DONE
        status = FYAJL_CONTINUE_PARSING
      else
        this%errmsg = 'positive value expected'
        status = FYAJL_TERMINATE_PARSING
      end if
    case (STATE_DX)
      this%dr(1) = this%csf * value
      this%state = STATE_DY
      status = FYAJL_CONTINUE_PARSING
    case (STATE_DY)
      this%dr(2) = this%csf * value
      this%state = STATE_DZ
      status = FYAJL_CONTINUE_PARSING
    case (STATE_DZ)
      this%dr(3) = this%csf * value
      this%state = STATE_DR_DONE
      status = FYAJL_CONTINUE_PARSING
    case (STATE_SPEED)
      if (value > 0) then
        this%speed = (this%csf/this%tsf) * value
        this%state = STATE_ACCEL
        status = FYAJL_CONTINUE_PARSING
      else
        this%errmsg = 'positive value expected'
        status = FYAJL_TERMINATE_PARSING
      end if
    case (STATE_ACCEL)
      if (value > 0) then
        this%accel = (this%csf/this%tsf**2) * value
        this%state = STATE_DECEL
        status = FYAJL_CONTINUE_PARSING
      else
        this%errmsg = 'positive value expected'
        status = FYAJL_TERMINATE_PARSING
      end if
    case (STATE_DECEL)
      if (value > 0) then
        this%decel = (this%csf/this%tsf**2) * value
        this%state = STATE_MOVEREL_DONE
        status = FYAJL_CONTINUE_PARSING
      else
        this%errmsg = 'positive value expected'
        status = FYAJL_TERMINATE_PARSING
      end if
    case default
      this%errmsg = 'invalid syntax: floating number not expected here'
      status = FYAJL_TERMINATE_PARSING
    end select
  end function double_value

  integer function string_value(this, value) result(status)
    class(path_builder) :: this
    character(*), intent(in) :: value
    select case (this%state)
    case (STATE_ACTION)
      select case (value)
      case ('setflag')
        this%state = STATE_SETFLAG1
        status = FYAJL_CONTINUE_PARSING
      case ('clrflag')
        this%state = STATE_CLRFLAG1
        status = FYAJL_CONTINUE_PARSING
      case ('dwell')
        this%state = STATE_DWELL_DT
        status = FYAJL_CONTINUE_PARSING
      case ('moverel')
        this%state = STATE_DR
        status = FYAJL_CONTINUE_PARSING
      case default
        this%errmsg = 'invalid action "' // value // '"'
        status = FYAJL_TERMINATE_PARSING
      end select
    case default
      this%errmsg = 'invalid syntax: string not expected here'
      status = FYAJL_TERMINATE_PARSING
    end select
  end function string_value

  !! The remaining parsing events are not part of the restricted path language.

  integer function start_map(this) result(status)
    class(path_builder) :: this
    this%errmsg = 'invalid syntax'
    status = FYAJL_TERMINATE_PARSING
  end function start_map

  integer function end_map(this) result(status)
    class(path_builder) :: this
    this%errmsg = 'invalid syntax'
    status = FYAJL_TERMINATE_PARSING
  end function end_map

  integer function map_key(this, value) result(status)
    class(path_builder) :: this
    character(*), intent(in) :: value
    this%errmsg = 'invalid syntax'
    status = FYAJL_TERMINATE_PARSING
  end function map_key

  integer function null_value(this) result(status)
    class(path_builder) :: this
    this%errmsg = 'invalid syntax'
    status = FYAJL_TERMINATE_PARSING
  end function null_value

  integer function logical_value(this, value) result(status)
    class(path_builder) :: this
    logical, intent(in) :: value
    this%errmsg = 'invalid syntax'
    status = FYAJL_TERMINATE_PARSING
  end function logical_value

end module toolpath_factory
