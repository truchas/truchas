!!
!! DED_HEAD_DRIVER
!!
!! Driver for the optional model describing the tool head of a LENS-like
!! directed energy deposition (DED) machine.  Conceptually this is part
!! of the Truchas physics driver (driver.F90), serving as an adapter to
!! the DED head object, which is held here as a module variable.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! May 2016
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module ded_head_driver

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use ded_head_namelist, only: read_ded_head_namelist, ded_params
  use ded_head_type
  use vector_func_class
  use sim_event_queue_type
  implicit none
  private

  public :: read_ded_head_namelist
  public :: ded_head_init
  public :: ded_head_path_events
  public :: ded_head_start_sim_phase
  public :: alloc_ded_head_laser_func

  type(ded_head), allocatable, public :: head

  !! VECTOR_FUNC class adapter for the laser source.  This enables the DED
  !! laser source to be used as a function for heat flux boundary conditions.
  type, extends(vector_func) :: laser_vector_func
  contains
    procedure :: eval => laser_func_eval
    procedure :: eval_comp => laser_func_eval_comp
  end type laser_vector_func

  integer, parameter, public :: DT_POLICY_NONE   = 0
  integer, parameter, public :: DT_POLICY_NEXT   = 1
  integer, parameter, public :: DT_POLICY_FACTOR = 2
  integer, parameter, public :: DT_POLICY_VALUE  = 3

  type, extends(event_action), public :: path_event
    private
    integer  :: dt_policy = DT_POLICY_NONE
    real(r8) :: c = 0.0_r8
  contains
    procedure :: init_dt => path_event_init_dt
  end type path_event

contains

  subroutine ded_head_init(t)

    real(r8), intent(in) :: t

    if (.not.associated(ded_params)) return ! move along, nothing to see here

    !! Instantiate the DED_HEAD module object using the input parameters,
    !! positioning the toolpath at the initial simulation time.
    allocate(head)
    call head%init(ded_params, t)

  end subroutine ded_head_init


  subroutine ded_head_path_events(eventq)

    type(sim_event_queue), intent(inout) :: eventq

    integer :: j
    real(r8), allocatable :: times(:)
    logical, allocatable :: discont(:)

    if (allocated(head)) then
      call head%tp%get_segment_starts(times, discont)
      do j = 1, size(times)
        if (discont(j)) then
          call eventq%add_event(times(j), path_event(DT_POLICY_FACTOR, 1.0e-3_r8))
        else
          call eventq%add_event(times(j), path_event(DT_POLICY_NEXT))
        end if
      end do
    end if

  end subroutine ded_head_path_events

  !! Called at the start of each simulation phase, which are assumed here to
  !! correspond to segments of the toolpath, to advance the toolpath to the
  !! next segment.  FIXME! THIS IS INCORRECT IF OTHER SIMULATION PHASES ARE
  !! DEFINED BY THE SIMULATION_CONTROL NAMELIST.  NEED TO REGISTER THE ACTIONS
  !! TO TAKE WITH EACH SIMULATION EVENT.

  subroutine ded_head_start_sim_phase(t)
    real(r8), intent(in) :: t
    if (allocated(head)) call head%next_tp_segment(t)
  end subroutine ded_head_start_sim_phase

!!!! VECTOR_FUNC ADAPTER FOR THE LASER SOURCE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine alloc_ded_head_laser_func(f)
    class(vector_func), allocatable, intent(out) :: f
    allocate(f, source=laser_vector_func(dim=3))
  end subroutine alloc_ded_head_laser_func

  !! Type bound EVAL procedure; delegates to the modules DED_HEAD object.
  function laser_func_eval(this, x) result(fx)
    class(laser_vector_func), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: fx(this%dim)
    integer :: i
    fx = 0
    fx(3) = -head%laser_irrad(x(1), x(2:4))
  end function laser_func_eval

  function laser_func_eval_comp(this, i, x) result(fx)
    class(laser_vector_func), intent(in) :: this
    integer, intent(in) :: i
    real(r8), intent(in) :: x(:)
    real(r8) :: fx
    real(r8) :: f(this%dim)
    f = this%eval(x)
    fx = f(i)
  end function laser_func_eval_comp

  pure function path_event_init_dt(this, dt_last, dt_next) result(dt)
    class(path_event), intent(in) :: this
    real(r8), intent(in) :: dt_last, dt_next
    real(r8) :: dt
    select case (this%dt_policy)
    case (DT_POLICY_NONE, DT_POLICY_NEXT)
      dt = dt_next
    case (DT_POLICY_FACTOR)
      dt = this%c*dt_last
    case (DT_POLICY_VALUE)
      dt = this%c
    case default
      dt = 0.0_r8 ! should never be here
    end select
  end function path_event_init_dt

end module ded_head_driver
