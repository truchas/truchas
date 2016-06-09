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
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module ded_head_driver

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use ded_head_namelist, only: read_ded_head_namelist, ded_params
  use ded_head_type
  use scalar_func_class
  implicit none
  private
  
  public :: read_ded_head_namelist
  public :: ded_head_init
  public :: ded_head_start_sim_phase
  public :: alloc_ded_head_laser_func

  type(ded_head), allocatable, public :: head
  
  !! SCALAR_FUNC class adapter for the laser source.  This enables the DED
  !! laser source to be used as a function for heat flux boundary conditions.
  !! NB: Function assumes the outward normal to the boundary is z-directed.
  type, extends(scalar_func) :: laser_scalar_func
  contains
    procedure :: eval => laser_func_eval
  end type laser_scalar_func

contains

  subroutine ded_head_init(t)

    use simulation_event_queue

    real(r8), intent(in) :: t

    integer :: j
    real(r8), allocatable :: times(:)
    logical, allocatable :: discont(:)

    if (.not.associated(ded_params)) return ! move along, nothing to see here
    
    !! Instantiate the DED_HEAD module object using the input parameters.
    allocate(head)
    call head%init(ded_params)

    !! Position the toolpath at the initial simulation time.
    call head%tp%set_segment(t)

    !! Add toolpath segment endpoints to the list of simulation events.
    call head%tp%get_segment_starts(times, discont)
    do j = 1, size(times)
      if (discont(j)) then
        call add_event(sim_event(times(j), DT_POLICY_FACTOR, 1.0e-3_r8))
      else
        call add_event(sim_event(times(j), DT_POLICY_NEXT))
      end if
    end do

  end subroutine ded_head_init

  !! Called at the start of each simulation phase, which are assumed here to
  !! correspond to segments of the toolpath, to advance the toolpath to the
  !! next segment.  FIXME! THIS IS INCORRECT IF OTHER SIMULATION PHASES ARE
  !! DEFINED BY THE SIMULATION_CONTROL NAMELIST.  NEED TO REGISTER THE ACTIONS
  !! TO TAKE WITH EACH SIMULATION EVENT.

  subroutine ded_head_start_sim_phase
    if (allocated(head)) call head%tp%next_segment()
  end subroutine ded_head_start_sim_phase

!!!! SCALAR_FUNC ADAPTER FOR THE LASER SOURCE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine alloc_ded_head_laser_func(f)
    class(scalar_func), allocatable, intent(out) :: f
    allocate(f, source=laser_scalar_func())
  end subroutine alloc_ded_head_laser_func
  
  !! Type bound EVAL procedure; delegates to the modules DED_HEAD object.
  function laser_func_eval(this, x) result(fx)
    class(laser_scalar_func), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: fx
    fx = -head%laser_irrad(x(1), x(2:4))
  end function laser_func_eval
  
end module ded_head_driver
