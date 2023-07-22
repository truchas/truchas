!!
!! TOOLHEAD_DRIVER
!!
!! Driver for the optional model describing the toolhead of an additive
!! manufacturing process. Conceptually this is part of the Truchas physics
!! driver (driver.F90), serving as an adapter to the toolhead object, which
!! is held here as a module variable.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! May 2016
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module toolhead_driver

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use toolhead_table_type
  use vector_func_class
  use sim_event_queue_type
  implicit none
  private

  public :: toolhead_init
  public :: add_toolhead_events

  type(toolhead_table), public :: table

contains

  subroutine toolhead_init(t)
    real(r8), intent(in) :: t
    call table%set_initial_state(t)
  end subroutine


  subroutine add_toolhead_events(eventq)
    use time_step_module, only: dt_init
    type(sim_event_queue), intent(inout) :: eventq
    call table%add_events(dt_init, eventq)
  end subroutine

end module toolhead_driver
