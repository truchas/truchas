!!
!! SIMULATION_TYPE
!!
!! This module defines the common lifecycle interface for self-contained
!! simulations.  Concrete simulations own their state and implementation;
!! this type only establishes the interface used by a generic driver.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, September 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

module simulation_type

  use parameter_list_type, only: parameter_list
  use simulation_environment_type, only: simulation_environment
  implicit none
  private

  type, abstract, public :: simulation
  contains
    procedure(simulation_init_interface), deferred :: init
    procedure(simulation_run_interface), deferred :: run
  end type simulation

  abstract interface
    subroutine simulation_init_interface(this, env, params, stat, errmsg)
      import simulation, simulation_environment, parameter_list
      class(simulation), intent(out) :: this
      type(simulation_environment), intent(inout) :: env
      type(parameter_list), intent(inout) :: params
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
    end subroutine simulation_init_interface

    subroutine simulation_run_interface(this, env, stat, errmsg)
      import simulation, simulation_environment
      class(simulation), intent(inout) :: this
      type(simulation_environment), intent(inout) :: env
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
    end subroutine simulation_run_interface
  end interface

end module simulation_type
