!!
!! SIMULATION_FACTORY
!!
!! This module constructs a two-dimensional simulation from its command-line
!! name.  The returned object has the common simulation lifecycle interface;
!! concrete simulations retain ownership of their state and implementation.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, September 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

module simulation_factory

  use simulation_type, only: simulation
  use ns_2d_sim_type, only: ns_2d_sim
  use ht_2d_sim_type, only: ht_2d_sim
  use ns_ht_2d_sim_type, only: ns_ht_2d_sim
  use vof_2d_sim_type, only: vof_2d_sim
  implicit none
  private

  public :: new_simulation

contains

  subroutine new_simulation(name, sim, stat, errmsg)

    character(*), intent(in) :: name
    class(simulation), allocatable, intent(out) :: sim
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    stat = 0
    select case (trim(name))
    case ('flow')
      allocate(ns_2d_sim :: sim)
    case ('thermal')
      allocate(ht_2d_sim :: sim)
    case ('flow_thermal')
      allocate(ns_ht_2d_sim :: sim)
    case ('vof')
      allocate(vof_2d_sim :: sim)
    case default
      stat = 1
      errmsg = 'unknown simulation: ' // trim(name)
    end select

  end subroutine new_simulation

end module simulation_factory
