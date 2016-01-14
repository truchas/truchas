!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module interface_energy_factory

  use interface_energy_class
  use interface_energy_gaussian_type
  use parameter_list_type
  use truchas_logging_services
  implicit none
  private

  public :: alloc_interface_energy

  contains


    subroutine alloc_interface_energy(this, params)

      class(interface_energy), allocatable, intent(out) :: this
      type(parameter_list), target :: params

      integer :: stat
      character(:), allocatable :: interface_energy_type, errmsg

      call params%get ('interface_energy_type', interface_energy_type, &
                       stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal ('ALLOCATE_INTERFACE_ENERGY: ' // errmsg)

      if (interface_energy_type == "gaussian") then
        allocate (interface_energy_gaussian :: this)
        call this%alloc()
      else
        call TLS_fatal ('ALLOCATE_INTERFACE_ENERGY: unknown "interface_energy_type": ' //  interface_energy_type)
      end if

    end subroutine alloc_interface_energy


end module interface_energy_factory
