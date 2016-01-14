!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module interface_mass_factory

  use interface_mass_class
  use interface_mass_gaussian_type
  use parameter_list_type
  use truchas_logging_services
  implicit none
  private

  public :: alloc_interface_mass

  contains


    subroutine alloc_interface_mass(this, params)

      class(interface_mass), allocatable, intent(out) :: this
      type(parameter_list), target :: params

      integer :: stat
      character(:), allocatable :: interface_mass_type, errmsg


      call params%get ('interface_mass_type', interface_mass_type, &
                       stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal ('ALLOCATE_INTERFACE_MASS: ' // errmsg)

      if (interface_mass_type == "gaussian") then
        allocate (interface_mass_gaussian :: this)
        call this%alloc()
      else
        call TLS_fatal ('ALLOCATE_INTERFACE_MASS: unknown "interface_mass_type": ' //  interface_mass_type)
      end if

    end subroutine alloc_interface_mass


end module interface_mass_factory
