!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module am_coord_system_factory

  use am_coord_system_class
  use am_coord_system_block_type
  use parameter_list_type
  use truchas_logging_services
  implicit none
  private

  public :: alloc_am_coord_system

  contains


    subroutine alloc_am_coord_system(this, params)

      class(am_coord_system), allocatable, intent(out) :: this
      type(parameter_list), target :: params

      integer :: stat
      character(:), allocatable :: am_coord_system_type, errmsg

      call params%get ('am_coord_system_type', am_coord_system_type, &
                       stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal ('ALLOCATE_COORD_SYSTEM: ' // errmsg)

      if (am_coord_system_type == "block_3D") then
        allocate (am_coord_system_block :: this)
      else
        call TLS_fatal ('ALLOCATE_COORD_SYSTEM: unknown "am_coord_system_type": ' //  am_coord_system_type)
      end if

    end subroutine alloc_am_coord_system


end module am_coord_system_factory
