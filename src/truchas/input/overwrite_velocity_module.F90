!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module overwrite_velocity_module
  !=======================================================================
  ! Purpose(s):
  !   Define procedures for the input of requesting velocity overwriting
  !
  ! Public Interface(s):
  !
  !   * call overwrite_velocity_input()
  !
  !     Read and process the INTERFACES namelist.
  !
  ! Contains: overwrite_velocity_default
  !           overwrite_velocity_input
  !           overwrite_velocity_parallel
  !
  ! Author(s): Robert Chiodi (robertchiodi@lanl.gov)
  !
  !=======================================================================
  use truchas_logging_services
  implicit none
  private

  public :: overwrite_velocity_input

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

contains

  subroutine overwrite_velocity_default ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Default OVERWRITE_VELOCITY  namelist.
    !
    !=======================================================================
    use vof_velocity_overwrite, only: velocity_overwrite_requested, &
                                      velocity_overwrite_case

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    velocity_overwrite_requested = .false.
    velocity_overwrite_case = 'NULL'


  end subroutine overwrite_velocity_default

  !-----------------------------------------------------------------------------

  subroutine overwrite_velocity_input (lun)
    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !   Read OVERWRITE_VELOCITY namelist. Sets use of overwriting velocity
    !   in vof_velocity_overwrite module and records the case name to be used
    !   later to set a known velocity.
    !
    !---------------------------------------------------------------------------
    use vof_velocity_overwrite, only: velocity_overwrite_requested, &
                                      velocity_overwrite_case
    use input_utilities,        only: seek_to_namelist
    use parallel_info_module,   only: p_info
    use parameter_module,       only: mbody
    use pgslib_module,          only: PGSLib_GLOBAL_ANY, pgslib_bcast
    use string_utilities,       only: lower_case
    use string_utilities, only: i_to_c

    integer, intent(in) :: lun

    ! local variables
    logical :: fatal, found
    integer :: ios, ib

    ! overwrite_velocity namelist
    namelist /OVERWRITE_VELOCITY/ velocity_overwrite_case

    !---------------------------------------------------------------------------

    ! assume success
    fatal = .false.

    call TLS_info ('')
    call TLS_info ('Reading OVERWRITE_VELOCITY namelist ...')

    ! prepare to read
    call overwrite_velocity_default  ()

    !! Locate the INTERFACES namelist (first occurence).
    found = .false.
    if (p_info%IOP) then
      rewind lun
      call seek_to_namelist (lun, 'OVERWRITE_VELOCITY', found, iostat=ios)
    end if
    if(found) then
       velocity_overwrite_requested = .true.
    else
       return
    end if
    
    call pgslib_bcast (ios)
    if (ios /= 0) call TLS_fatal ('error reading input file; iostat=' // i_to_c(ios))

    call pgslib_bcast (found)
    if (found) then
      if (p_info%IOP) read(lun, nml=overwrite_velocity, iostat=ios)
      call pgslib_bcast (ios)
      if (ios /= 0) call TLS_fatal ('error reading OVERWRITE_VELOCITY namelist; iostat=' // i_to_c(ios))
    else
      call TLS_info ('  OVERWRITE_VELOCITY  namelist not found; not overwriting velocities')
    end if

    ! broadcast data
    call overwrite_velocity_input_parallel ()

    call TLS_info ('Velocity will be overwritten for case "' // trim(velocity_overwrite_case) // '"')
    
  end subroutine overwrite_velocity_input

  !-----------------------------------------------------------------------------

  subroutine overwrite_velocity_input_parallel ()
    !======================================================================
    ! Purpose(s):
    !
    !   Broadcast components of the overwrite_velocity namelist to all PE's.
    !
    !======================================================================
    use vof_velocity_overwrite, only: velocity_overwrite_requested, &
                                      velocity_overwrite_case
    use parallel_info_module, only: p_info
    use pgslib_module,        only: PGSLib_BCAST

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Broadcast Data
    if (.not. p_info%UseGlobalServices) then
       call PGSLib_BCAST (velocity_overwrite_requested)
       call PGSLib_BCAST (velocity_overwrite_case)
    endif

  end subroutine overwrite_velocity_input_parallel

end module overwrite_velocity_module
