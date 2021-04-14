!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module advection_velocity_namelist

  use vector_func_class
  implicit none
  private

  public :: read_advection_velocity_namelist

  class(vector_func), allocatable, public :: adv_vel

contains

  subroutine read_advection_velocity_namelist(lun)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use parallel_communication, only: is_IOP, broadcast
    use string_utilities, only: i_to_c
    use input_utilities, only: seek_to_namelist, NULL_R, NULL_C
    use vector_func_factories
    use truchas_logging_services

    integer, intent(in) :: lun

    logical :: found
    integer :: ios
    character(128) :: iom

    real(r8) :: velocity_constant(3)
    namelist /advection_velocity/ velocity_constant

    call TLS_info('Reading ADVECTION_VELOCITY namelist ...')

    if (is_IOP) then
      rewind(lun)
      call seek_to_namelist(lun, 'ADVECTION_VELOCITY', found, iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))

    call broadcast(found)
    if (.not.found) call TLS_fatal('ADVECTION_VELOCITY namelist not found')

    if (is_IOP) then
      velocity_constant = NULL_R
      read(lun,nml=advection_velocity,iostat=ios,iomsg=iom)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading ADVECTION_VELOCITY namelist: ' // trim(iom))

    call broadcast(velocity_constant)

    if (all(velocity_constant == NULL_R)) call TLS_fatal('ADVECTION_VELOCITY not defined')
    if (any(velocity_constant == NULL_R)) call TLS_fatal('ADVECTION_VELOCITY not fully defined')
    call alloc_const_vector_func(adv_vel, velocity_constant)

  end subroutine read_advection_velocity_namelist

end module advection_velocity_namelist
