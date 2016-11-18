!!
!! TOOLPATH_NAMELIST
!!
!! This module provides a procedure for reading the TOOLPATH namelists and
!! populating the toolpath table with the toolpaths they define.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! November 2016
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  CALL READ_TOOLPATH_NAMELISTS (LUN) reads all the so-named namelists,
!!    creates the specified toolpath objects, and stores them in the global
!!    toolpath table.  See the TOOLPATH_TABLE module for subsequent access.
!!

#include "f90_assert.fpp"

module toolpath_namelist

  implicit none
  private

  public :: read_toolpath_namelists

contains

  subroutine read_toolpath_namelists(lun)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use parallel_communication, only: is_IOP, broadcast
    use string_utilities, only: i_to_c
    use input_utilities, only: seek_to_namelist, NULL_R, NULL_C
    use truchas_env, only: input_dir
    use toolpath_table, only: known_toolpath
    use truchas_logging_services

    integer, intent(in) :: lun

    logical :: found
    integer :: n, ios

    !! namelist variables
    character(31) :: name
    character(256) :: command_file
    character(1000) :: command_string
    real(r8) :: start_time, start_coord(3), time_scale_factor, coord_scale_factor

    namelist /toolpath/ name, start_time, start_coord, time_scale_factor, coord_scale_factor, &
                        command_string, command_file

    call TLS_info('')
    call TLS_info('Reading TOOLPATH namelists ...')

    if (is_IOP) rewind(lun)
    n = 0

    do  ! until all TOOLPATH namelists have been read or an error occurs

      if (is_IOP) call seek_to_namelist(lun, 'TOOLPATH', found, iostat=ios)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading input file')

      call broadcast(found)
      if (.not.found) return  ! no further TOOLPATH namelists found

      n = n + 1
      call TLS_info('  Reading TOOLPATH namelist #' // i_to_c(n))

      !! Read the namelist variables, assigning default values first.
      if (is_IOP) then
        name = NULL_C
        start_time = NULL_R
        start_coord = NULL_R
        time_scale_factor = NULL_R
        coord_scale_factor = NULL_R
        command_string = NULL_C
        command_file = NULL_C
        read(lun,nml=toolpath,iostat=ios)
      end if
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading TOOLPATH namelist')

      !! Broadcast the namelist variables.
      call broadcast(name)
      call broadcast(start_time)
      call broadcast(start_coord)
      call broadcast(time_scale_factor)
      call broadcast(coord_scale_factor)
      call broadcast(command_string)
      call broadcast(command_file)

      !! Check the name.
      if (name == NULL_C .or. name == '') &
          call TLS_fatal('NAME must be assigned a nonempty value')
      if (known_toolpath(name)) &
          call TLS_fatal('already read a TOOLPATH namelist with this name: ' // trim(name))

      !! Check that we got either COMMAND_STRING or COMMAND_FILE, but not both.
      if (command_string == NULL_C .and. command_file == NULL_C) &
          call TLS_fatal('either COMMAND_STRING or COMMAND_FILE must be specified')
      if (command_string /= NULL_C .and. command_file /= NULL_C) &
          call TLS_fatal('cannot specify both COMMAND_STRING and COMMAND_FILE')

      !! Adjust the file path and check it is accessible.
      if (command_file /= NULL_C) then
        if (command_file(1:1) /= '/') then ! not an absolute path
          command_file = trim(input_dir) // trim(command_file) ! make relative to input dir
        end if
        if (is_IOP) inquire(file=trim(command_file),exist=found)
        call broadcast(found)
        if (.not.found) call TLS_fatal('COMMAND_FILE not found: ' // trim(command_file))
      end if

      !! Check scaling factors if specified
      if (time_scale_factor /= NULL_R) then
        if (time_scale_factor <= 0) call TLS_fatal('TIME_SCALE_FACTOR must be > 0')
      end if
      if (coord_scale_factor /= NULL_R) then
        if (coord_scale_factor <= 0) call TLS_fatal('COORD_SCALE_FACTOR must be > 0')
      end if

      !! Check START_COORD is complete if specified
      if (any(start_coord /= NULL_R)) then
        if (any(start_coord == NULL_R)) call TLS_fatal('START_COORD not completely specified')
      end if

      call make_toolpath

    end do

  contains

    subroutine make_toolpath

      use parameter_list_type
      use toolpath_factory
      use toolpath_table, only: insert_toolpath

      type(parameter_list) :: params
      type(toolpath), allocatable :: path
      integer :: stat
      character(:), allocatable :: errmsg

      if (start_time /= NULL_R) call params%set('start-time', start_time)
      if (any(start_coord /= NULL_R)) call params%set('start-coord', start_coord)
      if (time_scale_factor /= NULL_R) call params%set('time-scale-factor', time_scale_factor)
      if (coord_scale_factor /= NULL_R) call params%set('coord-scale-factor', coord_scale_factor)
      if (command_string /= NULL_C) then
        call params%set('command-string', trim(command_string))
      else
        call params%set('command-file', trim(command_file))
      end if
      call alloc_toolpath(path, params, stat, errmsg)
      if (stat /= 0) call TLS_fatal('error creating toolpath: ' // errmsg)
      call insert_toolpath(trim(name), path)

    end subroutine make_toolpath

  end subroutine read_toolpath_namelists

end module toolpath_namelist
