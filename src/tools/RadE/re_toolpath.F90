!!
!! RE_TOOLPATH
!!
!! This module provides a procedure for reading the TOOLPATH namelists and
!! populating the toolpath table with the toolpaths they define.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! December 2016
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NB: This is a slightly modified version of TOOLPATH_NAMELIST from Truchas.
!! It is poor software practice, and I hated doing it, but I did not see a way
!! of avoiding it.  Genre uses a different MPI abstraction, and the Truchas
!! version relies on a few Truchas-specific utilities (message logging, file)
!! that are not usable outside of the Truchas runtime.
!!
!! PROGRAMMING INTERFACE
!!
!!  CALL READ_TOOLPATH_NAMELISTS (LUN) reads all the so-named namelists,
!!    creates the specified toolpath objects, and stores them in the global
!!    toolpath table.  See the TOOLPATH_TABLE module for subsequent access.
!!

#include "f90_assert.fpp"

module re_toolpath

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use re_utilities
  use toolpath_factory_type
  use scl
  implicit none
  private

  public :: read_toolpath_namelists

  type(toolpath_factory), public :: tp_fac

contains

  subroutine read_toolpath_namelists(lun)

    use string_utilities, only: i_to_c
    use input_utilities, only: seek_to_namelist, NULL_R, NULL_C

    integer, intent(in) :: lun

    logical :: is_IOP, found
    integer :: n, ios

    !! namelist variables
    character(31) :: name
    character(256) :: command_file, iom
    character(1000) :: command_string
    real(r8) :: start_time, start_coord(3), time_scale_factor, coord_scale_factor
    real(r8) :: plotfile_dt, partition_ds
    logical :: write_plotfile

    namelist /toolpath/ name, start_time, start_coord, time_scale_factor, coord_scale_factor, &
                        command_string, command_file, write_plotfile, plotfile_dt, partition_ds

    is_IOP = (scl_rank()==1)  ! process rank 1 does the reading

    if (is_IOP) rewind(lun)
    n = 0

    do  ! until all TOOLPATH namelists have been read or an error occurs

      if (is_IOP) call seek_to_namelist(lun, 'TOOLPATH', found, iostat=ios)
      call scl_bcast(ios)
      if (ios /= 0) call re_halt('error reading input file: iostat=' // i_to_c(ios))

      call scl_bcast(found)
      if (.not.found) return  ! no further TOOLPATH namelists found

      n = n + 1
      call re_info('Reading TOOLPATH namelist #' // i_to_c(n))

      !! Read the namelist variables, assigning default values first.
      name = NULL_C
      start_time = NULL_R
      start_coord = NULL_R
      time_scale_factor = NULL_R
      coord_scale_factor = NULL_R
      command_string = NULL_C
      command_file = NULL_C
      write_plotfile = .false.
      plotfile_dt = NULL_R
      partition_ds = NULL_R

      if (is_IOP) read(lun,nml=toolpath,iostat=ios,iomsg=iom)
      call scl_bcast(ios)
      if (ios /= 0) call re_halt('error reading TOOLPATH namelist' // trim(iom))

      !! Broadcast the namelist variables.
      call scl_bcast(name)
      call scl_bcast(start_time)
      call scl_bcast(start_coord)
      call scl_bcast(time_scale_factor)
      call scl_bcast(coord_scale_factor)
      call scl_bcast(command_string)
      call scl_bcast(command_file)
      call scl_bcast(write_plotfile)
      call scl_bcast(plotfile_dt)
      call scl_bcast(partition_ds)

      !! Check the name.
      if (name == NULL_C .or. name == '') &
          call re_halt('NAME must be assigned a nonempty value')
      if (tp_fac%known_toolpath(name)) &
          call re_halt('already read a TOOLPATH namelist with this name: ' // trim(name))

      !! Check that we got either COMMAND_STRING or COMMAND_FILE, but not both.
      if (command_string == NULL_C .and. command_file == NULL_C) &
          call re_halt('either COMMAND_STRING or COMMAND_FILE must be specified')
      if (command_string /= NULL_C .and. command_file /= NULL_C) &
          call re_halt('cannot specify both COMMAND_STRING and COMMAND_FILE')

      !! Adjust the file path and check it is accessible.
      if (command_file /= NULL_C) then
        if (is_IOP) inquire(file=trim(command_file),exist=found)
        call scl_bcast(found)
        if (.not.found) call re_halt('COMMAND_FILE not found: ' // trim(command_file))
      end if

      !! Check scaling factors if specified
      if (time_scale_factor /= NULL_R) then
        if (time_scale_factor <= 0) call re_halt('TIME_SCALE_FACTOR must be > 0')
      end if
      if (coord_scale_factor /= NULL_R) then
        if (coord_scale_factor <= 0) call re_halt('COORD_SCALE_FACTOR must be > 0')
      end if

      !! Check START_COORD is complete if specified
      if (any(start_coord /= NULL_R)) then
        if (any(start_coord == NULL_R)) call re_halt('START_COORD not completely specified')
      end if

      !! Check PLOTFILE_DT if needed
      if (write_plotfile) then
        if (plotfile_dt == NULL_R) call re_halt('PLOTFILE_DT not specified')
        if (plotfile_dt <= 0.0_r8) call re_halt('PLOTFILE_DT must be > 0')
        if (time_scale_factor /= NULL_R) plotfile_dt = time_scale_factor * plotfile_dt
      end if

      if (partition_ds /= NULL_R) then
        if (partition_ds <= 0.0_r8) call re_halt('PARTITION_DS must be > 0')
        if (coord_scale_factor /= NULL_R) partition_ds = coord_scale_factor * partition_ds
      end if

      call make_toolpath_plist

    end do

  contains

    subroutine make_toolpath_plist

      use parameter_list_type

      type(parameter_list), pointer :: plist
      character(:), allocatable :: plotfile

      plist => tp_fac%toolpath_plist(trim(name))

      !! Form the input parameter list for the toolpath factory.
      if (start_time /= NULL_R) call plist%set('start-time', start_time)
      if (any(start_coord /= NULL_R)) call plist%set('start-coord', start_coord)
      if (time_scale_factor /= NULL_R) call plist%set('time-scale-factor', time_scale_factor)
      if (coord_scale_factor /= NULL_R) call plist%set('coord-scale-factor', coord_scale_factor)
      if (command_string /= NULL_C) then
        call plist%set('command-string', trim(command_string))
      else
        call plist%set('command-file', trim(command_file))
      end if
      if (is_IOP .and. write_plotfile) then
        plotfile = 'toolpath-' // trim(name) // '.dat'
        call plist%set('plotfile', plotfile)
        call plist%set('plotfile-dt', plotfile_dt)
      end if
      if (partition_ds /= NULL_R) call plist%set('partition-ds', partition_ds)

    end subroutine

  end subroutine read_toolpath_namelists

end module re_toolpath
