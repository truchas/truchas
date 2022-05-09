!!
!! TOOLPATH_NAMELIST
!!
!! This module provides a procedure for reading the TOOLPATH namelists.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! November 2016
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module toolpath_namelist

  use toolpath_factory_type
  implicit none
  private

  public :: read_toolpath_namelists

contains

  subroutine read_toolpath_namelists(lun, tp_fac)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use parallel_communication, only: is_IOP, broadcast
    use string_utilities, only: i_to_c
    use input_utilities, only: seek_to_namelist, NULL_R, NULL_C
    use truchas_env, only: input_dir, output_dir
    use truchas_logging_services

    integer, intent(in) :: lun
    type(toolpath_factory), intent(inout) :: tp_fac

    logical :: found
    integer :: n, ios
    character(128) :: iom

    !! namelist variables
    character(31) :: name
    character(256) :: command_file
    character(1000) :: command_string
    character(:), allocatable :: label
    real(r8) :: start_time, start_coord(3), time_scale_factor, coord_scale_factor
    real(r8) :: plotfile_dt, partition_ds
    logical :: write_plotfile

    namelist /toolpath/ name, start_time, start_coord, time_scale_factor, coord_scale_factor, &
                        command_string, command_file, write_plotfile, plotfile_dt, partition_ds

    call TLS_info('Reading TOOLPATH namelists ...')

    if (is_IOP) rewind(lun)
    n = 0

    do  ! until all TOOLPATH namelists have been read or an error occurs

      if (is_IOP) call seek_to_namelist(lun, 'TOOLPATH', found, iostat=ios)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading input file')

      call broadcast(found)
      if (.not.found) exit  ! no further TOOLPATH namelists found

      n = n + 1
      label = 'TOOLPATH[' // i_to_c(n) // ']'

      !! Read the namelist variables, assigning default values first.
      if (is_IOP) then
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
        read(lun,nml=toolpath,iostat=ios,iomsg=iom)
      end if
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading ' // label // ' namelist: ' // trim(iom))

      !! Broadcast the namelist variables.
      call broadcast(name)
      call broadcast(start_time)
      call broadcast(start_coord)
      call broadcast(time_scale_factor)
      call broadcast(coord_scale_factor)
      call broadcast(command_string)
      call broadcast(command_file)
      call broadcast(write_plotfile)
      call broadcast(plotfile_dt)
      call broadcast(partition_ds)

      !! Check the name.
      if (name == NULL_C .or. name == '') &
          call TLS_fatal(label // ': NAME must be assigned a nonempty value')
      if (tp_fac%known_toolpath(trim(name))) &
          call TLS_fatal(label // ': another TOOLPATH namelist has this NAME: ' // trim(name))

      !! Check that we got either COMMAND_STRING or COMMAND_FILE, but not both.
      if (command_string == NULL_C .and. command_file == NULL_C) &
          call TLS_fatal(label // ': either COMMAND_STRING or COMMAND_FILE must be specified')
      if (command_string /= NULL_C .and. command_file /= NULL_C) &
          call TLS_fatal(label // ': cannot specify both COMMAND_STRING and COMMAND_FILE')

      !! Adjust the file path and check it is accessible.
      if (command_file /= NULL_C) then
        if (command_file(1:1) /= '/') then ! not an absolute path
          command_file = trim(input_dir) // trim(command_file) ! make relative to input dir
        end if
        if (is_IOP) inquire(file=trim(command_file),exist=found)
        call broadcast(found)
        if (.not.found) call TLS_fatal(label // ': COMMAND_FILE not found: ' // trim(command_file))
      end if

      !! Check scaling factors if specified
      if (time_scale_factor /= NULL_R) then
        if (time_scale_factor <= 0) call TLS_fatal(label // ': TIME_SCALE_FACTOR must be > 0')
      end if
      if (coord_scale_factor /= NULL_R) then
        if (coord_scale_factor <= 0) call TLS_fatal(label // ': COORD_SCALE_FACTOR must be > 0')
      end if

      !! Check START_COORD is complete if specified
      if (any(start_coord /= NULL_R)) then
        if (any(start_coord == NULL_R)) call TLS_fatal(label // ': START_COORD not completely specified')
      end if

      !! Check PLOTFILE_DT if needed
      if (write_plotfile) then
        if (plotfile_dt == NULL_R) call TLS_fatal(label // ': PLOTFILE_DT not specified')
        if (plotfile_dt <= 0.0_r8) call TLS_fatal(label // ': PLOTFILE_DT must be > 0')
        if (time_scale_factor /= NULL_R) plotfile_dt = time_scale_factor * plotfile_dt
      end if

      if (partition_ds /= NULL_R) then
        if (partition_ds <= 0.0_r8) call TLS_fatal(label // ': PARTITION_DS must be > 0')
        if (coord_scale_factor /= NULL_R) partition_ds = coord_scale_factor * partition_ds
      end if

      call TLS_info('  read namelist "' // trim(name) // '"')
      call make_toolpath_plist

    end do

    if (n == 0) call TLS_info('  none found')

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
      if (write_plotfile) then
        plotfile = trim(output_dir) // 'toolpath-' // trim(name) // '.dat'
        call plist%set('plotfile', plotfile)
        call plist%set('plotfile-dt', plotfile_dt)
      end if
      if (partition_ds /= NULL_R) call plist%set('partition-ds', partition_ds)

    end subroutine

  end subroutine read_toolpath_namelists

end module toolpath_namelist
