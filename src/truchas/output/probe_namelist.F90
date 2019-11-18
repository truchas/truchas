!!
!! PROBE_NAMELIST
!!
!! This module provides a procedure for reading the PROBE namelists.
!! The data is copied into a parameter list module variable for later use.
!!
!! Michael Hall <hall@lanl.gov>
!! July 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module probe_namelist

  use parameter_list_type
  implicit none
  private

  public :: read_probe_namelists

  type(parameter_list), public :: params

contains

  subroutine read_probe_namelists(lun)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use input_utilities, only: seek_to_namelist, NULL_C, NULL_R, NULL_I
    use parallel_communication, only: is_IOP, broadcast
    use string_utilities, only: i_to_c
    use truchas_logging_services

    integer, intent(in) :: lun

    integer :: n, ios
    logical :: found
    character(80) :: iom
    character(:), allocatable :: label
    type(parameter_list), pointer :: sublist

    !! Namelist variables
    character(31) :: data_file, data
    character(80) :: description
    real(r8) :: coord(3), coord_scale_factor
    integer  :: digits
    namelist /probe/ data_file, description, coord, coord_scale_factor, data, digits

    call TLS_info ('')
    call TLS_info ('Reading PROBE namelists ...')

    if (is_IOP) rewind(lun)

    n = 0 ! namelist counter
    do ! until all PROBE namelists have been read

      if (is_IOP) call seek_to_namelist(lun, 'PROBE', found, iostat=ios)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('Error reading input file: iostat=' // i_to_c(ios))
      call broadcast(found)
      if (.not.found) exit

      n = n + 1
      label = 'PROBE[' // i_to_c(n) // ']'

      data_file = NULL_C
      description = NULL_C
      coord = NULL_R
      coord_scale_factor = NULL_R
      data = NULL_C
      digits = NULL_I

      if (is_IOP) read(lun,nml=probe,iostat=ios,iomsg=iom)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading ' // label // ' namelist: ' // trim(iom))

      call broadcast(data_file)
      call broadcast(description)
      call broadcast(coord)
      call broadcast(coord_scale_factor)
      call broadcast(data)
      call broadcast(digits)

      !! Require a unique DATA_FILE name; also use for the sublist parameter name.
      if (data_file == NULL_C) then
        call TLS_fatal(label // ': DATA_FILE not specified')
      else if (params%is_sublist(data_file)) then
        call TLS_fatal(label // ': another PROBE uses this DATA_FILE: ' // trim(data_file))
      else
        sublist => params%sublist(trim(data_file))
        call sublist%set('data-file', trim(data_file))
      end if

      !! Optional DESCRIPTION string.
      if (description /= NULL_C) call sublist%set('description', trim(description))

      !! Probe DATA type is required. For now accept anything and
      !! defer checking for valid values to probe initialization.
      select case (data)
      case (NULL_C)
        call TLS_fatal(label // ': DATA not specified')
      case default
        call sublist%set('data', trim(data))
      end select

      !! Check required COORD.
      if (all(coord == NULL_R)) then
        call TLS_fatal(label // ': COORD not specified')
      else if (any(coord == NULL_R)) then
        call TLS_fatal(label // ': COORD requires 3 values')
      else
        call sublist%set('coord', coord)
      end if

      !! Check optional COORD_SCALE_FACTOR.
      if (coord_scale_factor /= NULL_R) then
        if (coord_scale_factor <= 0.0_r8) then
          call TLS_fatal(label // ': COORD_SCALE_FACTOR must be > 0.0')
        else
          call sublist%set('coord-scale-factor', coord_scale_factor)
        end if
      end if

      !! Check optional DIGITS.
      if (digits /= NULL_I) then
        if (digits < 2) then
          call TLS_fatal(label // ': DIGITS must be > 1')
        else
          call sublist%set('digits', digits)
        end if
      end if

    end do

    if (n == 0) then
      call TLS_info('  no PROBE namelists found')
    else
      call TLS_info('  read ' // i_to_c(n) // ' PROBE namelists')
    end if

  end subroutine read_probe_namelists

end module probe_namelist
