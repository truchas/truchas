!!
!! THERMAL_SOURCE_NAMELIST
!!
!! Provides a subroutine for reading the THERMAL_SOURCE namelists for a file and
!! storing the data in a parameter list for later use.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! February 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module thermal_source_namelist

  implicit none
  private

  public :: read_thermal_source_namelists

contains

  subroutine read_thermal_source_namelists(lun, params)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use parallel_communication, only: is_IOP, broadcast
    use input_utilities, only: seek_to_namelist, NULL_C, NULL_R, NULL_I
    use string_utilities, only: i_to_c
    use parameter_list_type
    use truchas_logging_services

    integer, intent(in) :: lun
    type(parameter_list), intent(inout) :: params

    integer :: n, ios
    logical :: found
    character(:), allocatable :: label
    character(128) :: iom
    type(parameter_list), pointer :: plist

    !! Namelist variables
    integer :: cell_set_ids(100)
    real(r8) :: prefactor, source
    character(32) :: name, prefactor_func, source_func
    character(128) :: data_file
    namelist /thermal_source/ name, data_file, prefactor, prefactor_func, &
        cell_set_ids, source, source_func

    call TLS_info('Reading THERMAL_SOURCE namelists ...')

    if (is_IOP) rewind(lun)

    n = 0 ! namelist counter
    do ! until all THERMAL_SOURCE namelists have been read or an error occurs

      if (is_IOP) call seek_to_namelist(lun, 'thermal_source', found, iostat=ios)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
      call broadcast(found)
      if (.not.found) exit

      n = n + 1
      label = 'THERMAL_SOURCE[' // i_to_c(n) // ']'

      name = NULL_C
      data_file = NULL_C
      prefactor = NULL_R
      prefactor_func = NULL_C
      cell_set_ids = NULL_I
      source = NULL_R
      source_func = NULL_C

      if (is_IOP) read(lun,nml=thermal_source,iostat=ios,iomsg=iom)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading ' // label // ' namelist: ' // trim(iom))

      call broadcast(name)
      call broadcast(data_file)
      call broadcast(prefactor)
      call broadcast(prefactor_func)
      call broadcast(cell_set_ids)
      call broadcast(source)
      call broadcast(source_func)

      !! A unique NAME is required; becomes the source sublist parameter name.
      if (name == NULL_C) then
        call TLS_fatal(label // ': NAME not specified')
      else if (params%is_sublist(name)) then
        call TLS_fatal(label // ': another THERMAL_SOURCE namelist has this NAME: ' // trim(name))
      else
        plist => params%sublist(trim(name))
      end if

      if (data_file /= NULL_C) then
        call plist%set('data-file', data_file)
        if (prefactor /= NULL_R .and. prefactor_func /= NULL_C) then
          call TLS_fatal(label // ': both PREFACTOR and PREFACTOR_FUNC specified')
        else if (prefactor /= NULL_R) then
          call plist%set('prefactor', prefactor)
        else if (prefactor_func /= NULL_C) then
          call plist%set('prefactor', trim(prefactor_func))
        end if
      else if (any(cell_set_ids /= NULL_I)) then
        call plist%set('cell-set-ids', pack(cell_set_ids, mask=(cell_set_ids/=NULL_I)))
        if (source /= NULL_R .and. source_func /= NULL_C) then
          call TLS_fatal(label // ': both SOURCE and SOURCE_FUNC specified')
        else if (source /= NULL_R) then
          call plist%set('source', source)
        else if (source_func /= NULL_C) then
          call plist%set('source', trim(source_func))
        end if
      else
        call TLS_fatal(label // ': neither CELL_SET_IDS nor DATA_FILE specified')
      end if

    end do

    select case (n)
    case (0)
      call TLS_info('  none found')
    case (1)
      call TLS_info('  read 1 THERMAL_SOURCE namelist')
    case default
      call TLS_info('  read ' // i_to_c(n) // ' THERMAL_SOURCE namelists')
    end select

  end subroutine read_thermal_source_namelists

end module thermal_source_namelist
