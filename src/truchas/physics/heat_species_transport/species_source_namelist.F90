!!
!! SPECIES_SOURCE_NAMELIST
!!
!! Provides a subroutine for reading the SPECIES_SOURCE namelists for a file
!! and storing the data in a parameter list for later use.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module species_source_namelist

  implicit none
  private

  public :: read_species_source_namelists

contains

  subroutine read_species_source_namelists(lun, params)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use parallel_communication, only: is_IOP, broadcast
    use input_utilities, only: seek_to_namelist, NULL_C, NULL_R, NULL_I
    use string_utilities, only: i_to_c
    use physics_module, only: number_of_species
    use parameter_list_type
    use truchas_logging_services

    integer, intent(in) :: lun
    type(parameter_list), intent(inout) :: params

    integer :: n, ios
    logical :: found
    character(:), allocatable :: label, src_par
    character(128) :: iom
    type(parameter_list), pointer :: plist

    !! Namelist variables
    integer :: comp_id, cell_set_ids(100)
    real(r8) :: source
    character(32) :: source_func
    namelist /species_source/ comp_id, cell_set_ids, source, source_func

    call TLS_info('Reading SPECIES_SOURCE namelists ...')

    if (is_IOP) rewind(lun)

    n = 0 ! namelist counter
    do ! until all SPECIES_SOURCE namelists have been read or an error occurs

      if (is_IOP) call seek_to_namelist(lun, 'species_source', found, iostat=ios)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
      call broadcast(found)
      if (.not.found) exit

      n = n + 1
      label = 'SPECIES_SOURCE[' // i_to_c(n) // ']'

      comp_id = 1
      cell_set_ids = NULL_I
      source = NULL_R
      source_func = NULL_C

      if (is_IOP) read(lun,nml=species_source,iostat=ios,iomsg=iom)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading ' // label // ' namelist: ' // trim(iom))

      call broadcast(comp_id)
      call broadcast(cell_set_ids)
      call broadcast(source)
      call broadcast(source_func)

      plist => params%sublist(label)

      if (comp_id == NULL_I) then
        call TLS_fatal(label // ': COMP_ID not specified')
      else if (comp_id < 1 .or. comp_id > number_of_species) then
        call TLS_fatal(label // ': invalid COMP_ID value')
      else
        src_par = 'source' // i_to_c(comp_id)
      end if

      if (any(cell_set_ids /= NULL_I)) then
        call plist%set('cell-set-ids', pack(cell_set_ids, mask=(cell_set_ids/=NULL_I)))
        if (source /= NULL_R .and. source_func /= NULL_C) then
          call TLS_fatal(label // ': both SOURCE and SOURCE_FUNC specified')
        else if (source /= NULL_R) then
          call plist%set(src_par, source)
        else if (source_func /= NULL_C) then
          call plist%set(src_par, trim(source_func))
        else
          call TLS_fatal(label // ': neither SOURCE nor SOURCE_FUNC specified')
        end if
      else
        call TLS_fatal(label // ': CELL_SET_IDS not specified')
      end if

    end do

    select case (n)
    case (0)
      call TLS_info('  none found')
    case (1)
      call TLS_info('  read 1 SPECIES_SOURCE namelist')
    case default
      call TLS_info('  read ' // i_to_c(n) // ' SPECIES_SOURCE namelists')
    end select

  end subroutine read_species_source_namelists

end module species_source_namelist
