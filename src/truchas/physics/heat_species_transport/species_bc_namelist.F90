!!
!! SPECIES_BC_NAMELIST
!!
!! Provides a subroutine for reading the SPECIES_BC namelists for a file and
!! storing the data in a parameter list for later use.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! February 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module species_bc_namelist

  implicit none
  private

  public :: read_species_bc_namelists

contains

  subroutine read_species_bc_namelists(lun, params)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use parallel_communication, only: is_IOP, broadcast
    use input_utilities, only: seek_to_namelist, NULL_C, NULL_I, NULL_R
    use string_utilities, only: i_to_c, lower_case
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
    integer :: face_set_ids(100), comp_id
    real(r8) :: conc, flux, mtc, ambient_conc
    character(32) :: name, type, conc_func, flux_func, mtc_func, ambient_conc_func
    namelist /species_bc/ name, face_set_ids, comp_id, type, conc, conc_func, flux, flux_func, &
        mtc, mtc_func, ambient_conc, ambient_conc_func


    call TLS_info('Reading SPECIES_BC namelists ...')

    if (is_IOP) rewind(lun)

    n = 0 ! namelist counter
    do ! until all SPECIES_BC namelists have been read or an error occurs

      if (is_IOP) call seek_to_namelist(lun, 'species_bc', found, iostat=ios)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
      call broadcast(found)
      if (.not.found) exit

      n = n + 1
      label = 'SPECIES_BC[' // i_to_c(n) // ']'

      name = NULL_C
      face_set_ids = NULL_I
      comp_id = NULL_I
      type = NULL_C
      conc = NULL_R
      conc_func = NULL_C
      flux = NULL_R
      flux_func = NULL_C
      mtc = NULL_R
      mtc_func = NULL_C
      ambient_conc = NULL_R
      ambient_conc_func = NULL_C

      if (is_IOP) read(lun,nml=species_bc,iostat=ios,iomsg=iom)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading ' // label // ' namelist: ' // trim(iom))

      call broadcast(name)
      call broadcast(face_set_ids)
      call broadcast(comp_id)
      call broadcast(type)
      call broadcast(conc)
      call broadcast(conc_func)
      call broadcast(flux)
      call broadcast(flux_func)
      call broadcast(mtc)
      call broadcast(mtc_func)
      call broadcast(ambient_conc)
      call broadcast(ambient_conc_func)

      !! A unique NAME is required; becomes the BC sublist parameter name.
      if (name == NULL_C) then
        call TLS_fatal(label // ': NAME not specified')
      else if (params%is_sublist(name)) then
        call TLS_fatal(label // ': another SPECIES_BC namelist has this NAME: ' // trim(name))
      else
        plist => params%sublist(trim(name))
      end if

      !! FACE_SET_IDS is required; cannot check that they are valid at this point.
      if (count(face_set_ids /= NULL_I) == 0) then
        call TLS_fatal(label // ': FACE_SET_IDS not specified')
      else
        call plist%set('face-set-ids', pack(face_set_ids, mask=(face_set_ids/=NULL_I)))
      end if

      !! COMP_ID is optional; will default to 1 if not set.
      if (comp_id /= NULL_I) then
        if (comp_id <= 0) then
          call TLS_fatal(label // ': COMP_ID must be > 0')
        else
          call plist%set('comp', comp_id)
        end if
      end if

      !! Check the required TYPE value.
      select case (lower_case(type))
      case ('concentration', 'flux', 'mtc', 'interface-mtc')
      case (NULL_C)
        call TLS_fatal(label // ': TYPE not specified')
      case default
        call TLS_fatal(label // ': unknown TYPE: ' // trim(type))
      end select
      call plist%set('type', trim(type))

      select case (lower_case(type))
      case ('concentration')

        if (conc /= NULL_R .and. conc_func /= NULL_C) then
          call TLS_fatal(label // ': both CONC and CONC_FUNC specified')
        else if (conc /= NULL_R) then
          call plist%set('conc', conc)
        else if (conc_func /= NULL_C) then
          call plist%set('conc', trim(conc_func))
        else
          call TLS_fatal(label // ': neither CONC or CONC_FUNC specified')
        end if

        !! NB: any other specified variables are silently ignored

      case ('flux')

        if (flux /= NULL_R .and. flux_func /= NULL_C) then
          call TLS_fatal(label // ': both FLUX and FLUX_FUNC specified')
        else if (flux /= NULL_R) then
          call plist%set('flux', flux)
        else if (flux_func /= NULL_C) then
          call plist%set('flux', trim(flux_func))
        else
          call TLS_fatal(label // ': neither FLUX or FLUX_FUNC specified')
        end if

        !! NB: any other specified variables are silently ignored

      case ('mtc')

        if (mtc /= NULL_R .and. mtc_func /= NULL_C) then
          call TLS_fatal(label // ': both MTC and MTC_FUNC specified')
        else if (mtc /= NULL_R) then
          call plist%set('mtc', mtc)
        else if (mtc_func /= NULL_C) then
          call plist%set('mtc', trim(mtc_func))
        else
          call TLS_fatal(label // ': neither MTC or MTC_FUNC specified')
        end if

        if (ambient_conc /= NULL_R .and. ambient_conc_func /= NULL_C) then
          call TLS_fatal(label // ': both AMBIENT_CONC and AMBIENT_CONC_FUNC specified')
        else if (ambient_conc /= NULL_R) then
          call plist%set('ambient-conc', ambient_conc)
        else if (ambient_conc_func /= NULL_C) then
          call plist%set('ambient-conc', trim(ambient_conc_func))
        else
          call TLS_fatal(label // ': neither AMBIENT_CONC or AMBIENT_CONC_FUNC specified')
        end if

        !! NB: any other specified variables are silently ignored

      case ('interface-mtc')

        if (mtc /= NULL_R .and. mtc_func /= NULL_C) then
          call TLS_fatal(label // ': both MTC and MTC_FUNC specified')
        else if (mtc /= NULL_R) then
          if (mtc < 0) call TLS_fatal(label // ': mtc < 0.0')
          call plist%set('mtc', mtc)
        else if (mtc_func /= NULL_C) then
          call plist%set('mtc', mtc_func)
        else
          call TLS_fatal(label // ': neither MTC or MTC_FUNC specified')
        end if

        !! NB: any other specified variables are silently ignored

      end select

    end do

    select case (n)
    case (0)
      call TLS_info('  none found')
    case (1)
      call TLS_info('  read 1 SPECIES_BC namelist')
    case default
      call TLS_info('  read ' // i_to_c(n) // ' SPECIES_BC namelists')
    end select

  end subroutine read_species_bc_namelists

end module species_bc_namelist
