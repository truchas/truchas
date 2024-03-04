!!
!! ELECTROMAGNETIC_BC_NAMELIST
!!
!! Zach Jibben <zjibben@lanl.gov>
!! August 2023
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module electromagnetic_bc_namelist

  implicit none
  private

  public :: read_electromagnetic_bc_namelists

contains

  subroutine read_electromagnetic_bc_namelists(lun, params)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use parameter_list_type
    use parallel_communication, only: is_IOP, broadcast
    use string_utilities, only: i_to_c, lower_case
    use input_utilities, only: seek_to_namelist, NULL_C, NULL_I, NULL_R
    use truchas_logging_services

    integer, intent(in) :: lun
    type(parameter_list), intent(inout) :: params

    integer :: n, ios
    logical :: found
    character(128) :: iom
    character(:), allocatable :: label
    type(parameter_list), pointer :: plist

    character(128) :: name, type
    integer :: face_set_ids(128)

    namelist /electromagnetic_bc/ name, type, face_set_ids

    call TLS_info('Reading ELECTROMAGNETIC_BC namelists ...')

    if (is_IOP) rewind(lun)

    n = 0 ! namelist counter
    do ! until all ELECTROMAGNETIC_BC namelists have been read or an error occurs
      if (is_IOP) call seek_to_namelist(lun, 'electromagnetic_bc', found, iostat=ios)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
      call broadcast(found)
      if (.not.found) exit

      n = n + 1
      label = 'ELECTROMAGNETIC_BC[' // i_to_c(n) // ']'

      !! Read the namelist
      name = NULL_C
      type  = NULL_C
      face_set_ids = NULL_I

      if (is_IOP) read(lun, nml=electromagnetic_bc, iostat=ios, iomsg=iom)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading ' // label // ' namelist: ' // trim(iom))

      call broadcast(name)
      call broadcast(type)
      call broadcast(face_set_ids)

      !! A unique NAME is required; becomes the BC sublist parameter name.
      if (name == NULL_C) then
        call TLS_fatal(label // ': NAME not specified')
      else if (params%is_sublist(name)) then
        call TLS_fatal(label // ': another ELECTROMAGNETIC_BC namelist has this NAME: ' // trim(name))
      else
        plist => params%sublist(trim(name))
      end if

      !! FACE_SET_IDS is required; cannot check that they are valid at this point.
      if (count(face_set_ids /= NULL_I) == 0) then
        call TLS_fatal(label // ': FACE_SET_IDS not specified')
      else
        call plist%set('face-set-ids', pack(face_set_ids, mask=(face_set_ids/=NULL_I)))
      end if

      !! Check the required TYPE value.
      if (type == NULL_C) call TLS_fatal(label // ': TYPE not specified')
      select case (lower_case(type))
      case ('pec')
        ! no parameters
      case ('ih-hfield')
        ! no parameters
      case default
        call TLS_fatal(label // ': unknown TYPE: ' // trim(type))
      end select
      call plist%set('type', trim(type))

      call TLS_info('  read namelist "' // trim(name) // '"')
    end do

    if (n == 0) call TLS_info('  none found')

  end subroutine read_electromagnetic_bc_namelists

end module electromagnetic_bc_namelist
