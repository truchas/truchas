!!
!! PBF_MATERIAL_NAMELIST
!!
!! Neil N. Carlson <neil.n.carlson@gmail.com>
!! June 2023
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module pbf_material_namelist

  implicit none
  private

  public :: read_pbf_material_namelist

contains

  subroutine read_pbf_material_namelist(lun, params)

    use parameter_list_type
    use input_utilities, only: seek_to_namelist, NULL_C
    use string_utilities, only: i_to_c
    use parallel_communication, only: is_IOP, broadcast
    use truchas_logging_services

    integer, intent(in) :: lun
    type(parameter_list), pointer :: params

    integer :: ios
    logical :: found
    character(128) :: iom

    character(31) :: material1, material2
    namelist /pbf_material/ material1, material2

    params => null()

    if (is_IOP) rewind(lun)

    if (is_IOP) call seek_to_namelist(lun, 'pbf_material', found, iostat=ios)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
    call broadcast(found)
    if (.not.found) return ! namelist is optional

    call TLS_info('Reading PBF_MATERIAL namelist ...')

    allocate(params)

    material1 = NULL_C
    material2 = NULL_C

    if (is_IOP) read(lun,nml=pbf_material,iostat=ios,iomsg=iom)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading PBF_MATERIAL namelist: ' // trim(iom))

    call broadcast(material1)
    call broadcast(material2)

    if (material1 /= NULL_C) call params%set('material1', material1)
    if (material2 /= NULL_C) call params%set('material2', material2)
    
  end subroutine read_pbf_material_namelist

end module pbf_material_namelist
