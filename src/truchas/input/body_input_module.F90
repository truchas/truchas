!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This module is used to set matnum, nbody, body_temp, body_phi, body_vel, and
! body_mass from interfaces_module, according to user input.

MODULE BODY_INPUT_MODULE

  use kinds, only: r8
  use interfaces_module
  use material_model_driver, only: matl_model
  use truchas_logging_services
  implicit none
  private

  public :: INTERFACES_INPUT

CONTAINS

  SUBROUTINE INTERFACES_INPUT (lun)

    integer, intent(in) :: lun

    ! defaults
    Body_mass = 0
    Body_Vel = 0

    call read_body_namelists(lun)

  END SUBROUTINE INTERFACES_INPUT


  subroutine read_body_namelists(lun)

    use parallel_communication, only: is_IOP, broadcast
    use input_utilities, only: seek_to_namelist, NULL_C, NULL_I, NULL_R
    use parameter_module, only: mbody, mphi
    use scalar_func_factories, only: alloc_const_scalar_func
    use scalar_func_table, only: lookup_func
    use string_utilities, only: i_to_c

    integer, intent(in) :: lun

    !! Namelist variables used here
    integer :: material_number
    real(r8) :: phi(mphi), temperature, velocity(3)
    character(31) :: temperature_function
    !! Namelist variables unneeded here (used by body_namelist.F90 to initialize body geometry)
    character(64) :: surface_name, axis, fill, material_name
    real(r8) :: height, length(3), radius, rotation_angle(3), rotation_pt(3), translation_pt(3)
    integer :: mesh_material_number(16)
    namelist /body/ surface_name, axis, height, radius, length, fill, &
        rotation_angle, rotation_pt, translation_pt, &
        material_name, phi, temperature, temperature_function, velocity, mesh_material_number

    integer :: ios
    logical :: found
    character(:), allocatable :: label
    character(128) :: iom

    call TLS_info('Reading BODY namelists (first pass) ...')

    if (is_IOP) rewind(lun)

    nbody = 0 ! namelist counter
    do ! until all BODY namelists have been read or an error occurs

      if (is_IOP) call seek_to_namelist(lun, 'body', found, iostat=ios)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
      call broadcast(found)
      if (.not.found) exit

      nbody = nbody + 1
      if (nbody > mbody) call TLS_fatal('too many BODY namelists')
      label = 'BODY[' // i_to_c(nbody) // ']'

      material_name = NULL_C
      temperature = NULL_R
      temperature_function = NULL_C
      phi = 0
      velocity = 0

      if (is_IOP) read(lun,nml=body,iostat=ios,iomsg=iom)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading ' // label // ' namelist: ' // trim(iom))

      call broadcast(material_name)
      call broadcast(temperature)
      call broadcast(temperature_function)
      call broadcast(phi)
      call broadcast(velocity)

      if (temperature == NULL_R .eqv. temperature_function == NULL_C) &
          call TLS_fatal(label // ': either TEMPERATURE or TEMPERATURE_FUNCTION must be specified')

      if (material_name == NULL_C) then
        call TLS_fatal(label // ': MATERIAL_NAME not specified')
      else
        matnum(nbody) = matl_model%phase_index(material_name) ! yes, really a phase name
        if (matnum(nbody) <= 0) then
          if (matl_model%has_matl(material_name)) then
            call TLS_fatal(label // ': must specify the phase of MATERIAL_NAME: ' // trim(material_name))
          else
            call TLS_fatal(label // ': unknown MATERIAL_NAME: ' // trim(material_name))
          end if
        end if
      end if

      body_phi(nbody,:) = phi
      body_vel(:,nbody) = velocity

      ! Generate and store the body temperature function.
      if (temperature /= NULL_R) then
        call alloc_const_scalar_func(body_temp(nbody)%f, temperature)
      else
        call lookup_func(temperature_function, body_temp(nbody)%f)
        if (.not.allocated(body_temp(nbody)%f)) &
            call TLS_fatal(label // ': unknown function name: ' // trim(temperature_function))
      end if
    end do

    select case (nbody)
    case (0)
      call TLS_fatal('no BODY namelist found')
    case (1)
      call TLS_info('  read 1 BODY namelist')
    case default
      call TLS_info('  read ' // i_to_c(nbody) // ' BODY namelists')
    end select

  end subroutine read_body_namelists

end module body_input_module
