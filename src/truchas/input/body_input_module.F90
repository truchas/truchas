!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This module is used to set matnum, nbody, body_temp, body_phi, body_vel, and
! body_mass from interfaces_module, according to user input.

MODULE BODY_INPUT_MODULE

  use kinds, only: r8
  use truchas_logging_services
  use interfaces_module
  use material_model_driver, only: matl_model
  implicit none
  private

  public :: INTERFACES_INPUT

CONTAINS

  SUBROUTINE INTERFACES_INPUT (lun)

    use parallel_info_module, only: p_info
    use parameter_module,     only: mbody

    integer, intent(in) :: lun

    logical :: body_namelist, fatal, found
    integer :: ios, ib
    character(128) :: message

    ! Initialize the number of bodies
    nbody = 0

    ! defaults
    Body_mass = 0
    Body_Vel = 0    

    if (p_info%IOP) rewind lun

    ! Read BODY geometry
    do ib = 1, mbody
      call BODY_INPUT (lun, found)
      if (.not. found) exit
    end do

    ! Make sure one body namelist has been read
    if (nbody == 0) call TLS_fatal ('no BODY namelists found')

  END SUBROUTINE INTERFACES_INPUT


  ! Read BODY namelist, put read data into place. If a namelist is found, then
  ! it is read and data is put into body data arrays. On succesful read, then
  ! nbody -> nbody + 1.
  SUBROUTINE BODY_INPUT (lun, body_namelist)

    use input_utilities,        only: seek_to_namelist, NULL_I, NULL_R, NULL_C
    use parallel_info_module,   only: p_info
    use parameter_module,       only: mbody, mphi
    use legacy_mesh_api,        only: ndim
    use pgslib_module,          only: PGSLIB_BCAST
    use scalar_func_factories,  only: alloc_const_scalar_func
    use scalar_func_table,      only: lookup_func

    ! Argument List
    integer, intent(in) :: lun
    logical :: body_namelist

    ! Local Variables
    logical :: fatal
    integer :: ioerror, is, n, material_number
    real(r8) :: phi(mphi), temperature, velocity(ndim)
    character(31) :: temperature_function
    character(128) :: message, fatal_error_string

    ! Namelist variables unneeded here (used by body_namelist.F90 to initialize body geometry)
    character(64) :: surface_name, axis, fill, material_name
    real(r8) :: height, length(3), radius, rotation_angle(3), rotation_pt(3), translation_pt(3)
    integer :: mesh_material_number(16)

    ! Define BODY Namelist
    namelist /body/ surface_name, axis, height, radius, length, fill, &
        rotation_angle, rotation_pt, translation_pt, &
        material_name, phi, temperature, temperature_function, velocity, mesh_material_number

    ! Initialize for error checking
    fatal = .false.
    fatal_error_string = 'BODY namelist input error!'

    ! Initialize the surface input
    material_name        = NULL_C
    temperature          = NULL_R
    temperature_function = NULL_C
    phi                  = 0
    Velocity             = 0

    ! Set error detection stuff
    IO_PE_ONLY: if (p_info%IOP) then
       ! Find namelist
       call seek_to_namelist (lun, 'BODY', found=body_namelist)

       ! Read namelist if found one
       if (body_namelist) then
          read (lun, NML = body, IOSTAT = ioerror)
          if (ioerror /= 0) then ! If read error, then didn't read namelist
             body_namelist = .false.
             fatal         = .true.
          end if
       end if
    end if IO_PE_ONLY

    ! Broadcast data just read in to all PE's.
    if (.not.p_info%UseGlobalServices) then
       call PGSLIB_BCAST (body_namelist)
       call PGSLIB_BCAST (phi)
       call PGSLIB_BCAST (temperature)
       call PGSLIB_BCAST (temperature_function)
       call PGSLIB_BCAST (Velocity)
       call PGSLIB_BCAST (material_name)
    end if

    ! Continue only if we found the namelist and read it without error.
    BODY_NML: if (body_namelist) then
       nbody = nbody + 1

       ! Read notice
       write (message, 15) nbody
15     format (' Reading BODY Namelist #',i2,' ...')
       call TLS_info ('')
       call TLS_info (message)

       ! Too many bodies is a fatal error
       if (nbody > mbody) then
          fatal = .true.
          write (fatal_error_string, 20) mbody
20        format('exceeded maximum number: mbody = ',i5)
       end if

       ! Check that either the temperature or a function name was read, and not both.
       if (temperature == NULL_R .eqv. temperature_function == NULL_C) then
          fatal_error_string = 'either TEMPERATURE or TEMPERATURE_FUNCTION must be specified'
          fatal = .true.
       end if

       ! Save body data if we still don't have any error
       if (.not.fatal) then
         if (material_name == NULL_C) then
            fatal_error_string = 'MATERIAL_NAME not specified'
            fatal = .true.
         else
            Matnum(nbody) = matl_model%phase_index(material_name) ! yes, really a phase name
            if (Matnum(nbody) <= 0) then
              fatal = .true.
              if (matl_model%has_matl(material_name)) then
                fatal_error_string = 'must specify the phase of MATERIAL_NAME: ' // trim(material_name)
              else
                fatal_error_string = 'unknown MATERIAL_NAME: ' // trim(material_name)
              end if
            end if
         end if
         Body_Phi(nbody,:) = phi
         Body_Vel(:,nbody) = Velocity
         
         ! Generate and store the body temperature function.
         if (temperature /= NULL_R) then
           call alloc_const_scalar_func (body_temp(nbody)%f, temperature)
         else
           call lookup_func (temperature_function, body_temp(nbody)%f)
           if (.not.allocated(body_temp(nbody)%f)) then
             fatal_error_string =  'unknown function name: ' // trim(temperature_function)
             fatal = .true.
           end if
         end if
       end if
    end if BODY_NML

    call TLS_fatal_if_any (fatal, fatal_error_string)

  END SUBROUTINE BODY_INPUT

END MODULE BODY_INPUT_MODULE
