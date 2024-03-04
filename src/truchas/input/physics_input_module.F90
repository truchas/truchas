!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module physics_input_module

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use truchas_logging_services
  implicit none
  private

  public :: physics_input

contains

  subroutine physics_check

    use physics_module

    logical :: fatal

    fatal = .false.

    ! Check the number_of_species value.
    if (species_transport) then
      if (number_of_species <= 0) then
        call TLS_error ('NUMBER_OF_SPECIES must be > 0')
        fatal = .true.
      end if
    end if

    if (fatal) call TLS_fatal('Errors found with PHYSICS namelists variables')

  end subroutine physics_check

  subroutine physics_input(lun)

    use physics_module
    use ih_driver,              only: SET_EM_SIMULATION_ON_OR_OFF
    use input_utilities,        only: seek_to_namelist, NULL_C
    use diffusion_solver_data,  only: ds_enabled, system_type, num_species
    use material_model_driver,  only: materials, nmat
    use parallel_communication, only: is_IOP, broadcast
    use string_utilities, only: i_to_c

    integer, intent(in) :: lun

    logical :: found
    integer :: ios
    character(80) :: iom

    namelist /physics/ heat_transport, species_transport, number_of_species, flow, &
                       body_force_density, prescribed_flow, &
                       electromagnetics, solid_mechanics, materials

    call TLS_info ('Reading PHYSICS namelist ...')

    !! Locate the PHYSICS namelist (required)
    if (is_IOP) then
      rewind(lun)
      call seek_to_namelist(lun, 'PHYSICS', found, iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
    call broadcast(found)
    if (.not.found) call TLS_fatal('PHYSICS namelist not found')

    !! Read the PHYSICS namelist, assigning default values first
    if (is_IOP) then
      flow = .false.
      heat_transport = .false.
      species_transport = .false.
      number_of_species = 0
      electromagnetics = .false.
      solid_mechanics = .false.
      body_force_density = 0.0_r8
      prescribed_flow = .false.
      materials = NULL_C
      read(lun,nml=physics,iostat=ios,iomsg=iom)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading PHYSICS namelist: ' // trim(iom))

    !! Broadcast all the namelist variables
    call broadcast(flow)
    call broadcast(body_force_density)
    call broadcast(heat_transport)
    call broadcast(species_transport)
    call broadcast(number_of_species)
    call broadcast(solid_mechanics)
    call broadcast(electromagnetics)
    call broadcast(prescribed_flow)
    call broadcast(materials)

    ! flow and prescribed_flow are mutually exclusive

    ! Check for stupid input errors.
    call physics_check

    !NNC: temporary test code
    if (prescribed_flow .and. (heat_transport .or. species_transport &
                                        .or. solid_mechanics .or. electromagnetics)) then
       call TLS_fatal('prescribed_flow is not compatible with any other physics')
    end if

    ! for now advection with a user-specified flow uses the flow solver
    flow = flow .or. prescribed_flow

    ! Electromagnetics is a local variable.  Now need to
    ! store the result in EM_Data_Store.
    ! Note that this handles the broadcast, also.
    call SET_EM_SIMULATION_ON_OR_OFF(electromagnetics)

    !! Configure the heat/species transport kernel.
    ds_enabled = heat_transport .or. species_transport
    if (species_transport) num_species = number_of_species
    if (heat_transport .and. species_transport) then
      system_type = 'thermal+species'
    else if (heat_transport) then
      system_type = 'thermal'
    else if (species_transport) then
      system_type = 'species'
    end if

    nmat = count(materials /= NULL_C)
    if (nmat == 0) call TLS_fatal('MATERIALS not specified')
    materials(1:nmat) = pack(materials, mask=(materials /= NULL_C))

  end subroutine physics_input

end module physics_input_module
