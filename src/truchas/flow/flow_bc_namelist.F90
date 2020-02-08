!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module flow_bc_namelist

  implicit none
  private

  public :: read_flow_bc_namelists

contains

  subroutine read_flow_bc_namelists(lun, params)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use parameter_list_type
    use parallel_communication, only: is_IOP, broadcast
    use input_utilities, only: seek_to_namelist, NULL_C, NULL_R, NULL_I
    use string_utilities, only: i_to_c, lower_case
    use material_model_driver, only: matl_model
    use truchas_logging_services

    integer, intent(in) :: lun
    type(parameter_list), intent(inout) :: params

    integer :: n, ios, data_size, matid
    logical :: found
    character(128) :: iom
    character(:), allocatable :: label
    type(parameter_list), pointer :: plist

    !! Namelist variables
    integer :: face_set_ids(32)
    real(r8) :: pressure, velocity(3), dsigma, inflow_temperature
    character(31) :: name, type, pressure_func, velocity_func, inflow_material
    namelist /flow_bc/ name, face_set_ids, type, pressure, pressure_func, velocity, &
        velocity_func, dsigma, inflow_material, inflow_temperature

    call TLS_info('')
    call TLS_info('Reading FLOW_BC namelists ...')

    if (is_IOP) rewind(lun)

    n = 0 ! namelist counter
    do ! until all FLOW_BC namelists have been read or an error occurs
      if (is_IOP) call seek_to_namelist(lun, 'flow_bc', found, iostat=ios)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
      call broadcast(found)
      if (.not.found) exit

      n = n + 1
      label = 'FLOW_BC[' // i_to_c(n) // ']'

      name = NULL_C
      face_set_ids = NULL_I
      type = NULL_C
      pressure = NULL_R
      velocity = NULL_R
      dsigma = NULL_R
      pressure_func = NULL_C
      velocity_func = NULL_C
      inflow_material = NULL_C
      inflow_temperature = NULL_R

      if (is_IOP) read(lun,nml=flow_bc,iostat=ios,iomsg=iom)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading ' // label // ' namelist: ' // trim(iom))

      call broadcast(name)
      call broadcast(face_set_ids)
      call broadcast(type)
      call broadcast(pressure)
      call broadcast(velocity)
      call broadcast(dsigma)
      call broadcast(pressure_func)
      call broadcast(velocity_func)
      call broadcast(inflow_material)
      call broadcast(inflow_temperature)

      !! A unique NAME is required; becomes the BC sublist parameter name.
      if (name == NULL_C) then
        call TLS_fatal(label // ': NAME not specified')
      else if (params%is_sublist(name)) then
        call TLS_fatal(label // ': another FLOW_BC namelist has this NAME: ' // trim(name))
      else
        plist => params%sublist(trim(name))
      end if

      !! FACE_SET_IDS is required.
      if (count(face_set_ids /= NULL_I) == 0) then
        call TLS_fatal(label // ': FACE_SET_IDS not specified')
      else
        call plist%set('face-set-ids', pack(face_set_ids, mask=(face_set_ids/=NULL_I)))
      end if

      !! Check the TYPE value (required)
      select case (lower_case(type))
      case (NULL_C)
        call TLS_fatal(label // ': TYPE not specified')
      case ('velocity')
        if (any(velocity /= NULL_R) .and. velocity_func /= NULL_C) then
          call TLS_fatal(label // ': cannot specify both VELOCITY and VELOCITY_FUNC')
        else if (any(velocity /= NULL_R)) then
          if (any(velocity == NULL_R)) call TLS_fatal(label // ': 3-vector required for VELOCITY')
          call plist%set('velocity', velocity)
        else if (velocity_func /= NULL_C) then
          !TODO: verify a vector function with this name exists
          call plist%set('velocity', trim(velocity_func))
        else
          call TLS_fatal(label // ': either VELOCITY or VELOCITY_FUNC is required')
        end if
      case ('pressure')
        if (pressure /= NULL_R .and. pressure_func /= NULL_C) then
          call TLS_fatal(label // ': cannot specify both PRESSURE and PRESSURE_FUNC')
        else if (pressure /= NULL_R) then
          call plist%set('pressure', pressure)
        else if (pressure_func /= NULL_C) then
          !TODO: verify a scalar function with this name exists
          call plist%set('pressure', trim(pressure_func))
        else
          call TLS_fatal(label // ': either PRESSURE or PRESSURE_FUNC is required')
        end if
      case ('free-slip')
      case ('no-slip')
      case ('marangoni')
        if (dsigma /= NULL_R) then
          call plist%set('dsigma', dsigma)
        else
          call TLS_fatal(label // ': DSIGMA is required')
        end if
      case default
        call TLS_fatal(label // ': unknown TYPE: "' // trim(type) // '"')
      end select
      call plist%set('type', trim(type))

      !! Additional inflow data
      select case (lower_case(type))
      case ('velocity', 'pressure')
        if (inflow_material /= NULL_C) then
          matid = matl_model%phase_index(inflow_material)
          if (matid == 0) then
            call TLS_fatal('unknown material for INFLOW_MATERIAL: ' // trim(inflow_material))
          else if (.not.matl_model%is_fluid(matid)) then
            call TLS_fatal('INFLOW_MATERIAL not a liquid: ' // trim(inflow_material))
          end if
          call plist%set('inflow-material', trim(inflow_material))
        end if
        if (inflow_temperature /= NULL_R) call plist%set('inflow-temperature', inflow_temperature)
      end select

    end do

  end subroutine read_flow_bc_namelists

end module flow_bc_namelist
