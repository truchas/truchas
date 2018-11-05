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
    use truchas_logging_services

    integer, intent(in) :: lun
    type(parameter_list), intent(inout) :: params

    integer :: n, ios, data_size
    logical :: found
    character(128) :: iom
    character(:), allocatable :: label
    type(parameter_list), pointer :: plist

    !! Namelist variables
    integer :: face_set_ids(32)
    real(r8) :: data_constant(3)
    character(31) :: name, condition, data_function
    namelist /flow_bc/ name, face_set_ids, condition, data_constant, data_function

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
      condition = NULL_C
      data_constant = NULL_R
      data_function = NULL_C

      if (is_IOP) read(lun,nml=flow_bc,iostat=ios,iomsg=iom)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading ' // label // ' namelist: ' // trim(iom))

      call broadcast(name)
      call broadcast(face_set_ids)
      call broadcast(condition)
      call broadcast(data_constant)
      call broadcast(data_function)

      !! A unique NAME is required; becomes the BC sublist parameter name.
      if (name == NULL_C) then
        call TLS_fatal(label // ': NAME not specified')
      else if (params%is_sublist(name)) then
        call TLS_fatal(label // ': another FLOW_BC namelist has this NAME: ' // trim(name))
      else
        plist => params%sublist(name)
      end if

      !! FACE_SET_IDS is required.
      if (count(face_set_ids /= NULL_I) == 0) then
        call TLS_fatal(label // ': FACE_SET_IDS not specified')
      else
        call plist%set('face sets', pack(face_set_ids, mask=(face_set_ids/=NULL_I)))
      end if

      !! Check the CONDITION value (required)
      select case (lower_case(condition))
      case (NULL_C)
        call TLS_fatal(label // ': CONDITION not specified')
      case ('velocity dirichlet')
        data_size = 2
      case ('pressure dirichlet')
        data_size = 1
      case ('pressure neumann')
        data_size = 0
      case ('slip')
        data_size = 0
      case ('no slip')
        data_size = 0
      case ('surface tension')
        data_size = 1
      case default
        call TLS_fatal(label // ': unknown CONDITION: "' // trim(condition) // '"')
      end select
      call plist%set('condition', trim(condition))

      !! Check the DATA_CONSTANT or DATA_FUNCTION values
      if (any(data_constant /= NULL_R) .and. data_function /= NULL_C) then
        call TLS_fatal(label // ': cannot specify both DATA_CONSTANT and DATA_FUNCTION')
      else if (any(data_constant /= NULL_R)) then
        select case (data_size)
        case (0) ! no data needed
          call TLS_warn(label // ': specified DATA_CONSTANT value will be ignored')
        case (1) ! scalar data needed
          if (any(data_constant(2:) /= NULL_R) .or. data_constant(1) == NULL_R) &
              call TLS_fatal(label // ': require scalar value for DATA_CONSTANT')
          call plist%set('data', data_constant(1))
        case (2) ! 3-vector data needed
          if (any(data_constant == NULL_R)) call TLS_fatal('require 3-vector for DATA_CONSTANT')
          call plist%set('data', data_constant)
        end select
      else if (data_function /= NULL_C) then
        !call TLS_fatal(label // ': DATA_FUNCTION not yet implemented')
        select case (data_size)
        case (0) ! no data needed
          call TLS_warn(label // ': specified DATA_FUNCTION value will be ignored')
        case (1) ! scalar data needed
          !TODO: verify a scalar function with this name exists
          call plist%set('data', trim(data_function))
        case (2) ! 3-vector data needed
          !TODO: verify a vector function with this name exists
          call plist%set('data', trim(data_function))
        end select
      else if (data_size > 0) then
        call TLS_fatal(label // ': either DATA_CONSTANT or DATA_FUNCTION is required')
      end if
    end do

  end subroutine read_flow_bc_namelists

end module flow_bc_namelist
