!!
!! Zach Jibben <zjibben@lanl.gov>
!! September 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module solid_mechanics_bc_namelist

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  public :: read_solid_mechanics_bc_namelists

contains

  subroutine read_solid_mechanics_bc_namelists(lun, params)

    use parameter_list_type
    use truchas_logging_services
    use parallel_communication, only: is_IOP, broadcast
    use string_utilities, only: i_to_c, lower_case
    use input_utilities, only: seek_to_namelist, NULL_C, NULL_R, NULL_I

    integer, intent(in) :: lun
    type(parameter_list), intent(inout) :: params

    integer :: n, ios
    logical :: found
    character(128) :: iom
    character(:), allocatable :: label
    type(parameter_list), pointer :: plist

    !! Namelist variables
    integer :: face_set_ids(32), node_set_ids(32)
    real(r8) :: displacement, traction
    character(32) :: name, type, displacement_func, traction_func
    namelist /solid_mechanics_bc/ name, face_set_ids, node_set_ids, type, &
        displacement, traction, displacement_func, traction_func

    call TLS_info('')
    call TLS_info('Reading SOLID_MECHANICS_BC namelists ...')

    if (is_IOP) rewind(lun)

    n = 0 ! namelist counter
    do ! until all SOLID_MECHANICS_BC namelists have been read or an error occurs
      if (is_IOP) call seek_to_namelist(lun, 'solid_mechanics_bc', found, iostat=ios)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
      call broadcast(found)
      if (.not.found) exit

      n = n + 1
      label = 'SOLID_MECHANICS_BC[' // i_to_c(n) // ']'

      name = NULL_C
      face_set_ids = NULL_I
      node_set_ids = NULL_I
      type = NULL_C
      displacement = NULL_R
      displacement_func = NULL_C
      traction = NULL_R
      traction_func = NULL_C

      if (is_IOP) read(lun,nml=solid_mechanics_bc,iostat=ios,iomsg=iom)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading ' // label // ' namelist: ' // trim(iom))

      call broadcast(name)
      call broadcast(face_set_ids)
      call broadcast(node_set_ids)
      call broadcast(type)
      call broadcast(displacement)
      call broadcast(displacement_func)
      call broadcast(traction)
      call broadcast(traction_func)

      !! A unique NAME is required; becomes the BC sublist parameter name.
      if (name == NULL_C) then
        call TLS_fatal(label // ': NAME not specified')
      else if (params%is_sublist(name)) then
        call TLS_fatal(label // ': another SOLID_MECHANICS_BC namelist has this NAME: ' // trim(name))
      else
        plist => params%sublist(trim(name))
      end if

      !! Check the TYPE value (required)
      select case (lower_case(type))
      case (NULL_C)
        call TLS_fatal(label // ': TYPE not specified')
      case ('displacement-x', 'displacement-y', 'displacement-z', 'displacement-n')
        if (displacement /= NULL_R .and. displacement_func /= NULL_C) then
          call TLS_fatal(label // ': cannot specify both DISPLACEMENT and DISPLACEMENT_FUNC')
        else if (displacement /= NULL_R) then
          call plist%set('displacement', displacement)
        else if (displacement_func /= NULL_C) then
          ! TODO: verify a function with this name exists
          call plist%set('displacement', trim(displacement_func))
        else
          call TLS_fatal(label // ': either DISPLACEMENT or DISPLACEMENT_FUNC is required')
        end if
      case ('traction-x', 'traction-y', 'traction-z', 'traction-n')
        if (traction /= NULL_R .and. traction_func /= NULL_C) then
          call TLS_fatal(label // ': cannot specify both TRACTION and TRACTION_FUNC')
        else if (traction /= NULL_R) then
          call plist%set('traction', traction)
        else if (traction_func /= NULL_C) then
          ! TODO: verify a function with this name exists
          call plist%set('traction', trim(traction_func))
        else
          call TLS_fatal(label // ': either TRACTION or TRACTION_FUNC is required')
        end if
      case ('gap-contact')
        ! no additional parameters
      case ('pinned-node')
        if (count(node_set_ids /= NULL_I) == 0) then
          call TLS_fatal(label // ': NODE_SET_IDS not specified')
        else
          call plist%set('node-set-ids', pack(node_set_ids, mask=(node_set_ids/=NULL_I)))
        end if
      case default
        call TLS_fatal(label // ': TYPE not recognized')
      end select
      call plist%set('type', trim(type))


      block
        logical :: displ_d, given_face_sets, given_node_sets
        displ_d = lower_case(type) == 'displacement-x' &
            .or.  lower_case(type) == 'displacement-y' &
            .or.  lower_case(type) == 'displacement-z'
        given_face_sets = count(face_set_ids /= NULL_I) /= 0
        given_node_sets = count(node_set_ids /= NULL_I) /= 0
        if (displ_d) then
          ! nodeset allowed, faceset or nodeset required
          if (.not.(given_face_sets .or. given_node_sets)) &
              call TLS_fatal(label // ': Neither FACE_SET_IDS nor NODE_SET_IDS specified')

          if (given_node_sets) &
              call plist%set('node-set-ids', pack(node_set_ids, mask=(node_set_ids/=NULL_I)))
        else
          ! faceset required
          if (given_node_sets) &
              call TLS_fatal(label // ': NODE_SET_IDS specified but incompatible with given type')
          if (.not.given_face_sets) &
              call TLS_fatal(label // ': FACE_SET_IDS not specified')
        end if
        if (given_face_sets) &
            call plist%set('face-set-ids', pack(face_set_ids, mask=(face_set_ids/=NULL_I)))
      end block
    end do

  end subroutine read_solid_mechanics_bc_namelists

end module solid_mechanics_bc_namelist
