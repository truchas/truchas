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
    complex(r8) :: alpha, g(3)
    character(32) :: g_func
    real(r8) :: sigma

    namelist /electromagnetic_bc/ name, type, face_set_ids, alpha, g, g_func, sigma

    !! Waveguide port feed BC parameters
    real(r8) :: center(3), x_axis(3), y_axis(3), x_width, y_width, power, e_mag
    namelist /electromagnetic_bc/ center, x_axis, y_axis, x_width, y_width, power, e_mag

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
      alpha = NULL_R
      g = NULL_R
      g_func = NULL_C
      sigma = NULL_R

      center = NULL_R
      x_axis = NULL_R
      y_axis = NULL_R
      x_width = NULL_R
      y_width = NULL_R
      power = NULL_R
      e_mag = NULL_R

      if (is_IOP) read(lun, nml=electromagnetic_bc, iostat=ios, iomsg=iom)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading ' // label // ' namelist: ' // trim(iom))

      call broadcast(name)
      call broadcast(type)
      call broadcast(face_set_ids)
      call broadcast(alpha)
      call broadcast(g)
      call broadcast(g_func)
      call broadcast(sigma)

      call broadcast(center)
      call broadcast(x_axis)
      call broadcast(y_axis)
      call broadcast(x_width)
      call broadcast(y_width)
      call broadcast(power)
      call broadcast(e_mag)

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
      case ('pec', 'pmc')
        ! no parameters
      case ('ih-hfield')
        ! no parameters
      case ('wg-port')
        if (all(center == NULL_R)) call TLS_fatal(label // ': CENTER not specified')
        if (any(center == NULL_R)) call TLS_fatal(label // ': 3-vector CENTER not completely specified')
        call plist%set('center', center)
        if (all(x_axis == NULL_R)) then
          call TLS_fatal(label // ': X_AXIS not specified')
        else if (any(x_axis == NULL_R)) then
          call TLS_fatal(label // ': 3-vector X_AXIS not completely specified')
        else
          call plist%set('x-axis', x_axis)
        end if
        if (all(y_axis == NULL_R)) then
          call TLS_fatal(label // ': Y_AXIS not specified')
        else if (any(y_axis == NULL_R)) then
          call TLS_fatal(label // ': 3-vector Y_AXIS not completely specified')
        else
          call plist%set('y-axis', y_axis)
        end if
        if (x_width == NULL_R) then
          call TLS_fatal(label // ': X_WIDTH not specified')
        else if (x_width <= 0) then
          call TLS_fatal(label // ': X_WIDTH <= 0.0')
        else
          call plist%set('x-width', x_width)
        end if
        if (y_width == NULL_R) then
          call TLS_fatal(label // ': Y_WIDTH not specified')
        else if (y_width <= 0) then
          call TLS_fatal(label // ': Y_WIDTH <= 0.0')
        else
          call plist%set('y-width', y_width)
        end if
        if (power /= NULL_R) then
          call plist%set('power', power)
          if (e_mag /= NULL_R) call TLS_fatal(label // ': POWER and E_MAG both specified')
        else
          if (e_mag /= NULL_R) call plist%set('e-mag', e_mag)
        end if
      case ('robin')
        if (alpha == NULL_R) then
          call TLS_fatal(label // ': ALPHA not specified')
        else
          call plist%set('alpha', alpha)
        end if
        if (any(g /= NULL_R) .and. g_func /= NULL_C) then
          call TLS_fatal(label // ': both G and G_FUNC specified')
        else if (any(g /= NULL_R)) then
          if (any(g == NULL_R)) then
            call TLS_fatal(label // ': G not fully specified')
          else
            call plist%set('g', g)
          end if
        else if (g_func /= NULL_C) then
          call plist%set('g', g_func)
        else
          call TLS_fatal(label // ': neither G nor G_FUNC specified')
        end if
      case ('impedance')
        if (sigma == NULL_R) then
          call TLS_fatal(label // ': SIGMA not specified')
        else
          call plist%set('sigma', sigma)
        end if
      case default
        call TLS_fatal(label // ': unknown TYPE: ' // trim(type))
      end select
      call plist%set('type', trim(type))

      call TLS_info('  read namelist "' // trim(name) // '"')
    end do

    if (n == 0) call TLS_info('  none found')

  end subroutine read_electromagnetic_bc_namelists

end module electromagnetic_bc_namelist
