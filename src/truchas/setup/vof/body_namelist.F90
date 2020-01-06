!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module body_namelist

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  implicit none
  private

  public :: read_body_namelists
  type(parameter_list), public :: bodies_params

contains

  subroutine read_body_namelists(lun)

    use parallel_communication, only: is_IOP, broadcast
    use input_utilities
    use string_utilities, only: i_to_c
    use truchas_logging_services

    integer, intent(in) :: lun

    integer :: ios, n
    logical :: found
    character(80) :: iom, pstr
    type(parameter_list), pointer :: plist

    real(r8) :: x(3)

    !! Namelist variables
    character(64) :: axis, surface_name
    real(r8) :: height, length(3), radius, rotation_angle(3), rotation_pt(3), translation_pt(3)

    ! These initial values should go elsewhere
    character(64) :: temperature_function
    integer :: material_number, mesh_material_number
    real(r8) :: velocity(3), temperature, phi(5)

    ! These are unused here and should be removed completely
    character(64) :: tabular_type, fill, file_name
    real(r8) :: Rtheta_tabular_pt(50), RZ_tabular_pt(50)
    
    namelist /body/ surface_name, axis, height, radius, length, &
        rotation_angle, rotation_pt, translation_pt, &
        material_number, phi, temperature, temperature_function, velocity, mesh_material_number, &
        rtheta_tabular_pt, rz_tabular_pt, tabular_type, file_name, fill

    call TLS_info('')
    call TLS_info('Reading BODY namelists ...')

    n = 0
    if (is_IOP) rewind(lun)
    do
      !! Locate the BODY namelist
      if (is_IOP) call seek_to_namelist(lun, 'body', found, iostat=ios)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
      call broadcast(found)
      if (.not.found) exit
      n = n + 1

      !! Default values
      surface_name = NULL_C
      fill = NULL_C
      axis = NULL_C
      height = NULL_R
      length = NULL_R
      radius = NULL_R
      material_number = NULL_I
      mesh_material_number = NULL_I
      rotation_angle = 0
      rotation_pt = 0
      translation_pt = 0

      !! Read the BODY namelist
      if (is_IOP) read(lun, nml=body, iostat=ios, iomsg=iom)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading BODY namelist: ' // trim(iom))

      !! Broadcast the namelist variables
      call broadcast(surface_name)
      call broadcast(axis)
      call broadcast(fill)
      call broadcast(height)
      call broadcast(length)
      call broadcast(radius)
      call broadcast(material_number)
      call broadcast(mesh_material_number)
      call broadcast(rotation_angle)
      call broadcast(rotation_pt)
      call broadcast(translation_pt)

      !! Check values and stuff into a parameter list for later use.
      write(pstr, '(a,i0)') 'body', n
      plist => bodies_params%sublist(trim(pstr))
      select case (trim(surface_name))
      case ('background', 'fill')
        call plist%set('type', 'background')

      case ('plane')
        call plist%set('type', trim(surface_name))
        !call plist%set('point-on-plane', reverse_transform(translation_pt, rotation_angle, rotation_pt))
        call plist%set('point-on-plane', translation_pt)

        ! Set the appropriate plane normal direction based on fill value.
        if (fill == NULL_C .or. trim(fill) == 'inside') then
          call plist%set('normal', -normal_vector(axis, rotation_angle, rotation_pt))
        else if (fill /= NULL_C .and. trim(fill) == 'outside') then
          call plist%set('normal', normal_vector(axis, rotation_angle, rotation_pt))
        else
          call TLS_fatal("error reading BODY namelist: Unexpected fill value. Expected 'inside' or 'outside'.")
        end if

      case ('box')
        call plist%set('type', trim(surface_name))
        call plist%set('upper', translation_pt + length / 2)
        call plist%set('lower', translation_pt - length / 2)

      case ('sphere')
        call plist%set('type', trim(surface_name))
        call plist%set('center', translation_pt)
        call plist%set('radius', radius)

      case ('from mesh file', 'element-block')
        call plist%set('type', 'element-block')
        call plist%set('blockid', mesh_material_number)

      case ('cylinder')
        ! old input centers on cylinder base, new centers on center
        x = normal_vector(axis, rotation_angle, rotation_pt)
        translation_pt = translation_pt + x * height / 2
        call plist%set('type', trim(surface_name))
        call plist%set('center', translation_pt)
        call plist%set('radius', radius)
        call plist%set('length', height)
        call plist%set('axis', x)
        
        ! If fill is unset or set to 'inside', we use the default value in our
        ! constructors (fill_inside = .true.). If set to 'outside', set
        ! fill_inside = .false. Otherwise, return an error.
        if (fill /= NULL_C .and. trim(fill) == 'outside') then
          call plist%set('fill_inside', .false.)
        else if (fill /= NULL_C .and. trim(fill) /= 'inside') then
          call TLS_fatal("error reading BODY namelist: Unexpected fill value. Expected 'inside' or 'outside'.")
        end if
        

      ! case ('ellipsoid')
      !   call plist%set('type', trim(surface_name))
      !   call plist%set('center', translation_pt)
      !   call plist%set('axes', length) ! TODO

      ! case ('box')
      !   call plist%set('type', trim(surface_name))

      ! case ('cone')
      !   call plist%set('type', trim(surface_name))

      case default
        call TLS_fatal('error reading BODY namelist: invalid surface_name')
      end select

      !! TODO: make an initial materials list
    end do

  end subroutine read_body_namelists


  !! Need to convert from axis + rotation to normal direction
  !! Rotation angles are provided in degrees.
  function normal_vector(axis, rotation_angle, rotation_point) result(normal)

    use constants_module, only: pi

    character(*), intent(in) :: axis
    real(r8), intent(in) :: rotation_angle(:), rotation_point(:)
    real(r8) :: normal(3)
    
    real(r8), parameter :: deg = pi / 180
    integer, parameter :: ndim = 3
    integer, parameter :: nrot = 3
    real(r8) :: c, s, tmp, rot(nrot)
    integer :: n, n1, n2, coeff, na

    INSIST(axis == 'x' .or. axis == 'y' .or. axis == 'z')
    normal = 0
    if (axis == 'x') normal(1) = 1
    if (axis == 'y') normal(2) = 1
    if (axis == 'z') normal(3) = 1
    normal = reverse_transform(normal, rotation_angle, rotation_point)
    normal = normal / norm2(normal)

  end function normal_vector


  
  function reverse_transform(q, rotation_angle, rotation_point) result(x)

    use constants_module, only: pi

    real(r8), intent(in) :: q(:), rotation_angle(:), rotation_point(:)
    real(r8) :: x(3)

    real(r8), parameter :: deg = pi / 180
    integer, parameter :: ndim = 3
    integer, parameter :: nrot = 3
    real(r8) :: c, s, tmp, rot(nrot)
    integer :: n, n1, n2, coeff, na

    x = q

    !! The legacy code worked by rotating input points to a coordinate system
    !! where a normal was axis-aligned. The new system works with normal vectors
    !! in the existing coordinate space, so we need to rotate in the opposite
    !! direction and order to get the correct surface. Here we also convert
    !! from degrees to radians.
    INSIST(size(rot) == size(rotation_angle))
    rot = deg * rotation_angle

    do n = nrot, 1, -1
      na = merge(3, n, ndim == 2)
      coeff = merge(-1, 1, na == 2)
      n1 = modulo(na, 3) + 1
      n2 = modulo(na+1, 3) + 1

      if (rot(n) /= 0) then
        c = cos(rot(n))
        s = sin(rot(n))
        tmp =   c*(x(n1)-rotation_point(n1)) - coeff*s*(x(n2)-rotation_point(n2)) + rotation_point(n1)
        x(n2) = c*(x(n2)-rotation_point(n2)) + coeff*s*(x(n1)-rotation_point(n1)) + rotation_point(n2)
        x(n1) = tmp
      end if
    end do

  end function reverse_transform

end module body_namelist