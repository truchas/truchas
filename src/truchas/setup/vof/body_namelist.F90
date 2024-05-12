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
    character(64) :: axis, surface_name, fill
    real(r8) :: height, length(3), radius, rotation_angle(3), rotation_pt(3), translation_pt(3)
    integer :: mesh_material_number(16)

    ! These values are read & initialized by body_input_module, not here.
    character(64) :: temperature_function, conc_func(5), material_name
    real(r8) :: velocity(3), temperature, conc(5)

    namelist /body/ surface_name, axis, height, radius, length, fill, &
        rotation_angle, rotation_pt, translation_pt, &
        material_name, conc, conc_func, temperature, temperature_function, velocity, mesh_material_number

    call TLS_info('Reading BODY namelists (second pass) ...')

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
      material_name = NULL_C
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
      call broadcast(material_name)
      call broadcast(mesh_material_number)
      call broadcast(rotation_angle)
      call broadcast(rotation_pt)
      call broadcast(translation_pt)

      ! If fill is unset or set to 'inside', we use the default value in our
      ! constructors (fill-inside = .true.). If set to 'outside', set
      ! fill-inside = .false. Otherwise, return an error.
      if ((trim(surface_name) /= 'background' .and. trim(surface_name) /= 'fill') &
          .and. fill /= NULL_C .and. trim(fill) /= 'inside' .and. trim(fill) /= 'outside') then
        call TLS_fatal("error reading BODY namelist: Unexpected fill value. Expected 'inside' or 'outside'.")
      end if

      !! Check values and stuff into a parameter list for later use.
      write(pstr, '(a,i0)') 'body', n
      plist => bodies_params%sublist(trim(pstr))
      select case (trim(surface_name))
      case ('background', 'fill')
        call plist%set('type', 'background')

      case ('plane')
        call plist%set('type', trim(surface_name))
        call plist%set('normal', -normal_vector(axis, rotation_angle, rotation_pt))
        call plist%set('point-on-plane', translation_pt)
        if (trim(fill) == 'outside') call plist%set('fill-inside', .false.)

      case ('box')
        call plist%set('type', trim(surface_name))
        call plist%set('upper-corner', translation_pt + length / 2)
        call plist%set('lower-corner', translation_pt - length / 2)
        if (trim(fill) == 'outside') call plist%set('fill-inside', .false.)

      case ('sphere')
        call plist%set('type', trim(surface_name))
        call plist%set('center', translation_pt)
        call plist%set('radius', radius)
        if (trim(fill) == 'outside') call plist%set('fill-inside', .false.)

      case ('from mesh file', 'element-block')
        call plist%set('type', 'element-block')
        call plist%set('blockids', mesh_material_number(:findloc(mesh_material_number, NULL_I, dim=1)-1))
        if (trim(fill) == 'outside') call plist%set('fill-inside', .false.)

      case ('cylinder')
        ! old input centers on cylinder base, new centers on center
        x = normal_vector(axis, rotation_angle, rotation_pt)
        translation_pt = translation_pt + x * height / 2
        call plist%set('type', trim(surface_name))
        call plist%set('center', translation_pt)
        call plist%set('radius', radius)
        call plist%set('length', height)
        call plist%set('axis', x)
        if (trim(fill) == 'outside') call plist%set('fill-inside', .false.)
        
      case ('ellipsoid')
        ! This type currently does not support rotation.
        call plist%set('type', trim(surface_name))
        call plist%set('center', translation_pt)
        call plist%set('coeffs', length)
        if (trim(fill) == 'outside') call plist%set('fill-inside', .false.)

      case ('ellipse')
        ! This type currently does not support rotation.
        call plist%set('type', trim(surface_name))
        call plist%set('center', translation_pt)
        call plist%set('coeffs', length(:2))
        if (trim(fill) == 'outside') call plist%set('fill-inside', .false.)

      ! case ('cone')
      !   call plist%set('type', trim(surface_name))

      case default
        call TLS_fatal('error reading BODY namelist: invalid surface_name')
      end select
    end do

    select case (n)
    case (0)
      call TLS_fatal('no BODY namelist found')
    case (1)
      call TLS_info('  read 1 BODY namelist')
    case default
      call TLS_info('  read ' // i_to_c(n) // ' BODY namelists')
    end select

  end subroutine read_body_namelists


  !! Need to convert from axis + rotation to normal direction
  !! Rotation angles are provided in degrees.
  function normal_vector(axis, rotation_angle, rotation_point) result(normal)

    use math_constants, only: pi

    character(*), intent(in) :: axis
    real(r8), intent(in) :: rotation_angle(:), rotation_point(:)
    real(r8) :: normal(3)
    
    real(r8), parameter :: deg = pi / 180
    integer, parameter :: ndim = 3
    integer, parameter :: nrot = 3

    INSIST(axis == 'x' .or. axis == 'y' .or. axis == 'z')
    normal = 0
    if (axis == 'x') normal(1) = 1
    if (axis == 'y') normal(2) = 1
    if (axis == 'z') normal(3) = 1
    normal = reverse_transform(normal, rotation_angle, rotation_point)
    normal = normal / norm2(normal)

  end function normal_vector


  
  function reverse_transform(q, rotation_angle, rotation_point) result(x)

    use math_constants, only: pi

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
