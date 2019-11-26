!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module body_factories

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use body_class
  use truchas_logging_services
  implicit none
  private

  !! higher-level procedure which takes a parameter list as input
  public :: alloc_body

  !! low-level factories
  public :: alloc_plane_body
  public :: alloc_box_body
  public :: alloc_sphere_body
  public :: alloc_ellipsoid_body
  ! public :: alloc_halfsphere_body
  ! public :: alloc_cylinder_body
  ! public :: alloc_ellipse_body
  public :: alloc_element_block_body
  public :: alloc_background_body

contains

  subroutine alloc_plane_body(r, n, p)
    use plane_body_type
    class(body), allocatable, intent(out) :: r
    real(r8), intent(in) :: n(:), p
    allocate(r, source=plane_body(n, p))
  end subroutine alloc_plane_body


  subroutine alloc_box_body(r, upper, lower)
    use box_body_type
    class(body), allocatable, intent(out) :: r
    real(r8), intent(in) :: upper(:), lower(:)
    allocate(r, source=box_body(upper, lower))
  end subroutine alloc_box_body


  subroutine alloc_sphere_body(r, xc, radius)
    use sphere_body_type
    class(body), allocatable, intent(out) :: r
    real(r8), intent(in) :: xc(:), radius
    allocate(r, source=sphere_body(xc, radius))
  end subroutine alloc_sphere_body


  subroutine alloc_ellipsoid_body(r, xc, axes)
    use ellipsoid_body_type
    class(body), allocatable, intent(out) :: r
    real(r8), intent(in) :: xc(:), axes(:)
    allocate(r, source=ellipsoid_body(xc, axes))
  end subroutine alloc_ellipsoid_body


  ! subroutine alloc_halfsphere_body(r, xc, radius, n)
  !   use halfsphere_body_type
  !   class(body), allocatable, intent(out) :: r
  !   real(r8), intent(in) :: xc(:), radius, n(:)
  !   allocate(r, source=halfsphere_body(xc, radius, n))
  ! end subroutine alloc_halfsphere_body


  ! subroutine alloc_cylinder_body(r, xc, axis, radius, halfheight)
  !   use cylinder_body_type
  !   class(body), allocatable, intent(out) :: r
  !   real(r8), intent(in) :: xc(:), axis(:), radius, halfheight
  !   allocate(r, source=cylinder_body(xc, axis, radius, halfheight))
  ! end subroutine alloc_cylinder_body


  ! subroutine alloc_ellipse_body(r, xc, axis, coeffs, halfheight)
  !   use ellipse_body_type
  !   class(body), allocatable, intent(out) :: r
  !   real(r8), intent(in) :: xc(:), axis(:), coeffs(:), halfheight
  !   allocate(r, source=ellipse_body(xc, axis, coeffs, halfheight))
  ! end subroutine alloc_ellipse_body


  ! subroutine alloc_sinusoid_body(r, coeffs)
  !   use sinusoid_body_type
  !   class(body), allocatable, intent(out) :: r
  !   real(r8), intent(in) :: coeffs(:)
  !   allocate(r, source=sinusoid_body(coeffs))
  ! end subroutine alloc_sinusoid_body


  subroutine alloc_element_block_body(r, mesh, cblockid)
    use element_block_body_type
    use unstr_mesh_type
    class(body), allocatable, intent(out) :: r
    type(unstr_mesh), target, intent(in) :: mesh
    integer, intent(in) :: cblockid
    allocate(r, source=element_block_body(mesh, cblockid))
  end subroutine alloc_element_block_body


  subroutine alloc_background_body(r)
    use background_body_type
    class(body), allocatable, intent(out) :: r
    allocate(background_body :: r)
  end subroutine alloc_background_body


  subroutine alloc_body(mesh, params, r)

    use unstr_mesh_type
    use parameter_list_type
    use cell_geometry, only: normalized

    type(unstr_mesh), intent(in) :: mesh
    type(parameter_list), intent(inout) :: params
    class(body), allocatable, intent(out) :: r

    real(r8), allocatable :: x(:), coeffs(:)
    real(r8) :: p
    character(:), allocatable :: rtype, context, errmsg
    integer :: i, stat

    context = 'processing ' // params%name() // ': '
    call params%get('type', rtype, stat=stat, errmsg=errmsg)
    if (stat /= 0) call TLS_fatal(context//errmsg)
    select case (rtype)
    case ('sphere')
      call params%get('center', x, stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
      call params%get('radius', p, stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
      ASSERT(size(x)==3)
      call alloc_sphere_body(r, x, p)
    case ('ellipsoid')
      call params%get('center', x, stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
      call params%get('axes', coeffs, stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
      ASSERT(size(x)==3)
      call alloc_ellipsoid_body(r, x, coeffs)
    ! case ('sinusoid')
    !   call params%get('coeffs', x, stat=stat, errmsg=errmsg)
    !   if (stat /= 0) call TLS_fatal(context//errmsg)
    !   ASSERT(size(x)==7)
    !   call alloc_sinusoid_body(r, x)
    case ('plane')
      call params%get('normal', coeffs, stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
      ASSERT(size(coeffs)==3)
      coeffs = normalized(coeffs) ! normalize the user provided normal direction
      call params%get('point-on-plane', x, stat=stat)
      call params%get('plane-const', p, stat=i, errmsg=errmsg)
      if (stat == 0) then
        ASSERT(size(x)==3)
        if (i == 0) call TLS_fatal('plane body may only have defined either point-on-plane OR plane-const.')
        p = dot_product(coeffs,x)
      else if (i /= 0) then
        call TLS_fatal(context//errmsg)
      end if
      call alloc_plane_body(r, coeffs, p)
    case ('box')
      call params%get('upper', x, stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
      call params%get('lower', coeffs, stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
      ASSERT(size(x)==3)
      ASSERT(size(coeffs)==3)
      call alloc_box_body(r, x, coeffs)
    ! case ('halfsphere', 'half-sphere')
    !   call params%get('normal', normal, stat=stat, errmsg=errmsg)
    !   if (stat /= 0) call TLS_fatal(context//errmsg)
    !   call params%get('center', x, stat=stat, errmsg=errmsg)
    !   if (stat /= 0) call TLS_fatal(context//errmsg)
    !   call params%get('radius', p, stat=stat, errmsg=errmsg)
    !   if (stat /= 0) call TLS_fatal(context//errmsg)
    !   ASSERT(size(x)==3)
    !   ASSERT(size(normal)==3)
    !   normal = normalized(normal)
    !   call alloc_halfsphere_body(r, x, p, normal)
    ! case ('cylinder')
    !   call params%get('center', x, stat=stat, errmsg=errmsg)
    !   if (stat /= 0) call TLS_fatal(context//errmsg)
    !   call params%get('axis', coeffs, stat=stat, errmsg=errmsg)
    !   if (stat /= 0) call TLS_fatal(context//errmsg)
    !   call params%get('radius', p, stat=stat, errmsg=errmsg)
    !   if (stat /= 0) call TLS_fatal(context//errmsg)
    !   call params%get('halfheight', d, stat=stat, errmsg=errmsg)
    !   if (stat /= 0) call TLS_fatal(context//errmsg)
    !   ASSERT(size(x)==3)
    !   ASSERT(size(coeffs)==3)
    !   call alloc_cylinder_body(r, x, coeffs, p, d)
    ! case ('ellipse')
    !   call params%get('center', x, stat=stat, errmsg=errmsg)
    !   if (stat /= 0) call TLS_fatal(context//errmsg)
    !   call params%get('axis', coeffs, stat=stat, errmsg=errmsg)
    !   if (stat /= 0) call TLS_fatal(context//errmsg)
    !   call params%get('coeffs', coeffs, stat=stat, errmsg=errmsg)
    !   if (stat /= 0) call TLS_fatal(context//errmsg)
    !   call params%get('halfheight', d, stat=stat, errmsg=errmsg)
    !   if (stat /= 0) call TLS_fatal(context//errmsg)
    !   call alloc_ellipse_body(r, x, coeffs, coeffs, d)
    case ('element-block')
      call params%get('blockid', i, stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
      call alloc_element_block_body(r, mesh, i)
    case ('background')
      call alloc_background_body(r)
    case default
      call TLS_fatal(context//'unknown "type" value: '//rtype)
    end select
    ASSERT(allocated(r))

  end subroutine alloc_body

end module body_factories
