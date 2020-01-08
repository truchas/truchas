!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module body_factories

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use body_class
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

  subroutine alloc_plane_body(r, n, p, fill_inside)
    use plane_body_type
    class(body), allocatable, intent(out) :: r
    real(r8), intent(in) :: n(:), p
    logical, intent(in) :: fill_inside
    allocate(r, source=plane_body(n, p, fill_inside))
  end subroutine alloc_plane_body


  subroutine alloc_box_body(r, upper, lower, fill_inside)
    use box_body_type
    class(body), allocatable, intent(out) :: r
    real(r8), intent(in) :: upper(:), lower(:)
    logical, intent(in) :: fill_inside
    allocate(r, source=box_body(upper, lower, fill_inside))
  end subroutine alloc_box_body


  subroutine alloc_sphere_body(r, xc, radius, fill_inside)
    use sphere_body_type
    class(body), allocatable, intent(out) :: r
    real(r8), intent(in) :: xc(:), radius
    logical, intent(in) :: fill_inside
    allocate(r, source=sphere_body(xc, radius, fill_inside))
  end subroutine alloc_sphere_body


  subroutine alloc_cylinder_body(r, xc, axis, length, radius, fill_inside)
    use cylinder_body_type
    class(body), allocatable, intent(out) :: r
    real(r8), intent(in) :: xc(:), axis(:), length, radius
    logical, intent(in) :: fill_inside
    allocate(r, source=cylinder_body(xc, axis, length, radius, fill_inside))
  end subroutine alloc_cylinder_body


  subroutine alloc_ellipsoid_body(r, xc, axes, fill_inside)
    use ellipsoid_body_type
    class(body), allocatable, intent(out) :: r
    real(r8), intent(in) :: xc(:), axes(:)
    logical, intent(in) :: fill_inside
    allocate(r, source=ellipsoid_body(xc, axes, fill_inside))
  end subroutine alloc_ellipsoid_body


  subroutine alloc_ellipse_body(r, xc, coeffs, fill_inside)
    use ellipse_body_type
    class(body), allocatable, intent(out) :: r
    real(r8), intent(in) :: xc(:), coeffs(:)
    logical, intent(in) :: fill_inside
    allocate(r, source=ellipse_body(xc, coeffs, fill_inside))
  end subroutine alloc_ellipse_body


  ! subroutine alloc_sinusoid_body(r, coeffs)
  !   use sinusoid_body_type
  !   class(body), allocatable, intent(out) :: r
  !   real(r8), intent(in) :: coeffs(:)
  !   allocate(r, source=sinusoid_body(coeffs))
  ! end subroutine alloc_sinusoid_body


  subroutine alloc_element_block_body(r, mesh, cblockids, fill_inside)
    use element_block_body_type
    use unstr_mesh_type
    class(body), allocatable, intent(out) :: r
    type(unstr_mesh), target, intent(in) :: mesh
    integer, intent(in) :: cblockids(:)
    logical, intent(in) :: fill_inside
    allocate(r, source=element_block_body(mesh, cblockids, fill_inside))
  end subroutine alloc_element_block_body


  subroutine alloc_background_body(r)
    use background_body_type
    class(body), allocatable, intent(out) :: r
    allocate(background_body :: r)
  end subroutine alloc_background_body


  subroutine alloc_body(mesh, params, r, stat, errmsg)

    use unstr_mesh_type
    use parameter_list_type
    use cell_geometry, only: normalized

    type(unstr_mesh), target, intent(in) :: mesh
    type(parameter_list), intent(inout) :: params
    class(body), allocatable, intent(out) :: r
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8), allocatable :: x(:), coeffs(:)
    real(r8) :: p, l
    character(:), allocatable :: rtype, context
    integer :: i
    integer, allocatable :: ids(:)
    logical :: flag

    stat = 0
    context = 'processing ' // params%name() // ': '
    call params%get('type', rtype, stat=stat, errmsg=errmsg)
    if (stat /= 0) return

    select case (rtype)
    case ('sphere')
      call params%get('center', x, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      call params%get('radius', p, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      call params%get('fill-inside', flag, default=.true., stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      ASSERT(size(x)==3)
      ASSERT(p > 0)
      call alloc_sphere_body(r, x, p, flag)

    ! case ('sinusoid')
    !   call params%get('coeffs', x, stat=stat, errmsg=errmsg)
    !   if (stat /= 0) return
    !   ASSERT(size(x)==7)
    !   call alloc_sinusoid_body(r, x)

    case ('plane')
      call params%get('normal', coeffs, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      ASSERT(size(coeffs)==3)
      ASSERT(any(coeffs /= 0))
      coeffs = normalized(coeffs) ! normalize the user provided normal direction
      call params%get('point-on-plane', x, stat=stat)
      call params%get('plane-const', p, stat=i, errmsg=errmsg)
      if (stat == 0) then
        if (i == 0) then
          stat = 1
          errmsg = context // 'cannot define both point-on-plane and plane-const for plane body'
          return
        end if
        ASSERT(size(x)==3)
        p = dot_product(coeffs,x)
      else if (i /= 0) then
        stat = i
        errmsg = context // 'must define either point-on-plane or plane-const for plane body'
        return
      end if
      call params%get('fill-inside', flag, default=.true., stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      call alloc_plane_body(r, coeffs, p, flag)

    case ('box')
      call params%get('upper-corner', x, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      call params%get('lower-corner', coeffs, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      call params%get('fill-inside', flag, default=.true., stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      ASSERT(size(x)==3)
      ASSERT(size(coeffs)==3)
      call alloc_box_body(r, x, coeffs, flag)

    case ('cylinder')
      call params%get('center', x, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      call params%get('axis', coeffs, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      call params%get('radius', p, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      call params%get('length', l, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      call params%get('fill-inside', flag, default=.true., stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      ASSERT(size(x)==3)
      ASSERT(size(coeffs)==3)
      ASSERT(p > 0 .and. l > 0)
      call alloc_cylinder_body(r, x, coeffs, l, p, flag)

    case ('ellipsoid')
      call params%get('center', x, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      call params%get('coeffs', coeffs, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      call params%get('fill-inside', flag, default=.true., stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      ASSERT(size(x)==3)
      ASSERT(all(coeffs > 0))
      call alloc_ellipsoid_body(r, x, coeffs, flag)

    case ('ellipse')
      call params%get('center', x, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      call params%get('coeffs', coeffs, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      call params%get('fill-inside', flag, default=.true., stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      call alloc_ellipse_body(r, x, coeffs, flag)

    case ('element-block')
      call params%get('blockids', ids, stat=stat, errmsg=errmsg)
      print *, 'ids: ', ids
      if (stat /= 0) return
      call params%get('fill-inside', flag, default=.true., stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      call alloc_element_block_body(r, mesh, ids, flag)

    case ('background')
      call alloc_background_body(r)

    case default
      stat = 1
      errmsg = context // 'unknown "type" value: ' // rtype
      return
    end select
    ASSERT(allocated(r))

  end subroutine alloc_body

end module body_factories
