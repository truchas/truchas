!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module GeometricModeler

  !use kind_parameters
  implicit none
  private
  
  public :: AddPlane, AddCylinder, AddCone, DestroyGeometricModel, GMTune
  public :: SurfaceExists, OnSurface
  public :: dump_geometric_model
  
  public :: same_normal ! added 10 Jul 2003, NNC
  
  integer, parameter :: r8 = selected_real_kind(10,50)  ! IEEE 64 bit floating point
  
  type, public :: GeometricModel
    private
    type(gm_surf), pointer :: surface(:) => null()  ! object array
    integer,       pointer :: link(:) => null()     ! link array
    integer :: nsurf = 0        ! Number of defined surfaces
    integer :: first = 0        ! Array index of the first object
    integer :: next = 0         ! First free array index
    integer :: init_size = 10   ! Initial array size
    integer :: incr_size = 5    ! Size increment when resizing array
    real(kind=r8) :: tol = 1.0e-3
  end type GeometricModel
  
  integer, parameter, private :: GMSURFNULL = 0, GMSURFPLANE = 1, GMSURFCYLINDER = 2, GMSURFCONE = 3
  
  type, private :: gm_surf
    integer :: type = GMSURFNULL
    type(gm_plane), pointer :: plane => null()
    type(gm_cylinder), pointer :: cylinder => null()
    type(gm_cone), pointer :: cone => null()
    !type(gm_sphere), pointer :: sphere => null()
  end type gm_surf
  
  type, private :: gm_plane
    real(kind=r8) :: point(3)   ! point on plane
    real(kind=r8) :: normal(3)  ! unit normal to plane
  end type gm_plane
  
  type, private :: gm_cylinder
    real(kind=r8) :: point(3) ! point on axis of cylinder
    real(kind=r8) :: axis(3)  ! unit axis vector
    real(kind=r8) :: radius   ! cylinder radius
  end type gm_cylinder
  
  type, private :: gm_cone
    real(kind=r8) :: vertex(3)  ! vertex of the cone
    real(kind=r8) :: axis(3)    ! unit axis vector
    real(kind=r8) :: slope      ! radius/axial distance from vertex
  end type gm_cone
  
  interface OnSurface
    module procedure on_surface_by_id_one, on_surface_by_id_many, on_surface_list_many
  end interface
  
  real(kind=r8), parameter, public :: XHAT(3) = (/ 1.0_r8, 0.0_r8, 0.0_r8 /)
  real(kind=r8), parameter, public :: YHAT(3) = (/ 0.0_r8, 1.0_r8, 0.0_r8 /)
  real(kind=r8), parameter, public :: ZHAT(3) = (/ 0.0_r8, 0.0_r8, 1.0_r8 /)
  
contains

  subroutine GMTune (gm, tol, init_size, incr_size)
    type(GeometricModel), intent(inout) :: gm
    real(kind=r8), intent(in), optional :: tol
    integer, intent(in), optional :: init_size, incr_size
    if (present(tol)) then
      if (tol > 0.0_r8) gm%tol = tol
    end if
    if (present(init_size)) then
      if (init_size > 0) gm%init_size = init_size
    end if
    if (present(incr_size)) then
      if (incr_size > 0) gm%incr_size = incr_size
    end if
  end subroutine GMTune

  elemental logical function SurfaceExists (gm, id)
    type(GeometricModel), intent(in) :: gm
    integer, intent(in) :: id
    integer :: n
    n = gm%first
    SurfaceExists = .false.
    do while (n /= 0 .and. id /= n)
      n = gm%link(n)
    end do
    if (n /= 0) SurfaceExists = .true.
  end function SurfaceExists
  
  
  subroutine AddPlane (gm, gmo, point, normal)
  
    type(GeometricModel), intent(inout) :: gm
    integer, intent(out) :: gmo
    real(kind=r8), intent(in) :: point(:)
    real(kind=r8), intent(in) :: normal(:)
    
    integer :: n
    
    call new_index (gm, n)
    gm%nsurf = gm%nsurf + 1
    gm%surface(n)%type = GMSURFPLANE
    allocate(gm%surface(n)%plane)
    gm%surface(n)%plane%point = point
    gm%surface(n)%plane%normal = normalize(normal)
    gmo = n
    
  end subroutine AddPlane
  
  subroutine AddCylinder (gm, gmo, point, axis, radius)
  
    type(GeometricModel), intent(inout) :: gm
    integer, intent(out) :: gmo
    real(kind=r8), intent(in) :: point(:), axis(:), radius
    
    integer :: n
    
    ASSERT( size(point) == 3 )
    ASSERT( size(axis) == 3 )
    
    call new_index (gm, n)
    gm%nsurf = gm%nsurf + 1
    gm%surface(n)%type = GMSURFCYLINDER
    allocate(gm%surface(n)%cylinder)
    gm%surface(n)%cylinder%point = point
    gm%surface(n)%cylinder%axis = normalize(axis)
    gm%surface(n)%cylinder%radius = radius
    gmo = n
    
  end subroutine AddCylinder
  
  subroutine AddCone (gm, gmo, vertex, axis, slope)
  
    type(GeometricModel), intent(inout) :: gm
    integer, intent(out) :: gmo
    real(kind=r8), intent(in) :: vertex(:), axis(:), slope
    
    integer :: n
    
    ASSERT( size(vertex) == 3 )
    ASSERT( size(axis) == 3 )
    ASSERT( slope > 0.0_r8 )
    
    call new_index (gm, n)
    gm%nsurf = gm%nsurf + 1
    gm%surface(n)%type = GMSURFCONE
    allocate(gm%surface(n)%cone)
    gm%surface(n)%cone%vertex = vertex
    gm%surface(n)%cone%axis = normalize(axis)
    gm%surface(n)%cone%slope = slope
    gmo = n
    
  end subroutine AddCone
  
  logical function on_surface_by_id_one (gm, id, point)
    type(GeometricModel), intent(in) :: gm
    integer, intent(in) :: id
    real(kind=r8), intent(in) :: point(:)
    ASSERT( size(point) == 3 )
    if (SurfaceExists(gm, id)) then
      on_surface_by_id_one = on_surface_one(gm%surface(id), point, gm%tol)
    else
      on_surface_by_id_one = .false.
    end if
  end function on_surface_by_id_one
  
  logical function on_surface_by_id_many (gm, id, point)
    type(GeometricModel), intent(in) :: gm
    integer, intent(in) :: id
    real(kind=r8), intent(in) :: point(:,:)
    ASSERT( size(point,dim=1) == 3 )
    if (SurfaceExists(gm, id)) then
      on_surface_by_id_many = on_surface_many(gm%surface(id), point, gm%tol)
    else
      on_surface_by_id_many = .false.
    end if
  end function on_surface_by_id_many
  
  function on_surface_list_many (gm, point) result (list)
    type(GeometricModel), intent(in) :: gm
    real(kind=r8), intent(in) :: point(:,:)
    integer, pointer :: list(:)
    integer :: n, list_size, l(gm%nsurf)
    ASSERT( size(point,dim=1) == 3 )
    list_size = 0
    n = gm%first
    do while (n /= 0)
      if (on_surface_many(gm%surface(n), point, gm%tol)) then
        list_size = list_size + 1
        l(list_size) = n
      end if
      n = gm%link(n)
    end do
    allocate(list(list_size))
    if (list_size > 0) list = l(1:list_size)
  end function on_surface_list_many
    
    
  logical function on_surface_one (surf, point, tol)
    type(gm_surf), intent(in) :: surf
    real(kind=r8), intent(in) :: point(:)
    real(kind=r8), intent(in) :: tol
    select case (surf%type)
    case (GMSURFPLANE)
      on_surface_one = on_plane(surf%plane, point, tol)
    case (GMSURFCYLINDER)
      on_surface_one = on_cylinder(surf%cylinder, point, tol)
    case (GMSURFCONE)
      on_surface_one = on_cone(surf%cone, point, tol)
    case default
      on_surface_one = .false.
    end select
  end function on_surface_one
  
  logical function on_surface_many (surf, point, tol)
    type(gm_surf), intent(in) :: surf
    real(kind=r8), intent(in) :: point(:,:)
    real(kind=r8), intent(in) :: tol
    integer :: n
    on_surface_many = .false.
    select case (surf%type)
    case (GMSURFPLANE)
      n = size(point,dim=2)
      do while (on_plane(surf%plane, point(:,n), tol))
        n = n - 1
        if (n == 0) exit
      end do
      if (n == 0) on_surface_many = .true.
    case (GMSURFCYLINDER)
      n = size(point,dim=2)
      do while (on_cylinder(surf%cylinder, point(:,n), tol))
        n = n - 1
        if (n == 0) exit
      end do
      if (n == 0) on_surface_many = .true.
    case (GMSURFCONE)
      n = size(point,dim=2)
      do while (on_cone(surf%cone, point(:,n), tol))
        n = n - 1
        if (n == 0) exit
      end do
      if (n == 0) on_surface_many = .true.
    end select
  end function on_surface_many
  
  logical function on_plane (plane, point, tol)
    type(gm_plane), intent(in) :: plane
    real(kind=r8),  intent(in) :: point(:)
    real(kind=r8),  intent(in) :: tol
    on_plane = abs(dot_product(plane%normal,point-plane%point)) < tol
  end function on_plane
  
  logical function on_cylinder (cyl, point, tol)
    type(gm_cylinder), intent(in) :: cyl
    real(kind=r8), intent(in) :: point(:), tol
    real :: dist
    dist = sqrt(sum((point - cyl%point - dot_product(cyl%axis, point - cyl%point) * cyl%axis)**2))
    on_cylinder = abs(dist - cyl%radius) < tol
  end function on_cylinder
  
  logical function on_cone (cone, point, tol)
    type(gm_cone), intent(in) :: cone
    real(kind=r8), intent(in) :: point(:), tol
    real(kind=r8) :: dist, z(3)
    z = point - cone%vertex
    dist = sqrt(sum((z - dot_product(cone%axis, z) * cone%axis)**2))
    on_cone = abs(dist - cone%slope*abs(dot_product(cone%axis, z))) < tol*sqrt(1.0_r8+cone%slope**2)
  end function on_cone
  
  subroutine new_index (gm, n)
  
    type(GeometricModel), intent(inout) :: gm
    integer, intent(out) :: n
    
    integer :: j, old_size, new_size
    integer, pointer :: old_link(:)
    
    !! Check for free space; resize arrays if necessary.
    if (gm%next == 0) then  ! no free space
      if (associated(gm%link)) then ! increase size
        old_size = size(gm%link)
        new_size = old_size + gm%incr_size
        old_link => gm%link
        allocate(gm%link(new_size))
        gm%link(1:old_size) = old_link
        deallocate(old_link)
      else  ! initial allocatation
        old_size = 0
        new_size = gm%init_size
        allocate(gm%link(new_size))
      end if
      !! Place new indices onto the free list
      gm%next = old_size + 1
      do j = old_size+1, new_size-1
        gm%link(j) = j + 1
      end do
      gm%link(new_size) = 0
      call resize_surf (gm%surface, new_size) ! resize object array
    end if
    
    !! Move index from free list onto used list.
    n = gm%next
    gm%next = gm%link(n)
    gm%link(n) = gm%first
    gm%first = n
    
  end subroutine new_index
  
  
  subroutine resize_surf (array, new_size)
  
    type(gm_surf), pointer :: array(:)
    integer, intent(in) :: new_size
    
    integer :: old_size
    type(gm_surf), pointer :: old_array(:)
    
    if (associated(array)) then
      if (size(array) < new_size) then
        old_size = size(array)
        old_array => array
        allocate(array(new_size))
        array(1:old_size) = old_array
        deallocate(old_array)
      end if
    else
      allocate(array(new_size))
    end if
    
  end subroutine resize_surf
  
  pure function normalize (x) result (y)
  
    real(kind=r8), intent(in) :: x(:)
    real(kind=r8) :: y(3)
    
    integer :: n
    real(kind=r8) :: u, v, r

    n = 1
    if (abs(x(2)) > abs(x(n))) n = 2
    if (abs(x(3)) > abs(x(n))) n = 3

    select case (n)
    case (1)
      u = x(2) / abs(x(1))
      v = x(3) / abs(x(1))
      r = 1.0_r8 / sqrt(1.0_r8 + (u**2 + v**2))
      y(1) = sign(r, x(1))
      y(2) = u * r
      y(3) = v * r
    case (2)
      u = x(1) / abs(x(2))
      v = x(3) / abs(x(2))
      r = 1.0_r8 / sqrt(1.0_r8 + (u**2 + v**2))
      y(1) = u * r
      y(2) = sign(r, x(2))
      y(3) = v * r
    case (3)
      u = x(1) / abs(x(3))
      v = x(2) / abs(x(3))
      r = 1.0_r8 / sqrt(1.0_r8 + (u**2 + v**2))
      y(1) = u * r
      y(2) = v * r
      y(3) = sign(r, x(3))
    end select
    
  end function normalize
  
  subroutine dump_surface (s)
    type(gm_surf), intent(in) :: s
    write(unit=*,fmt='(3x,a,i3)') 'TYPE=', s%type
    if (associated(s%plane)) then
      write(unit=*,fmt='(3x,a)') 'Pointer PLANE is associated:'
      write(unit=*,fmt='(6x,a,3es13.5)') 'PLANE%POINT  =', s%plane%point
      write(unit=*,fmt='(6x,a,3es13.5)') 'PLANE%NORMAL =', s%plane%normal
    else
      write(unit=*,fmt='(a)') 'pointer PLANE is null'
    end if
    if (associated(s%cylinder)) then
      write(unit=*,fmt='(3x,a)') 'Pointer CYLINDER is associated:'
      write(unit=*,fmt='(6x,a,3es13.5)') 'CYLINDER%POINT  =', s%cylinder%point
      write(unit=*,fmt='(6x,a,3es13.5)') 'CYLINDER%AXIS   =', s%cylinder%axis
      write(unit=*,fmt='(6x,a,es13.5)')  'CYLINDER%RADIUS =', s%cylinder%radius
    else
      write(unit=*,fmt='(a)') 'pointer CYLINDER is null'
    end if
    if (associated(s%cone)) then
      write(unit=*,fmt='(3x,a)') 'Pointer CONE is associated:'
      write(unit=*,fmt='(6x,a,3es13.5)') 'CONE%VERTEX =', s%cone%vertex
      write(unit=*,fmt='(6x,a,3es13.5)') 'CONE%AXIS   =', s%cone%axis
      write(unit=*,fmt='(6x,a,es13.5)')  'CONE%SLOPE  =', s%cone%slope
    else
      write(unit=*,fmt='(a)') 'pointer CONE is null'
    end if
  end subroutine dump_surface
  
  subroutine dump_geometric_model (g)
    type(GeometricModel), intent(in) :: g
    integer :: n
    write(unit=*,fmt='(a,es13.5)') 'TOL=', g%tol
    write(unit=*,fmt='(a,i4)') 'INIT_SIZE=', g%init_size
    write(unit=*,fmt='(a,i4)') 'INCR_SIZE=', g%incr_size
    write(unit=*,fmt='(a,i3)') 'NEXT=', g%next
    write(unit=*,fmt='(a,i3)') 'FIRST=', g%first
    if (associated(g%link)) then
      write(unit=*,fmt='(a,10i3,:,/,(5x,10i3))') 'LINK=', g%link
    else
      write(unit=*,fmt='(a)') 'pointer LINK is null'
    end if
    if (associated(g%surface)) then
      n = g%first
      do while (n /= 0)
        write(unit=*,fmt='(a,i3.3,a)') 'SURFACE(', n, '):'
        call dump_surface (g%surface(n))
        n = g%link(n)
      end do
    else
      write(unit=*,fmt='(a)') 'pointer SURFACE is null'
    end if
  end subroutine dump_geometric_model
  
  !!
  !! Added 10 Jul 2003, NNC
  !!
  !! This procedure determines whether a facet is coherently oriented or oppositely
  !! oriented with a surface.  It assumes that the facet lies on the surface (as
  !! indicated by OnSurfaces).  'Coherence' is an ill-defined concept unless the
  !! facet well-approximates the surface (as it will in a reasonable mesh);
  !! nevertheless this routine will always return a result, meaningful or not.
  !! NB: This has been added in great haste; further thought and refinement is
  !! certainly needed.
  !!
  
  logical function same_normal (g, id, face)
  
    type(GeometricModel), intent(in) :: g
    integer, intent(in) :: id
    real(kind=r8), intent(in) :: face(:,:)
    
    real(kind=r8) :: n(3), p(3)
    type(gm_cylinder), pointer :: cyl
    type(gm_cone), pointer :: cone
    
    ASSERT( size(face,dim=1) == 3 )
    
    select case (size(face,dim=2))
    case (3)  ! triangular face
      n = face_normal(face)
    case default
      print *, 'same_normal: PANIC! not implemented for non-triangular faces'
      stop
    end select 
    
    select case (g%surface(id)%type)
    case (GMSURFPLANE)
      same_normal = (dot_product(n, g%surface(id)%plane%normal) > 0)
    case (GMSURFCYLINDER)
      cyl => g%surface(id)%cylinder
      p = sum(face,dim=2) / 3.0_r8 - cyl%point
      p = p - dot_product(p, cyl%axis) * cyl%axis
      same_normal = (dot_product(n, p) > 0)
    case (GMSURFCONE)
      cone => g%surface(id)%cone
      p = sum(face,dim=2) / 3.0_r8 - cone%vertex
      p = p - dot_product(p, cone%axis) * cone%axis
      same_normal = (dot_product(n, p) > 0)
    end select
    
  end function same_normal
  
  !! Copied from MeshSupport
  pure function face_normal (x) result (p)
    real(kind=r8), intent(in) :: x(:,:)
    real(kind=r8) :: p(3)
    p(:) = 0.5_r8 * cross_product(x(:,1) - x(:,2), x(:,2) - x(:,3))
  end function face_normal
  
  !! Copied from MeshSupport
  pure function cross_product (a, b) result (axb)
    real(kind=r8), intent(in) :: a(:), b(:)
    real(kind=r8) :: axb(3)
    axb(1) = a(2) * b(3) - a(3) * b(2)
    axb(2) = a(3) * b(1) - a(1) * b(3)
    axb(3) = a(1) * b(2) - a(2) * b(1)
  end function cross_product

  subroutine DestroyGeometricModel (gm)
    type(GeometricModel), intent(inout) :: gm
    type(GeometricModel) :: default
    integer :: j
    if (associated(gm%surface)) then
      do j = 1, gm%nsurf
        select case (gm%surface(j)%type)
        case (GMSURFPLANE)
          deallocate(gm%surface(j)%plane)
        case (GMSURFCYLINDER)
          deallocate(gm%surface(j)%cylinder)
        case (GMSURFCONE)
          deallocate(gm%surface(j)%cone)
        end select
      end do
      deallocate(gm%surface, gm%link)
    end if
    gm = default  ! set GM to the default initialization values
  end subroutine DestroyGeometricModel

end module GeometricModeler
