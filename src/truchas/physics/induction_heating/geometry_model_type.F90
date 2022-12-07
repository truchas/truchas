!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module geometry_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  real(r8), parameter, public :: XHAT(3) = [1.0_r8, 0.0_r8, 0.0_r8]
  real(r8), parameter, public :: YHAT(3) = [0.0_r8, 1.0_r8, 0.0_r8]
  real(r8), parameter, public :: ZHAT(3) = [0.0_r8, 0.0_r8, 1.0_r8]

  type, public :: geometry_model
    private
    type(surface_box), allocatable :: array(:)  ! object array
    integer :: nsurf = 0        ! Number of defined surfaces
    integer :: init_size = 10   ! Initial array size
    integer :: incr_size = 5    ! Size increment when resizing array
    real(r8) :: tol = 1.0e-3
  contains
    procedure :: add_plane, add_cylinder, add_cone
    procedure :: surface_exists
    procedure :: tune
    procedure :: get_on_surface_list
    generic :: on_surface => on_surface_by_id_one, on_surface_by_id_many
    procedure, private :: on_surface_by_id_one, on_surface_by_id_many
    procedure :: co_oriented
   procedure, private :: new_index
  end type geometry_model

  type :: surface_box
    class(surface), allocatable :: surf
  end type

  type, abstract :: surface
  contains
    procedure(on_surface), deferred :: on_surface
    procedure(same_nrml), deferred :: co_oriented
  end type

  abstract interface
    pure logical function on_surface(this, point, tol)
      import surface, r8
      class(surface), intent(in) :: this
      real(r8), intent(in) :: point(:), tol
    end function
    pure logical function same_nrml(this, face)
      import surface, r8
      class(surface), intent(in) :: this
      real(r8), intent(in) :: face(:,:)
    end function
  end interface

  type, extends(surface) :: plane_surface
    real(r8) :: point(3)   ! point on plane
    real(r8) :: normal(3)  ! unit normal to plane
  contains
    procedure :: on_surface => on_plane
    procedure :: co_oriented => co_oriented_plane
  end type

  type, extends(surface) :: cylinder_surface
    real(r8) :: point(3) ! point on axis of cylinder
    real(r8) :: axis(3)  ! unit axis vector
    real(r8) :: radius   ! cylinder radius
  contains
    procedure :: on_surface => on_cylinder
    procedure :: co_oriented => co_oriented_cylinder
  end type

  type, extends(surface) :: cone_surface
    real(r8) :: vertex(3)  ! vertex of the cone
    real(r8) :: axis(3)    ! unit axis vector
    real(r8) :: slope      ! radius/axial distance from vertex
  contains
    procedure :: on_surface => on_cone
    procedure :: co_oriented => co_oriented_cone
  end type

contains

  !! Redefines internal parameters.
  subroutine tune(this, tol, init_size, incr_size)
    class(geometry_model), intent(inout) :: this
    real(r8), intent(in), optional :: tol
    integer, intent(in), optional :: init_size, incr_size
    if (present(tol)) then
      if (tol > 0.0_r8) this%tol = tol
    end if
    if (present(init_size)) then
      if (init_size > 0) this%init_size = init_size
    end if
    if (present(incr_size)) then
      if (incr_size > 0) this%incr_size = incr_size
    end if
  end subroutine

  !! Returns true of a surface with the given ID exists; otherwise false.
  elemental logical function surface_exists(this, id)
    class(geometry_model), intent(in) :: this
    integer, intent(in) :: id
    surface_exists = (id <= this%nsurf)
  end function

  !! Add the specified plane to the geometry model
  subroutine add_plane(this, point, normal, id)
    use cell_geometry, only: normalized
    class(geometry_model), intent(inout) :: this
    real(r8), intent(in) :: point(:), normal(:)
    integer, intent(out) :: id
    type(plane_surface), allocatable :: surf
    allocate(surf)
    surf%point = point
    surf%normal = normalized(normal)
    id = this%new_index()
    call move_alloc(surf, this%array(id)%surf)
  end subroutine

  !! Add the specified cylinder to the geometry model
  subroutine add_cylinder(this, point, axis, radius, id)
    use cell_geometry, only: normalized
    class(geometry_model), intent(inout) :: this
    real(r8), intent(in) :: point(:), axis(:), radius
    integer, intent(out) :: id
    type(cylinder_surface), allocatable :: surf
    allocate(surf)
    surf%point = point
    surf%axis = normalized(axis)
    surf%radius = radius
    id = this%new_index()
    call move_alloc(surf, this%array(id)%surf)
  end subroutine

  !! Add the specified cone to the geometry model
  subroutine add_cone(this, vertex, axis, slope, id)
    use cell_geometry, only: normalized
    class(geometry_model), intent(inout) :: this
    real(r8), intent(in) :: vertex(:), axis(:), slope
    integer, intent(out) :: id
    type(cone_surface), allocatable :: surf
    allocate(surf)
    surf%vertex = vertex
    surf%axis = normalized(axis)
    surf%slope = slope
    id = this%new_index()
    call move_alloc(surf, this%array(id)%surf)
  end subroutine

  !! Return true if POINT lies on the surface with ID; otherwise false.
  logical function on_surface_by_id_one(this, id, point)
    class(geometry_model), intent(in) :: this
    integer, intent(in) :: id
    real(r8), intent(in) :: point(:)
    ASSERT(size(point) == 3)
    if (this%surface_exists(id)) then
      on_surface_by_id_one = this%array(id)%surf%on_surface(point, this%tol)
    else
      on_surface_by_id_one = .false.
    end if
  end function

  !! Return true if every point in the array POINT lies on the surface with ID; otherwise false.
  logical function on_surface_by_id_many(this, id, point)
    class(geometry_model), intent(in) :: this
    integer, intent(in) :: id
    real(r8), intent(in) :: point(:,:)
    ASSERT( size(point,dim=1) == 3 )
    if (this%surface_exists(id)) then
      on_surface_by_id_many = on_surface_many(this%array(id)%surf, point, this%tol)
    else
      on_surface_by_id_many = .false.
    end if
  end function

  !! Returns a LIST of surface IDs that contain every point in the array POINT of points.
  subroutine get_on_surface_list(this, point, list)
    class(geometry_model), intent(in) :: this
    real(r8), intent(in) :: point(:,:)
    integer, allocatable, intent(out) :: list(:)
    integer :: n, list_size, l(this%nsurf)
    list_size = 0
    do n = 1, this%nsurf
      if (on_surface_many(this%array(n)%surf, point, this%tol)) then
        list_size = list_size + 1
        l(list_size) = n
      end if
    end do
    list = l(1:list_size)
  end subroutine

  !! Auxiliary function that returns true if every point in the array POINT
  !! lies on the surface SURF; otherwise returns false.
  logical function on_surface_many(surf, point, tol)
    class(surface), intent(in) :: surf
    real(r8), intent(in) :: point(:,:)
    real(r8), intent(in) :: tol
    integer :: n
    on_surface_many = .false.
    do n = 1, size(point,dim=2)
      if (.not.surf%on_surface(point(:,n), tol)) return
    end do
    on_surface_many = .true.
  end function

  !! Returns true if the given oriented FACE has the same orientation as the
  !! surface with the given ID; otherwise it returns false. This assumes that
  !! the face lies on the surface (within tolerance as indicated by on_surface).

  logical function co_oriented(this, id, face)
    class(geometry_model), intent(in) :: this
    integer, intent(in) :: id
    real(r8), intent(in) :: face(:,:)
    co_oriented = this%array(id)%surf%co_oriented(face)
  end function

  !! Auxiliary function that returns a new surface index, automatically resizing
  !! the surface array when necessary.
  integer function new_index(this) result(n)
    class(geometry_model), intent(inout) :: this
    this%nsurf = this%nsurf + 1
    if (allocated(this%array)) then
      if (this%nsurf > size(this%array)) then ! resize THIS%ARRAY
        block
          type(surface_box), allocatable :: old_array(:)
          integer :: j
          call move_alloc(this%array, old_array)
          allocate(this%array(size(old_array)+this%incr_size))
          do j = 1, size(old_array)
            call move_alloc(old_array(j)%surf, this%array(j)%surf)
          end do
        end block
      end if
    else
      allocate(this%array(this%init_size))
    end if
    n = this%nsurf
  end function

!!!! PLANE SURFACE TYPE BOUND PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure logical function on_plane(this, point, tol)
    class(plane_surface), intent(in) :: this
    real(r8), intent(in) :: point(:), tol
    on_plane = abs(dot_product(this%normal,point-this%point)) < tol
  end function

  pure logical function co_oriented_plane(this, face) result(same)
    use cell_geometry, only: face_normal
    class(plane_surface), intent(in) :: this
    real(r8), intent(in) :: face(:,:)
    same = (dot_product(face_normal(face), this%normal) > 0)
  end function

!!!! CYLINDER SURFACE TYPE BOUND PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure logical function on_cylinder (this, point, tol)
    class(cylinder_surface), intent(in) :: this
    real(r8), intent(in) :: point(:), tol
    real :: dist
    dist = sqrt(sum((point - this%point - dot_product(this%axis, point - this%point) * this%axis)**2))
    on_cylinder = abs(dist - this%radius) < tol
  end function

  pure logical function co_oriented_cylinder(this, face) result(same)
    use cell_geometry, only: face_normal, polygon_center
    class(cylinder_surface), intent(in) :: this
    real(r8), intent(in) :: face(:,:)
    real(r8) :: p(3)
    p = polygon_center(face) - this%point
    p = p - dot_product(p, this%axis) * this%axis
    same = (dot_product(face_normal(face), p) > 0)
  end function

!!!! CONE SURFACE TYPE BOUND PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure logical function on_cone(this, point, tol)
    class(cone_surface), intent(in) :: this
    real(r8), intent(in) :: point(:), tol
    real(r8) :: dist, z(3)
    z = point - this%vertex
    dist = sqrt(sum((z - dot_product(this%axis, z) * this%axis)**2))
    on_cone = abs(dist - this%slope*abs(dot_product(this%axis, z))) < tol*sqrt(1.0_r8+this%slope**2)
  end function

  pure logical function co_oriented_cone(this, face) result(same)
    use cell_geometry, only: face_normal, polygon_center
    class(cone_surface), intent(in) :: this
    real(r8), intent(in) :: face(:,:)
    real(r8) :: p(3)
    p = polygon_center(face) - this%vertex
    p = p - dot_product(p, this%axis) * this%axis
    same = (dot_product(face_normal(face), p) > 0)
  end function

end module geometry_model_type
