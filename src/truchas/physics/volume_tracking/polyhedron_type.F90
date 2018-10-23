!!
!! polyhedron_type
!!
!! This module defines an arbitrary polyhedron type, along with routines for
!! calculating volume, splitting polyhedra, locating intersections, etc.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! October 2015
!!
!! References:
!!     1. Hopcroft and Kahn. A Paradigm for Robust Geometric Algorithms. Algorithmica, 1992
!!

#include "f90_assert.fpp"

module polyhedron_type

  use kinds,  only: r8
  use truchas_logging_services
  use polygon_type
  use pure_polyhedron_type
  implicit none
  private

  type, public :: polyhedron
    !private
    type(pure_polyhedron) :: parent
    type(pure_polyhedron), allocatable :: convex_polyhedron(:)
    integer :: nchildren = 0
  contains
    procedure, private :: init_polyhedron
    procedure, private :: init_pure_polyhedron
    procedure, private :: init_polyhedron_null
    procedure, private :: init_tet
    procedure, private :: init_tesselated
    generic   :: init => init_polyhedron, init_polyhedron_null, init_tet, init_pure_polyhedron
    procedure :: volume
    !procedure :: intersection_verts
    procedure :: split
    procedure :: volume_behind_plane
    procedure :: print_data
    procedure :: tesselate
    procedure :: is_inside
    procedure, private :: tesselated_split
    procedure, private :: tesselated_volume_behind_plane
  end type polyhedron

contains

  subroutine init_polyhedron (this, ierr, x, face_v, edge_v, face_normal, vol, tesselate)

    class(polyhedron),  intent(out) :: this
    integer,            intent(out) :: ierr
    real(r8),           intent(in)  :: x(:,:)
    integer,            intent(in)  :: face_v(:,:), edge_v(:,:)
    real(r8), optional, intent(in)  :: face_normal(:,:), vol
    logical, optional, intent(in) :: tesselate

    logical :: tesselateh

    if (present(tesselate)) then
      tesselateh = tesselate
    else
      tesselateh = .true.
    end if

    call this%parent%init(ierr, x, face_v, edge_v, face_normal, vol)
    if (tesselateh) then
      if (this%parent%has_nonplanar_face()) &
          call this%tesselate()
    end if

  end subroutine init_polyhedron

  subroutine init_pure_polyhedron(this, p)
    class(polyhedron), intent(out) :: this
    type(pure_polyhedron), intent(in) :: p
    this%parent = p
  end subroutine init_pure_polyhedron

  subroutine init_tesselated(this, polyhedra)

    class(polyhedron), intent(out) :: this
    type(pure_polyhedron), intent(in) :: polyhedra(:)

    call this%parent%init()
    this%convex_polyhedron = polyhedra
    this%nchildren = size(polyhedra)

  end subroutine init_tesselated

  subroutine init_tet (this, ierr, x, face_normal, vol, set_face_normals)
    class(polyhedron),  intent(out) :: this
    integer,            intent(out) :: ierr
    real(r8),           intent(in)  :: x(:,:)
    real(r8), optional, intent(in)  :: face_normal(:,:), vol
    logical, optional, intent(in) :: set_face_normals
    call this%parent%init(ierr, x, face_normal, vol, set_face_normals)
  end subroutine init_tet

  subroutine init_polyhedron_null (this)
    class(polyhedron), intent(out) :: this
    call this%parent%init()
  end subroutine init_polyhedron_null

  real(r8) function volume(this)

    class(polyhedron), intent(inout) :: this

    integer :: t

    if (this%nchildren > 0) then
      volume = 0
      do t = 1,this%nchildren
        volume = volume + this%convex_polyhedron(t)%volume()
      end do
    else
      volume = this%parent%volume()
    end if

  end function volume

  subroutine tesselate(this)

    class(polyhedron), intent(inout) :: this

    if (this%nchildren == 0) then
      this%convex_polyhedron = this%parent%tesselated()
      this%nchildren = size(this%convex_polyhedron)
    end if

  end subroutine tesselate

  subroutine split(this,P,split_poly,ierr)

    use plane_type

    class(polyhedron), intent(in) :: this
    class(plane),      intent(in) :: P
    type(polyhedron), intent(inout) :: split_poly(:)
    !type(polygon_box), intent(out) :: interface_polygons ! interface polygon
    integer,          intent(out) :: ierr

    type(pure_polyhedron) :: spp(2)

    if (this%nchildren > 0) then
      call this%tesselated_split(P, split_poly, ierr)
    else
      call this%parent%split(P, spp, ierr)
      call split_poly(1)%init(spp(1))
      call split_poly(2)%init(spp(2))
    end if

  end subroutine split

  subroutine tesselated_split (this,P,split_poly,ierr)

    use plane_type

    class(polyhedron), intent(in) :: this
    class(plane),      intent(in) :: P
    type(polyhedron), intent(inout) :: split_poly(:)
    !type(polygon_box), intent(out) :: interface_polygons
    integer,          intent(out) :: ierr

    type(pure_polyhedron) :: tmp(size(this%convex_polyhedron),2)
    integer :: t
    ! type(polygon_box) :: intpoly
    ! type(polygon), allocatable :: tmp_polygon(:)


    ! allocate(interface_polygons%elements(this%nchildren))
    ! interface_polygons%n_elements = 0

    do t = 1,this%nchildren
      call this%convex_polyhedron(t)%split(P, tmp(t,:), ierr)
      ! if (intpoly%n_elements == 1) then
      !   interface_polygons%n_elements = interface_polygons%n_elements + 1
      !   interface_polygons%elements(interface_polygons%n_elements) = intpoly%elements(1)
      ! end if
    end do

    ! ! resize interface_polygons list
    ! tmp_polygon = interface_polygons%elements(:interface_polygons%n_elements)
    ! call move_alloc(tmp_polygon, interface_polygons%elements)

    ! give split sub-polyhedra to split_poly
    call split_poly(1)%init_tesselated(tmp(:,1))
    call split_poly(2)%init_tesselated(tmp(:,2))

  end subroutine tesselated_split

  real(r8) function tesselated_volume_behind_plane(this,P,ierr)

    use plane_type

    class(polyhedron), intent(inout) :: this
    class(plane), intent(in) :: P
    integer, intent(out) :: ierr

    integer :: t

    tesselated_volume_behind_plane = 0
    do t = 1,this%nchildren
      tesselated_volume_behind_plane = tesselated_volume_behind_plane + &
          this%convex_polyhedron(t)%volume_behind_plane(P,ierr)
      if (ierr /= 0) exit
    end do

  end function tesselated_volume_behind_plane

  real(r8) function volume_behind_plane(this, P, ierr)

    use plane_type

    class(polyhedron), intent(inout) :: this
    class(plane),      intent(in) :: P
    integer,           intent(out) :: ierr

    if (this%nchildren > 0) then
      volume_behind_plane = this%tesselated_volume_behind_plane(P,ierr)
    else
      volume_behind_plane = this%parent%volume_behind_plane(P,ierr)
    end if

  end function volume_behind_plane

  ! check if the given point is inside the polyhedron
  logical function is_inside(this, x)

    class(polyhedron), intent(in) :: this
    real(r8), intent(in) :: x(:)

    integer :: t

    if (this%nchildren > 0) then
      is_inside = .false.
      do t = 1,this%nchildren
        if (this%convex_polyhedron(t)%is_inside(x)) then
          is_inside = .true.
          exit
        end if
      end do
    else
      is_inside = this%parent%is_inside(x)
    end if

  end function is_inside

  subroutine print_data (this,normalized)
    class(polyhedron), intent(in) :: this
    logical, optional, intent(in) :: normalized
    call this%parent%print_data(normalized)
  end subroutine print_data

end module polyhedron_type
