!!
!! This module provides a type for handling boundary conditions on the
!! integration points of the solid mechanics discretization. There are one of
!! these points for each node of each face. While the time-dependent data will
!! only be associated with the face, there are geometric factors associated with
!! each integration point that necessitate this independent type.
!!
!! Zach Jibben <zjibben@lanl.gov>
!! October 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module bndry_ip_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use bndry_func1_class
  use bndry_face_func_type
  implicit none
  private

  type, extends(bndry_func1), public :: bndry_ip_func
    private
    real(r8), public, allocatable :: factor(:)

    type(bndry_face_func), allocatable :: bff
    class(unstr_mesh), pointer :: mesh => null() ! reference only - do not own
  contains
    procedure :: init
    procedure :: compute
  end type bndry_ip_func

contains

  subroutine init(this, mesh, ig, bff)

    use integration_geometry_type

    class(bndry_ip_func), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(integration_geometry), intent(in) :: ig
    type(bndry_face_func), intent(inout), allocatable :: bff

    integer :: i, f, xf, xn, n, j
    real(r8) :: area(4)

    this%mesh => mesh
    call move_alloc(bff, this%bff)

    associate (faces => this%bff%index)
      !! Get the size of index (number of on-rank boundary integration points in
      !! this face set)
      j = 0
      do i = 1, size(faces)
        f = faces(i)
        do xn = this%mesh%xfnode(f), this%mesh%xfnode(f+1)-1
          n = this%mesh%fnode(xn)
          if (n <= this%mesh%nnode_onP) j = j+1
        end do
      end do

      !! Populate the index and constant geometric value factors
      !! NB: This is slightly wasteful, as these geometric factors are duplicated
      !!     for every boundary condition associated with a face set. They are
      !!     only computed once, but they will consume some memory.
      allocate(this%index(j), this%value(j), this%factor(j))
      j = 1
      do i = 1, size(faces)
        f = faces(i)
        call area_boundary_ips(this, ig, f, area)
        associate (xn => this%mesh%fnode(this%mesh%xfnode(f):this%mesh%xfnode(f+1)-1))
          do xf = 1, size(xn)
            n = xn(xf)
            if (n > this%mesh%nnode_onP) cycle
            this%index(j) = n
            this%factor(j) = area(xf)
            j = j+1
          end do
        end associate
      end do
    end associate

  end subroutine init


  !! Return the areas for each integration surface on a face. Each surface
  !! is associated with a node.
  subroutine area_boundary_ips(this, ig, f, area)

    use integration_cell_type
    use integration_geometry_type

    type(bndry_ip_func), intent(in) :: this
    type(integration_geometry), intent(in) :: ig
    integer, intent(in) :: f
    real(r8), intent(out) :: area(:)

    integer :: j, xf
    integer, pointer :: cn(:) => null()
    type(integration_cell) :: ic

    xf = ig%fface(1,f) ! local face index
    j = this%mesh%fcell(1,f) ! cell index
    cn => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1)
    call ic%init(this%mesh%x(:,cn))

    !area = [(norm2(ic%normal_boundary(xf,j)), j = 1, ic%fsize(xf))] ! NAG incorrect evaluation
    ASSERT(size(area) >= ic%fsize(xf))
    do j = 1, ic%fsize(xf)
      area(j) = norm2(ic%normal_boundary(xf,j))
    end do

  end subroutine area_boundary_ips


  !! Computes the underlying boundary face function, and maps the result onto
  !! the integration points, applying the constant geometric factor along the
  !! way.
  subroutine compute(this, t)

    class(bndry_ip_func), intent(inout) :: this
    real(r8), intent(in) :: t

    integer :: i, j, f, xf, n
    real(r8) :: v

    call this%bff%compute(t)

    j = 1
    do i = 1, size(this%bff%value)
      f = this%bff%index(i)
      v = this%bff%value(i)
      associate (xn => this%mesh%fnode(this%mesh%xfnode(f):this%mesh%xfnode(f+1)-1))
        do xf = 1, size(xn)
          n = xn(xf)
          if (n > this%mesh%nnode_onP) cycle
          this%value(j) = v
          j = j + 1
        end do
      end associate
    end do

  end subroutine compute

end module bndry_ip_func_type
