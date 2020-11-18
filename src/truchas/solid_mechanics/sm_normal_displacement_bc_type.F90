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

module sm_normal_displacement_bc_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use bndry_func1_class
  use bndry_face_func_type
  use integration_geometry_type
  implicit none
  private

  type, extends(bndry_func1), public :: sm_normal_displacement_bc
    private
    real(r8), public, allocatable :: rotation_matrix(:,:,:)

    type(bndry_face_func), allocatable :: bff
    real(r8), allocatable :: area_ip(:)
    integer, allocatable :: fini(:), xfini(:)
  contains
    procedure :: init
    procedure :: compute
  end type sm_normal_displacement_bc

contains

  subroutine init(this, mesh, ig, bff)

    use sm_bc_utilities

    class(sm_normal_displacement_bc), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(integration_geometry), intent(in) :: ig
    type(bndry_face_func), intent(inout), allocatable :: bff

    integer :: xni, n, ni, fi
    real(r8), allocatable :: normal_ip(:,:), normal_node(:,:)

    ASSERT(allocated(bff))

    call move_alloc(bff, this%bff)
    call compute_index_connectivity(mesh, this%bff%index, this%fini, this%xfini, this%index)

    n = this%xfini(size(this%xfini))-1
    allocate(normal_ip(3,n), this%area_ip(n))
    allocate(this%rotation_matrix(3,3,size(this%index)), this%value(size(this%index)))
    allocate(normal_node(3,size(this%index)))
    call compute_ip_normals(this%bff%index, this%xfini, mesh, ig, normal_ip)
    call compute_node_normals(this%xfini, this%fini, normal_ip, normal_node)

    do n = 1, size(normal_ip, dim=2)
      this%area_ip(n) = norm2(normal_ip(:,n))
    end do

    do ni = 1, size(this%index)
      this%rotation_matrix(:,:,ni) = rotation_matrix(normal_node(:,ni))
    end do

  end subroutine init


  !! The value function corresponds to a displacement, for both the bff and this
  !! derived type. We map from face-centers to integration points via direct
  !! copy, then integration points to nodes via an area-weighted average.
  subroutine compute(this, t)

    class(sm_normal_displacement_bc), intent(inout) :: this
    real(r8), intent(in) :: t

    integer :: fi, xni, ni
    real(r8) :: v, weight(size(this%value))
    
    call this%bff%compute(t)
    this%value = 0
    weight = 0
    do fi = 1, size(this%bff%value)
      v = this%bff%value(fi)
      do xni = this%xfini(fi), this%xfini(fi+1)-1
        ni = this%fini(xni)
        this%value(ni) = this%value(ni) + this%area_ip(ni) * v
        weight(ni) = weight(ni) + this%area_ip(ni)
      end do
    end do

    do ni = 1, size(this%value)
      this%value(ni) = this%value(ni) / weight(ni)
    end do

  end subroutine compute

end module sm_normal_displacement_bc_type
