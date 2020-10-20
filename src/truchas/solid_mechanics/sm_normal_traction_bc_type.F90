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

module sm_normal_traction_bc_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use bndry_vfunc_class
  use bndry_face_func_type
  use integration_geometry_type
  implicit none
  private

  type, extends(bndry_vfunc), public :: sm_normal_traction_bc
    private
    type(bndry_face_func), allocatable :: bff
    integer, allocatable :: fini(:), xfini(:)
    real(r8), allocatable :: normal_ip(:,:)
  contains
    procedure :: init
    procedure :: compute
  end type sm_normal_traction_bc

contains

  subroutine init(this, mesh, ig, bff)

    use sm_bc_utilities

    class(sm_normal_traction_bc), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(integration_geometry), intent(in) :: ig
    type(bndry_face_func), intent(inout), allocatable :: bff

    ASSERT(allocated(bff))

    call move_alloc(bff, this%bff)
    call compute_index_connectivity(mesh, this%bff%index, this%fini, this%xfini, this%index)
    allocate(this%normal_ip(3,this%xfini(size(this%xfini))-1), this%value(3,size(this%index)))
    call compute_ip_normals(this%bff%index, this%xfini, mesh, ig, this%normal_ip)

  end subroutine init


  !! Computes the underlying boundary face function, and maps the result onto
  !! the integration points, applying the constant geometric factor along the
  !! way.
  !!
  !! NB: Note this mapping from face to integration points is effectively
  !! first-order, so anything which is not a constant normal traction will incur
  !! some error. To fix this, we need an alternative to bndry_face_func_type for
  !! arbitrary integration points, so user-provided functions are evalutated at
  !! each IP. I imagine non-constant normal traction boundary conditions are not
  !! common. -zjibben
  subroutine compute(this, t)

    class(sm_normal_traction_bc), intent(inout) :: this
    real(r8), intent(in) :: t

    integer :: fi, xni, ni
    real(r8) :: v

    call this%bff%compute(t)
    this%value = 0
    do fi = 1, size(this%bff%value)
      v = this%bff%value(fi)
      do xni = this%xfini(fi), this%xfini(fi+1)-1
        ni = this%fini(xni)
        this%value(:,ni) = this%value(:,ni) + v * this%normal_ip(:,xni)
      end do
    end do

  end subroutine compute

end module sm_normal_traction_bc_type
