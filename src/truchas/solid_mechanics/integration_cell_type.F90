!!
!! TODO -- documentation
!!
!! By convention, the IP normal vector is oriented in the same direction as the
!! edge as defined in cell_topology. I.e., from *_edge(1,e) to *_edge(2,e). This
!! is enforced by the ordering of the nodes in the *_ipface arrays.
!!
!! Zach Jibben <zjibben@lanl.gov>
!! August 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module integration_cell_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use cell_topology
  implicit none
  private

  type, public :: integration_cell
    private
    integer :: nnode, nface, nedge
    real(r8), allocatable :: x(:,:)
    integer, pointer :: subcell(:) => null()
    integer, pointer :: xsubcell(:) => null()
    integer, pointer :: ipface(:,:) => null()
    integer, pointer :: faces(:) => null()
    integer, pointer :: xface(:) => null()
    integer, pointer, public :: fsize(:) => null()
    integer, pointer :: edges(:,:) => null()
  contains
    procedure :: init
    procedure :: subvolume
    procedure :: normal
    procedure :: normal_boundary
  end type integration_cell

  ! TODO-WARN
  integer, target :: tet4_subcell(40), tet4_xsubcell(5), tet4_ipface(4,6)
  ! data TET4_SUBCELL //
  ! data TET4_XSUBCELL //
  ! data TET4_IPFACE //

  integer, target :: pyr5_subcell(0), pyr5_xsubcell(6), pyr5_ipface(4,8)
  ! data PYR5_SUBCELL //
  ! data PYR5_XSUBCELL //
  ! data PYR5_IPFACE //

  integer, target :: wed6_subcell(0), wed6_xsubcell(7), wed6_ipface(4,9)
  ! data WED6_SUBCELL //
  ! data WED6_XSUBCELL //
  ! data WED6_IPFACE //

  integer, target :: hex8_subcell(64), hex8_xsubcell(9), hex8_ipface(4,12)
  data HEX8_SUBCELL /1,15,13,16,17,9,27,12, &
      &              2,18,13,15,19,10,27,9, &
      &              3,20,13,18,21,11,27,10, &
      &              4,16,13,20,22,12,27,11, &
      &              5,24,14,23,17,12,27,9, &
      &              6,23,14,25,19,9,27,10, &
      &              7,25,14,26,21,10,27,11, &
      &              8,26,14,24,22,11,27,12/
  data HEX8_XSUBCELL /1,9,17,25,33,41,49,57,65/
  data HEX8_IPFACE /15,13,27,9, &
      &             16,12,27,13, &
      &             17,9,27,12, &
      &             18,13,27,10, &
      &             19,10,27,9, &
      &             20,13,27,11, &
      &             21,11,27,10, &
      &             22,12,27,11, &
      &             23,9,27,14, &
      &             24,14,27,12, &
      &             25,10,27,14, &
      &             26,11,27,14/

contains

  subroutine init(this, x)

    class(integration_cell), intent(out) :: this
    real(r8), intent(in), target :: x(:,:)

    integer :: j, f, e

    this%nnode = size(x, dim=2)

    select case (this%nnode)
    case (4)
      this%subcell => tet4_subcell
      this%xsubcell => tet4_xsubcell
      this%ipface => tet4_ipface
      this%faces => tet4_faces
      this%xface => tet4_xface
      this%fsize => tet4_fsize
      this%edges => tet4_edges
      INSIST(.false.) ! TODO-WARN: not implemented yet
    case (5)
      this%subcell => pyr5_subcell
      this%xsubcell => pyr5_xsubcell
      this%ipface => pyr5_ipface
      this%faces => pyr5_faces
      this%xface => pyr5_xface
      this%fsize => pyr5_fsize
      this%edges => pyr5_edges
      INSIST(.false.) ! TODO-WARN: not implemented yet
    case (6)
      this%subcell => wed6_subcell
      this%xsubcell => wed6_xsubcell
      this%ipface => wed6_ipface
      this%faces => wed6_faces
      this%xface => wed6_xface
      this%fsize => wed6_fsize
      this%edges => wed6_edges
      INSIST(.false.) ! TODO-WARN: not implemented yet
    case (8)
      this%subcell => hex8_subcell
      this%xsubcell => hex8_xsubcell
      this%ipface => hex8_ipface
      this%faces => hex8_faces
      this%xface => hex8_xface
      this%fsize => hex8_fsize
      this%edges => hex8_edges
    case default
      ASSERT(.false.)
    end select
    this%nface = size(this%xface)-1
    this%nedge = size(this%edges, dim=2)

    allocate(this%x(3,this%nnode + this%nface + this%nedge + 1))

    this%x(:,:this%nnode) = x

    do f = 1, this%nface
      j = this%nnode + f
      this%x(:,j) = sum(this%x(:,this%faces(this%xface(f):this%xface(f+1)-1)), dim=2) / this%fsize(f)
    end do

    do e = 1, this%nedge
      j = this%nnode + this%nface + e
      this%x(:,j) = sum(this%x(:,this%edges(:,e)), dim=2) / 2
    end do

    j = this%nnode+this%nface+this%nedge+1
    this%x(:,j) = sum(this%x(:,:this%nnode), dim=2) / this%nnode

  end subroutine init


  !! Returns the area-weighted normal for the surface associated with
  !! integration point p. By convention, the normal vector is oriented in the
  !! same direction as the edge associated with p.
  function normal(this, p)
    use cell_geometry, only: face_normal
    class(integration_cell), intent(in) :: this
    integer, intent(in) :: p
    real(r8) :: normal(3)
    normal = face_normal(this%x(:,this%ipface(:,p)))
  end function normal


  !! Returns the area-weighted normal for the boundary surface associated with
  !! face f and node k local to that face. The normal vector is oriented outward
  !! with respect to the cell/domain. This integration point surface is
  !! defined by the loop from the referenced cell, to the next edge center (when
  !! looping through the face in a right-handed fashion), to the face center, to
  !! the previous edge center.
  function normal_boundary(this, f, k) result(normal)

    use cell_geometry, only: face_normal

    class(integration_cell), intent(in) :: this
    integer, intent(in) :: f, k
    real(r8) :: normal(3)

    real(r8) :: x(3,4)
    integer :: knext, kprev, j, jnext, jprev

    knext = modulo(k, this%fsize(f))+1
    kprev = modulo(k-2, this%fsize(f))+1
    j = this%faces(this%xface(f)+k-1)
    jnext = this%faces(this%xface(f)+knext-1)
    jprev = this%faces(this%xface(f)+kprev-1)

    x(:,1) = this%x(:,j)
    x(:,2) = (this%x(:,j) + this%x(:,jnext)) / 2
    x(:,3) = this%x(:,this%nnode+f)
    x(:,4) = (this%x(:,j) + this%x(:,jprev)) / 2

    normal = face_normal(x)

  end function normal_boundary


  !! Compute the subvolume of the region associated with node n of this cell.
  !! This subvolume is bounded by the adjacent integration points.
  real(r8) function subvolume(this, n)
    use cell_geometry, only: cell_volume
    class(integration_cell), intent(in) :: this
    integer, intent(in) :: n
    subvolume = cell_volume(this%x(:,this%subcell(this%xsubcell(n):this%xsubcell(n+1)-1)))
  end function subvolume

end module integration_cell_type
