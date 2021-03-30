!!
!! SM_BC_NODE_LIST_TYPE
!!
!! TODO
!!
!! Zach Jibben <zjibben@lanl.gov>
!! March 2021
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module sm_bc_node_list_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use truchas_logging_services
  implicit none
  private

  type, public :: sm_bc_node_list
    real(r8), allocatable :: normal(:,:)
    integer, allocatable :: node(:) ! stores the node id
    integer, allocatable :: bcid(:), xbcid(:) ! stores bc values for each node
  contains
    procedure :: init
  end type sm_bc_node_list

contains

  !! Given is a list of faces, each with a list of unique BCs associated with
  !! that face. From that we produce a list of nodes (attached to those faces),
  !! each with a list of unique BCs associated with that node. The main
  !! challenge is in not double-counting BCs, as a node can neighbor multiple
  !! faces with the same BC.
  !!
  !! Also a challenge is the accumulation of area-weighted normals. Normals are
  !! face-normals scaled by the area of the integration point associated with
  !! the node-face pair. And a normal associated with a particular BC at a
  !! particular node is the sum of the neighboring area-scaled normals at
  !! integration points associated with that BC.
  subroutine init(this, mesh, ig, bc_face_list)

    use unstr_mesh_type
    use integration_geometry_type
    use sm_bc_face_list_type
    use sm_bc_utilities, only: compute_index_connectivity, compute_ip_normals

    class(sm_bc_node_list), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(integration_geometry), intent(in) :: ig
    type(sm_bc_face_list), intent(in) :: bc_face_list

    character(32) :: msg
    integer :: nnode, fi, ni, xni, n, bcid, xbcid
    integer, allocatable :: fini(:), xfini(:)
    logical, allocatable :: bc_is_active(:,:)
    real(r8), allocatable :: normal(:,:,:), normal_ip(:,:)

    call compute_index_connectivity(mesh, bc_face_list%face, fini, xfini, this%node)
    allocate(normal_ip(3,xfini(size(xfini))-1))
    call compute_ip_normals(bc_face_list%face, xfini, mesh, ig, normal_ip)

    nnode = size(this%node)
    allocate(normal(3,bc_face_list%nbc,nnode), bc_is_active(bc_face_list%nbc,nnode))
    normal = 0
    bc_is_active = .false.
    do fi = 1, size(bc_face_list%face)
      do xni = xfini(fi), xfini(fi+1)-1
        ni = fini(xni)
        n = this%node(ni)

        do xbcid = bc_face_list%xbcid(fi), bc_face_list%xbcid(fi+1)-1
          bcid = bc_face_list%bcid(xbcid)

          normal(:,bcid,ni) = normal(:,bcid,ni) + normal_ip(:,xni)
          bc_is_active(bcid,ni) = .true.
        end do
      end do
    end do

    allocate(this%xbcid(nnode+1))
    this%xbcid(1) = 1
    n = count(bc_is_active)
    allocate(this%bcid(n), this%normal(3,n))
    do ni = 1, size(this%node)
      xbcid = this%xbcid(ni)
      do bcid = 1, bc_face_list%nbc
        if (bc_is_active(bcid,ni)) then
          this%bcid(xbcid) = bcid
          this%normal(:,xbcid) = normal(:,bcid,ni)
          xbcid = xbcid + 1
        end if
      end do
      this%xbcid(ni+1) = this%xbcid(ni) + count(bc_is_active(:,ni))
    end do

    nnode = size(this%xbcid)-1
    write(msg,"('SM BC nodes: ',i6)") nnode
    call TLS_info(trim(msg))

    ! ! ! 1. Get a map from nodes to faces
    ! ! call compute_inverted_connectivity(mesh%fnode, mesh%xfnode, mesh%nnode, nface, xnface)

    ! do fi = 1, size(bc_face_list%face)
    !   f = bc_face_list%face(fi)
    !   do xni = mesh%xfnode(f), mesh%xfnode(f+1)-1
    !     n = mesh%fnode(xni)
    !     nbc(n) = nbc(n) + 1
    !   end do
    ! end do

  end subroutine init

end module sm_bc_node_list_type
