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
  implicit none
  private

  type, public :: sm_bc_face_list
    public
    real(r8), allocatable :: normal(:,:)
    integer, allocatable :: node(:) ! stores the node id
    integer, allocatable :: bcid(:), offset(:) ! stores bc values for each node
  contains
    procedure :: init
  end type sm_bc_face_list

contains

  subroutine init(this, mesh, bc_face_list, stat, errmsg)

    use unstr_mesh_type
    use sm_bc_face_list_type
    use sm_bc_utilities, only: compute_inverted_connectivity

    class(sm_bc_node_list), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(sm_bc_face_list), intent(in) :: bc_face_list
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer, allocatable :: nface(:), xnface(:)

    stat = 0

    ! ! 1. Get a map from nodes to faces
    ! call compute_inverted_connectivity(mesh%fnode, mesh%xfnode, mesh%nnode, nface, xnface)

    do fi = 1, size(bc_face_list%face)
      f = bc_face_list%face(fi)
      do xni = mesh%xfnode(f), mesh%xfnode(f+1)-1
        n = mesh%fnode(xni)
        nbc(n) = nbc(n) + 1
      end do
    end do

  end subroutine init

end module sm_bc_node_list_type
