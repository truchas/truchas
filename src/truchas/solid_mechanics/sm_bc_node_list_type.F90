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
    integer, allocatable :: bcid(:), xbcid(:) ! stores bc identifiers for each node

    ! Stores the ID of the node opposite a given contact BC. Indexed
    ! via xbcid. non-contact BCs have a linked_node value of -1.
    integer, allocatable :: linked_node(:)

    integer, allocatable :: link(:,:) ! stores indexes to the node array above at links
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
  subroutine init(this, mesh, ig, bc_list, bc_face_list)

    use parallel_communication, only: global_sum
    use unstr_mesh_type
    use integration_geometry_type
    use sm_bc_list_type
    use sm_bc_face_list_type
    use sm_bc_utilities, only: compute_index_connectivity, compute_ip_normals

    class(sm_bc_node_list), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(integration_geometry), intent(in) :: ig
    type(sm_bc_list), intent(in) :: bc_list
    type(sm_bc_face_list), intent(in) :: bc_face_list

    character(32) :: msg
    integer :: nnode, fi, ni, xni, n, bcid, xbcid
    integer, allocatable :: fini(:), xfini(:), link_setid(:)
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
    allocate(this%bcid(n), this%normal(3,n), this%linked_node(n))
    do ni = 1, size(this%node)
      xbcid = this%xbcid(ni)
      do bcid = 1, bc_face_list%nbc
        if (bc_is_active(bcid,ni)) then
          this%bcid(xbcid) = bcid
          this%normal(:,xbcid) = select_normal(bc_list%bc_type(bcid), normal(:,bcid,ni))
          xbcid = xbcid + 1
        end if
      end do
      this%xbcid(ni+1) = this%xbcid(ni) + count(bc_is_active(:,ni))
    end do

    call compute_links
    call compute_linked_nodes

    call append_nodeset_bcs

    nnode = count(this%node <= mesh%nnode_onP)
    n = count(this%link(1,:) <= mesh%nnode_onP .or. this%link(2,:) <= mesh%nnode_onP)
    nnode = global_sum(nnode)
    n = global_sum(n)
    write(msg,"('SM BC nodes/links: ',2i6)") nnode, n
    call TLS_info(trim(msg))

  contains

    ! Use definitions from sm_bc_list_type.F90 to identify the direction of this BC.
    function select_normal(type, face_normal) result(normal)
      integer, intent(in) :: type
      real(r8), intent(in) :: face_normal(:)
      real(r8) :: normal(3)
      select case (type)
      case (SMBCL_N)
        normal = face_normal
      case (SMBCL_X)
        normal = [1.0_r8, 0.0_r8, 0.0_r8]
      case (SMBCL_Y)
        normal = [0.0_r8, 1.0_r8, 0.0_r8]
      case (SMBCL_Z)
        normal = [0.0_r8, 0.0_r8, 1.0_r8]
      case default
        INSIST(.false.)
      end select
    end function select_normal

    ! Generate a list of linked nodes for gaps. The mesh structure holds a list
    ! of linked faces, which the integration geometry converts to a list of
    ! linked nodes. Here we generate a list of linked node *indexes* to
    ! this%node(:). That is, ig%lnode(:,n) = this%node(this%link(:,n)), for
    ! every n in the length of this%link. Note there might be links in the mesh
    ! which are not assigned contact conditions.
    subroutine compute_links()

      integer :: l, nlink, node_index(mesh%nnode)

      node_index = 0
      do ni = 1, size(this%node)
        n = this%node(ni)
        node_index(n) = ni
      end do

      nlink = 0
      do l = 1, size(ig%lnode,dim=2)
        if (any(node_index(ig%lnode(:,l)) > 0)) then
          ASSERT(all(node_index(ig%lnode(:,l)) > 0)) ! we need off-rank nodes to be accounted for
          nlink = nlink + 1
        end if
      end do

      allocate(this%link(2,nlink), link_setid(nlink))
      nlink = 0
      do l = 1, size(ig%lnode,dim=2)
        if (all(node_index(ig%lnode(:,l)) > 0)) then
          nlink = nlink + 1
          this%link(:,nlink) = node_index(ig%lnode(:,l))
          link_setid(nlink) = ig%link_setid(l)
        end if
      end do


    end subroutine compute_links


    ! Generate a list of linked nodes corresponding to contact BCs
    subroutine compute_linked_nodes()

      integer :: li, lni, ni, opposite_ni, xbcid, bcid

      this%linked_node = -1

      do li = 1, size(this%link, dim=2)
        do lni = 1, 2
          ni = this%link(lni,li)
          opposite_ni = this%link(modulo(lni,2)+1,li)
          do xbcid = this%xbcid(ni), this%xbcid(ni+1)-1
            bcid = this%bcid(xbcid)
            if (bcid >= bc_list%xcontact) then
              if (bc_list%contact(bcid-bc_list%xcontact+1)%setid == link_setid(li)) &
                  this%linked_node(xbcid) = this%node(opposite_ni)
            end if
          end do
        end do
      end do

    end subroutine compute_linked_nodes


    subroutine append_nodeset_bcs()

      logical :: applied
      integer :: nbc, nbcnode, n, ki, k, nnode
      integer, allocatable :: nodeset_index(:), nodeset_bcid(:), xbcid(:), bcid(:), node(:)
      real(r8), allocatable :: normal(:,:)

      nnode = size(this%node)

      ! get a list of the BCIDs for nodeset BCs, and a list of the associated nodesets.
      nodeset_index = pack(bc_list%displacement(:)%setid, mask=bc_list%displacement(:)%nodeset)
      do ki = 1, size(nodeset_index)
        nodeset_index(ki) = findloc(mesh%node_set_id, nodeset_index(ki), dim=1)
      end do
      allocate(nodeset_bcid(size(nodeset_index)))
      k = 0
      do ki = 1, bc_list%xcontact-1
        if (bc_list%displacement(ki)%nodeset) then
          k = k + 1
          nodeset_bcid(k) = ki
        end if
      end do

      ! count the number of newly applied boundary conditions across the mesh
      nbc = 0
      ni = 0
      do n = 1, mesh%nnode_onP
        applied = .false.
        do ki = 1, size(nodeset_index)
          k = nodeset_index(ki)
          if (btest(mesh%node_set_mask(n),k)) then
            nbc = nbc + 1
            applied = .true.
          end if
        end do
        if (applied) ni = ni + 1
      end do

      ! list the new BCs in a ragged array
      allocate(bcid(nbc), xbcid(ni+1), node(ni), normal(3,nbc))
      xbcid(1) = 1
      ni = 1
      do n = 1, mesh%nnode_onP
        applied = .false.
        nbc = 0
        do ki = 1, size(nodeset_index)
          k = nodeset_index(ki)
          if (btest(mesh%node_set_mask(n),k)) then
            bcid(xbcid(ni)+nbc) = nodeset_bcid(ki)
            normal(:,xbcid(ni)+nbc) = select_normal(bc_list%bc_type(nodeset_bcid(ki)), [0d0,0d0,0d0])
            nbc = nbc + 1
            applied = .true.
          end if
        end do
        if (applied) then
          node(ni) = n
          xbcid(ni+1) = xbcid(ni) + nbc
          ni = ni + 1
        end if
      end do

      call concatenate_ragged_arrays(this%node, this%xbcid, this%bcid, this%normal, &
          node, xbcid, bcid, normal)

    end subroutine append_nodeset_bcs

  end subroutine init


  subroutine concatenate_ragged_arrays(id1, xarr1, arr1, rarr1, id2, xarr2, arr2, rarr2)

    integer, intent(inout), allocatable :: id1(:), xarr1(:), arr1(:)
    real(r8), intent(inout), allocatable :: rarr1(:,:)
    integer, intent(in) :: id2(:), xarr2(:), arr2(:)
    real(r8), intent(in) :: rarr2(:,:)

    integer :: ni, ni2, n, k, nnode
    integer, allocatable :: id3(:), xarr3(:), arr3(:), ni_to_ni2(:)
    real(r8), allocatable :: rarr3(:,:)
    logical, allocatable :: counted(:)

    allocate(ni_to_ni2(size(id1)), counted(size(id2)))
    counted = .false.
    ni_to_ni2 = 0
    nnode = size(id1) + size(id2)
    do ni = 1, size(id1)
      ni2 = findloc(id2, id1(ni), dim=1)
      if (ni2 > 0 .and. ni2 <= size(id2)) then
        ni_to_ni2(ni) = ni2
        nnode = nnode - 1
        counted(ni2) = .true.
      end if
    end do

    allocate(id3(nnode), xarr3(nnode+1), arr3(size(arr1) + size(arr2)), &
        rarr3(3,size(arr1) + size(arr2)))
    xarr3(1) = 1
    do ni = 1, size(id1)
      id3(ni) = id1(ni)
      xarr3(ni+1) = xarr3(ni) + (xarr1(ni+1) - xarr1(ni))
      arr3(xarr3(ni):xarr3(ni+1)-1) = arr1(xarr1(ni):xarr1(ni+1)-1)
      rarr3(:,xarr3(ni):xarr3(ni+1)-1) = rarr1(:,xarr1(ni):xarr1(ni+1)-1)

      ni2 = ni_to_ni2(ni)
      if (ni2 > 0) then
        k = xarr3(ni+1) + (xarr2(ni2+1) - xarr2(ni2))
        arr3(xarr3(ni+1):k) = arr2(xarr2(ni2):xarr2(ni2+1)-1)
        rarr3(:,xarr3(ni+1):k) = rarr2(:,xarr2(ni2):xarr2(ni2+1)-1)
        xarr3(ni+1) = k
      end if
    end do

    do ni2 = 1, size(id2)
      if (counted(ni2)) cycle
      ni = size(id1) + ni2
      id3(ni) = id2(ni2)
      xarr3(ni+1) = xarr3(ni) + (xarr2(ni2+1) - xarr2(ni2))
      arr3(xarr3(ni):xarr3(ni+1)-1) = arr2(xarr2(ni2):xarr2(ni2+1)-1)
      rarr3(:,xarr3(ni):xarr3(ni+1)-1) = rarr2(:,xarr2(ni2):xarr2(ni2+1)-1)
    end do

    id1 = id3
    xarr1 = xarr3
    arr1 = arr3
    rarr1 = rarr3

  end subroutine concatenate_ragged_arrays

end module sm_bc_node_list_type
