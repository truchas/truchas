!!
!! SM_BC_FACE_LIST_TYPE
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

module sm_bc_face_list_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use truchas_logging_services
  implicit none
  private

  type, public :: sm_bc_face_list
    integer :: nbc
    integer, allocatable :: face(:) ! stores the face id
    integer, allocatable :: bcid(:), xbcid(:) ! stores bc values for each face
  contains
    procedure :: init
  end type sm_bc_face_list

contains

  subroutine init(this, mesh, bc_list)

    use parallel_communication, only: global_sum
    use unstr_mesh_type
    use sm_bc_list_type
    use bitfield_type

    class(sm_bc_face_list), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(sm_bc_list), intent(in) :: bc_list

    character(32) :: msg
    integer :: i, j, k, n, nface, bcid
    integer, allocatable :: nbc(:)
    type(bitfield), allocatable :: face_set_mask(:)

    this%nbc = count(.not.bc_list%displacement(:)%nodeset) + size(bc_list%contact)
    call get_regularized_face_set_mask(mesh, face_set_mask)

    ! 1. Count the BCs on each face
    allocate(nbc(mesh%nface))
    nbc = 0
    do n = 1, size(bc_list%displacement)
      if (bc_list%displacement(n)%nodeset) cycle ! nodeset displacement BCs added in the node list
      k = findloc(mesh%face_set_id, bc_list%displacement(n)%setid, dim=1)
      do j = 1, mesh%nface
        if (btest(mesh%face_set_mask(j),k)) nbc(j) = nbc(j) + 1
      end do
    end do
    do n = 1, size(bc_list%contact)
      k = findloc(mesh%face_set_id, bc_list%contact(n)%setid, dim=1)
      ASSERT(k > 0)
      do j = 1, mesh%nface
        ! use the local regularized face set mask to catch
        ! face set matches on both sides of links
        if (btest(face_set_mask(j),k)) nbc(j) = nbc(j) + 1
      end do
      ! k = findloc(mesh%link_set_id, bc_list%contact(n)%setid, dim=1)
      ! ASSERT(k > 0)
      ! do j = 1, mesh%nlink
      !   if (btest(mesh%link_set_mask(j),k)) then
      !     nbc(mesh%lface(1,j)) = nbc(mesh%lface(1,j)) + 1
      !     nbc(mesh%lface(2,j)) = nbc(mesh%lface(2,j)) + 1
      !   end if
      ! end do
    end do

    ! 2. Make a list of faces with BCs
    nface = count(nbc > 0)
    allocate(this%xbcid(nface+1), this%face(nface))
    nface = 1
    this%xbcid(1) = 1
    do j = 1, mesh%nface
      if (nbc(j) > 0) then
        this%face(nface) = j
        this%xbcid(nface+1) = this%xbcid(nface) + nbc(j)
        nface = nface + 1
      end if
    end do

    ! 3. List the BCs on faces having any (the list made in step 2)
    deallocate(nbc)
    allocate(nbc(size(this%face)), this%bcid(this%xbcid(size(this%xbcid))))
    nbc = 0
    do n = 1, size(bc_list%displacement)
      bcid = n
      if (bc_list%displacement(n)%nodeset) cycle ! nodeset displacement BCs added in the node list
      k = findloc(mesh%face_set_id, bc_list%displacement(n)%setid, dim=1)
      do i = 1, size(this%face)
        j = this%face(i)
        if (btest(mesh%face_set_mask(j),k)) then
          this%bcid(this%xbcid(i)+nbc(i)) = bcid
          nbc(i) = nbc(i) + 1
        end if
      end do
    end do
    do n = 1, size(bc_list%contact)
      bcid = size(bc_list%displacement) + n
      k = findloc(mesh%face_set_id, bc_list%contact(n)%setid, dim=1)
      do i = 1, size(this%face)
        j = this%face(i)
        ! use the local regularized face set mask to catch
        ! face set matches on both sides of links
        if (btest(face_set_mask(j),k)) then
          this%bcid(this%xbcid(i)+nbc(i)) = bcid
          nbc(i) = nbc(i) + 1
        end if
      end do
    end do

    nface = count(this%face <= mesh%nface_onP)
    nface = global_sum(nface)
    write(msg,"('SM BC faces: ',i6)") nface
    call TLS_info(trim(msg))

  end subroutine init


  !! The mesh component mesh%face_set_mask marks faces matching given face set
  !! IDs. Internal face set IDs can be used to generate gaps in the mesh, where
  !! the mesh is broken apart and a new set of faces introduced to pair with the
  !! existing faces. These pairs are listed in links, and mesh%link_set_mask
  !! lists those pairs of faces for the given face set. BUT, the new faces
  !! created aren't associated with the original face set. I.e., you can get to
  !! both sides of a link for a given face set ID through link_set_mask, but not
  !! through face_set_mask.
  !!
  !! In the above init routine, we need to iterate through a compressed list of
  !! faces and find the boundary conditions applied on each face. To do this, we
  !! need the new link-faces to be associated with the face set ID that their
  !! sibling faces are associated with.
  !!
  !! This routine generates such a "regularized" face set mask from the mesh
  !! structure to be used locally.
  subroutine get_regularized_face_set_mask(mesh, face_set_mask)

    use unstr_mesh_type
    use bitfield_type

    type(unstr_mesh), intent(in) :: mesh
    type(bitfield), intent(out), allocatable :: face_set_mask(:)

    integer :: kl, fid, kf, l
    integer :: face_set_index(size(mesh%link_set_id))

    face_set_mask = mesh%face_set_mask

    do kl = 1, size(mesh%link_set_id)
      fid = mesh%link_set_id(kl)
      kf = findloc(mesh%face_set_id, fid, dim=1)
      ASSERT(kf > 0)
      face_set_index(kl) = kf
    end do

    do l = 1, size(mesh%lface,dim=2)
      do kl = 1, size(face_set_index)
        if (btest(mesh%link_set_mask(l),kl)) then
          kf = face_set_index(kl)
          face_set_mask(mesh%lface(1,l)) = ibset(face_set_mask(mesh%lface(1,l)),kf)
          face_set_mask(mesh%lface(2,l)) = ibset(face_set_mask(mesh%lface(2,l)),kf)
        end if
      end do
    end do

  end subroutine get_regularized_face_set_mask

end module sm_bc_face_list_type
