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
  implicit none
  private

  type, public :: sm_bc_face_list
    public
    integer :: nbc
    real(r8), allocatable :: normal(:,:)
    integer, allocatable :: face(:) ! stores the face id
    integer, allocatable :: fbc(:), xfbc(:) ! stores bc values for each face
  contains
    procedure :: init
  end type sm_bc_face_list

contains

  subroutine init(this, mesh, bc_list, stat, errmsg)

    use unstr_mesh_type
    use sm_bc_list_type

    class(sm_bc_face_list), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(sm_bc_list), intent(in) :: bc_list
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: j, n, k, nface, bcid
    integer, allocatable :: nbc(:)

    stat = 0

    ! 1. Count the BCs on each face
    allocate(nbc(mesh%nface))
    nbc = 0
    do n = 1, size(bc_list%displacement)
      k = findloc(mesh%face_set_id, bc_list%displacement(n)%setid)
      do j = 1, mesh%nface
        if (btest(mesh%face_set_mask(j),k)) nbc(j) = nbc(j) + 1
      end do
    end do
    do n = 1, size(bc_list%contact)
      k = findloc(mesh%face_set_id, bc_list%contact(n)%setid)
      do j = 1, mesh%nface
        if (btest(mesh%face_set_mask(j),k)) nbc(j) = nbc(j) + 1
      end do
    end do

    ! 2. Make a list of faces with BCs
    nface = count(nbc > 0)
    allocate(this%xfbc(nface+1), this%face(nface))
    nface = 1
    this%xfbc(1) = 1
    do j = 1, mesh%nface
      if (nbc(j) > 0) then
        this%index(nface) = j
        this%xfbc(nface+1) = this%xfbc(nface) + nbc(j)
        nface = nface + 1
      end if
    end do

    ! 3. List the BCs on faces having any (the list made in step 2)
    deallocate(nbc)
    allocate(nbc(size(this%face)), this%fbc(this%xfbc(size(this%xfbc))))
    nbc = 0
    do n = 1, size(bc_list%displacement)
      bcid = n
      k = findloc(mesh%face_set_id, bc_list%displacement(n)%setid)
      do i = 1, size(this%face)
        j = this%face(i)
        if (btest(mesh%face_set_mask(j),k)) then
          this%fbc(this%xfbc(i)+nbc(i)) = bcid
          nbc(i) = nbc(i) + 1
        end if
      end do
    end do
    do n = 1, size(bc_list%contact)
      bcid = size(bc_list%displacement) + n
      k = findloc(mesh%face_set_id, bc_list%contact(n)%setid)
      do i = 1, size(this%face)
        j = this%face(i)
        if (btest(mesh%face_set_mask(j),k)) then
          this%fbc(this%xfbc(i)+nbc(i)) = bcid
          nbc(i) = nbc(i) + 1
        end if
      end do
    end do

  end subroutine init

end module sm_bc_face_list_type
