!!
!! MESH_FACE_SET_IMPL
!!
!! Implements the legacy MESH_FACE_SET array.
!!

#include "f90_assert.fpp"

module mesh_face_set_impl

  implicit none
  private

  public :: init_mesh_face_set_impl

  integer, allocatable, public, protected :: mesh_face_set(:,:,:)

contains

  subroutine init_mesh_face_set_impl

    use common_impl, only: ncells, mesh => new_mesh, pcell_old_to_new
    use common_impl, only: NEW_TET_SIDE_MAP, NEW_PYR_SIDE_MAP, NEW_PRI_SIDE_MAP, NEW_HEX_SIDE_MAP
    use parallel_permutations, only: rearrange
    use bitfield_type

    integer :: i, j, k, kk
    type(bitfield) :: ssmask
    integer, allocatable :: src(:,:,:), buffer(:)

    ssmask = not(ibset(ZERO_BITFIELD,pos=0))

    allocate(src(size(mesh%face_set_id),6,mesh%ncell_onP))
    allocate(buffer(size(src,1)))

    src = 0
    do j = 1, mesh%ncell_onP
      associate (cface => mesh%cface(mesh%xcface(j):mesh%xcface(j+1)-1))
        do k = 1, size(cface)
          if (iand(ssmask,mesh%face_set_mask(cface(k))) == ZERO_BITFIELD) cycle
          buffer = 0
          do i = 1, size(mesh%face_set_id)
            if (btest(mesh%face_set_mask(cface(k)),pos=i)) buffer(i) = mesh%face_set_id(i)
          end do
          select case (mesh%xcnode(j+1)-mesh%xcnode(j))
          case (4)
            kk = NEW_TET_SIDE_MAP(k)
          case (5)
            kk = NEW_PYR_SIDE_MAP(k)
          case (6)
            kk = NEW_PRI_SIDE_MAP(k)
          case (8)
            kk = NEW_HEX_SIDE_MAP(k)
          case default
            INSIST(.false.)
          end select
          src(:,kk,j) = buffer
        end do
      end associate
    end do

    allocate(mesh_face_set(size(src,1),size(src,2),ncells))
    do i = 1, size(src,1)
      call rearrange (pcell_old_to_new, mesh_face_set(i,:,:), src(i,:,:), default=0)
    end do

  end subroutine init_mesh_face_set_impl

end module mesh_face_set_impl
