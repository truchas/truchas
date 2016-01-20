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

    call test_mesh_face_set ! Caution: see remarks below

  end subroutine init_mesh_face_set_impl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! TESTING CODE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! There is a 1-1 mapping between exodus side sets and original mesh_face_sets,
  !! but some side set info is lost when they are converted into new mesh face
  !! sets.  Thus it is not generally possible to fully recover the side sets and
  !! compare to the original mesh_face_sets array.  However in practice it seems
  !! we can, and if there are differences, which will only occur for side sets
  !! internal to the mesh, that difference is likely not significant practically.
  !! The issue is that a (cell,side) pair maps to a unique face, but a face can
  !! map to two (cell,side) pairs (except on the boundary).  Exodus side sets
  !! may be one sided or two sided (or even some combination).  The above code
  !! takes the approach of generating 2-sided side sets from the face sets; the
  !! original may not be.

  subroutine test_mesh_face_set

    use mesh_module, only: old_mesh_face_set => mesh_face_set
    use parallel_communication, only: global_any
    use truchas_logging_services

    integer :: j, k, unit
    logical :: error

    unit = TLS_debug_unit()

    error = .false.
    do j = 1, size(mesh_face_set,dim=3)
      do k = 1, 6
        if (all(mesh_face_set(:,k,j) == old_mesh_face_set(:,k,j))) cycle
        if (.not.error) then
          write(unit,'(/,a)') 'MESH_FACE_SET **********************'
          error = .true.
        end if
        write(unit,'(2(a,i0),a,*(1x,i0))') 'old[', k, ',', j, ']=', old_mesh_face_set(:,k,j)
        write(unit,'(2(a,i0),a,*(1x,i0))') 'old[', k, ',', j, ']=', mesh_face_set(:,k,j)
      end do
    end do
    if (global_any(error)) call TLS_fatal ('error validating mesh_face_set')

  end subroutine test_mesh_face_set

end module mesh_face_set_impl
