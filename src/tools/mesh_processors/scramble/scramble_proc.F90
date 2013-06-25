!!
!! SCRAMBLE_PROC
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! 2 December 2006
!!
!! This module provides the mesh transformation procedure for
!! the utility program scramble.
!!

#include "f90_assert.fpp"

module scramble_proc

  implicit none
  private

  public :: scramble_elem_nodes

contains

  subroutine scramble_elem_nodes (seed, mesh)

  use string_utilities
  use exodus_mesh_type

  integer, intent(in) :: seed
  type(exodus_mesh), intent(inout) :: mesh

  integer :: n, j, p, nperm, offset
  integer, pointer :: vperm(:,:), fperm(:,:)
  integer, allocatable :: put(:)
  real :: r

  type :: elem_face_map
    integer, pointer :: ptr(:) => null()
  end type
  type(elem_face_map) :: new_faces(mesh%num_elem)

  !! Equivalent relabelings of the cell vertices (new-to-old).
  integer, target, save :: HEX8_VPERM(8,0:23)
  data HEX8_VPERM /1,2,3,4, 5,6,7,8,  1,5,6,2, 4,8,7,3,  1,4,8,5, 2,3,7,6, &
                   2,3,4,1, 6,7,8,5,  2,6,7,3, 1,5,8,4,  2,1,5,6, 3,4,8,7, &
                   3,4,1,2, 7,8,5,6,  3,7,8,4, 2,6,5,1,  3,2,6,7, 4,1,5,8, &
                   4,1,2,3, 8,5,6,7,  4,8,5,1, 3,7,6,2,  4,3,7,8, 1,2,6,5, &
                   5,8,7,6, 1,4,3,2,  5,1,4,8, 6,2,3,7,  5,6,2,1, 8,7,3,4, &
                   6,5,8,7, 2,1,4,3,  6,2,1,5, 7,3,4,8,  6,7,3,2, 5,8,4,1, &
                   7,6,5,8, 3,2,1,4,  7,3,2,6, 8,4,1,5,  7,8,4,3, 6,5,1,2, &
                   8,7,6,5, 4,3,2,1,  8,4,3,7, 5,1,2,6,  8,5,1,4, 7,6,2,3/

  !! Corresponding relabelings of the local face numbers (old-to-new).
  integer, target, save :: HEX8_FPERM(6,0:23)
  data HEX8_FPERM /1,2,3,4,5,6,  5,3,6,1,4,2,  4,6,2,5,1,3, &
                   4,1,2,3,5,6,  1,5,3,6,4,2,  5,4,6,2,1,3, &
                   3,4,1,2,5,6,  6,1,5,3,4,2,  2,5,4,6,1,3, &
                   2,3,4,1,5,6,  3,6,1,5,4,2,  6,2,5,4,1,3, &
                   4,3,2,1,6,5,  1,6,3,5,2,4,  5,2,6,4,3,1, &
                   1,4,3,2,6,5,  5,1,6,3,2,4,  4,5,2,6,3,1, &
                   2,1,4,3,6,5,  3,5,1,6,2,4,  6,4,5,2,3,1, &
                   3,2,1,4,6,5,  6,3,5,1,2,4,  2,6,4,5,3,1/

  !! Set the random number seed.
  call random_seed (size=n)
  allocate(put(n))
  put = seed
  call random_seed (put=put)
  deallocate(put)

  offset = 0
  do n = 1, mesh%num_eblk

    !! Setup permutation data for this block's element type.
    select case (mesh%eblk(n)%elem_type)
    case ('HEX', 'HEX8')
      nperm = size(HEX8_VPERM,dim=2)
      vperm => HEX8_VPERM
      fperm => HEX8_FPERM
    case default
      !! Scrambling other element types not implemented.
      nperm = 0
      write(*,fmt='(a)') 'Not relabeling element nodes for ' // trim(mesh%eblk(n)%elem_type) // &
        'elements in element block ' // i_to_c(mesh%eblk(n)%ID)
    end select

    do j = 1, mesh%eblk(n)%num_elem
      call random_number (r)
      p = int(r*nperm)
      if (p == 0) cycle ! identity permutation
      new_faces(j+offset)%ptr => fperm(:,p)
      mesh%eblk(n)%connect(:,j) = mesh%eblk(n)%connect(vperm(:,p),j)
    end do

    offset = offset + mesh%eblk(n)%num_elem

  end do

  !! Map sideset face data to the new face numbers.
  do n = 1, mesh%num_sset
    do j = 1, mesh%sset(n)%num_side
      if (associated(new_faces(mesh%sset(n)%elem(j))%ptr)) &
        mesh%sset(n)%face(j) = new_faces(mesh%sset(n)%elem(j))%ptr(mesh%sset(n)%face(j))
    end do
  end do

  end subroutine scramble_elem_nodes

end module scramble_proc
