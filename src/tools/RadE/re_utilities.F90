!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module re_utilities

  use scl
  implicit none
  private

  public :: re_info, re_halt
  public :: compute_face_area

contains

  subroutine re_info (mesg)
    character(len=*), intent(in) :: mesg
    if (scl_rank()==1) write(*,'(a)') mesg
  end subroutine re_info

  subroutine re_halt (errmsg)
    character(len=*), intent(in) :: errmsg
    if (scl_rank()==1) write(*,'(2a)') 'FATAL: ', errmsg
    call scl_finalize ()
    stop 1
  end subroutine re_halt

  !! Returns the area of enclosure faces on rank 1.
  subroutine compute_face_area(xface, fnode, x, area)

    use kinds, only: r8
    use cell_geometry, only: face_normal, vector_length

    integer, allocatable, intent(in) :: xface(:), fnode(:)
    real(r8), allocatable, intent(in) :: x(:,:)
    real(r8), allocatable, intent(out) :: area(:)

    integer :: i, n, nface

    nface = size(xface)-1
    n = merge(nface, 0, scl_rank()==1)
    allocate(area(n))

    if (scl_rank() == 1) then
      do i = 1, nface
        associate (face_nodes => fnode(xface(i):xface(i+1)-1))
          area(i) = vector_length(face_normal(x(:,face_nodes)))
        end associate
      end do
    end if

  end subroutine compute_face_area

end module re_utilities
