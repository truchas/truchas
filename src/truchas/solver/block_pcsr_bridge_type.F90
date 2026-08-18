!!
!! BLOCK_PCSR_BRIDGE_TYPE
!!
!! This module defines BLOCK_PCSR_BRIDGE, which adapts a fixed-block sparse
!! matrix to scalar PCSR storage. It owns the expanded scalar index map needed
!! by the PCSR graph, and caches the scalar matrix structure between updates.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module block_pcsr_bridge_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication, only: is_IOP
  use index_map_type
  use pbsr_matrix_type
  use pcsr_matrix_type
  implicit none
  private

  type, public :: block_pcsr_bridge
    private
    type(index_map), pointer :: scalar_imap => null()
    type(pcsr_graph), pointer :: block_graph => null()  ! unowned reference
    type(pcsr_matrix) :: scalar_matrix
    integer :: block_size = 0
  contains
    procedure :: update
    procedure :: matrix
    final :: delete
  end type

contains

  subroutine delete(this)
    type(block_pcsr_bridge), intent(inout) :: this

    if (associated(this%scalar_imap)) deallocate(this%scalar_imap)
  end subroutine


  !! Update the cached scalar PCSR matrix from A. The PCSR graph is rebuilt
  !! only when A's graph or block size changes.
  subroutine update(this, A)
    class(block_pcsr_bridge), intent(inout) :: this
    type(pbsr_matrix), intent(in) :: A

    integer :: i, ii, ic, xj, j, jj, jc

    ASSERT(associated(A%graph))
    if (.not.associated(this%scalar_imap) .or. .not.associated(this%block_graph, A%graph) .or. &
        this%block_size /= A%bsize) then
      call rebuild(this, A)
    end if

    call this%scalar_matrix%set_all(0.0_r8)
    do i = 1, A%nrow
      do xj = A%graph%xadj(i), A%graph%xadj(i+1)-1
        j = A%graph%adjncy(xj)
        do jj = 1, A%bsize
          jc = A%bsize*(j - 1) + jj
          do ii = 1, A%bsize
            ic = A%bsize*(i - 1) + ii
            call this%scalar_matrix%set(ic, jc, A%values(ii,jj,xj))
          end do
        end do
      end do
    end do
  end subroutine


  subroutine rebuild(this, A)
    class(block_pcsr_bridge), intent(inout) :: this
    type(pbsr_matrix), intent(in) :: A

    type(pcsr_graph), pointer :: graph
    integer, allocatable :: nvars(:)
    integer :: i, ii, ic, xj, j, jj, jc

    if (associated(this%scalar_imap)) deallocate(this%scalar_imap)
    allocate(this%scalar_imap, graph)
    allocate(nvars(merge(A%graph%row_imap%global_size, 0, is_IOP)))
    nvars = A%bsize
    call this%scalar_imap%init(A%graph%row_imap, nvars)

    call graph%init(this%scalar_imap)
    do i = 1, A%nrow
      do xj = A%graph%xadj(i), A%graph%xadj(i+1)-1
        j = A%graph%adjncy(xj)
        do jj = 1, A%bsize
          jc = A%bsize*(j - 1) + jj
          do ii = 1, A%bsize
            ic = A%bsize*(i - 1) + ii
            call graph%add_edge(ic, jc)
          end do
        end do
      end do
    end do
    call graph%add_complete()
    call this%scalar_matrix%init(graph, take_graph=.true.)
    this%block_graph => A%graph
    this%block_size = A%bsize
  end subroutine


  function matrix(this)
    class(block_pcsr_bridge), intent(in), target :: this
    type(pcsr_matrix), pointer :: matrix

    ASSERT(associated(this%scalar_imap))
    matrix => this%scalar_matrix
  end function

end module block_pcsr_bridge_type
